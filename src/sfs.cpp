#include <getopt.h>
#include <string>
#include <algorithm>
#include <vector>
#include <iterator>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <unordered_map>
#include <set>
#include <cstdlib>
#include <utility>
#include <cmath>
#include <random>
#include <math.h>
#include <zlib.h>
#include <htslib/vcf.h>

using std::cout;
using std::endl;
using namespace std;

/**
 * Print a help message to the terminal and exit.
 */
void help(int code){
    fprintf(stderr, "sfs [OPTIONS]\n");
    fprintf(stderr, "Given a VCF and a file mapping individuals in the VCF to populations\n");
    fprintf(stderr, "Computes site frequency spectrum (SFS)-based statistics in each population.\n");
    fprintf(stderr, "Statistics include:\n");
    fprintf(stderr, "Tajima's D (Tajima 1989; 10.1093/genetics/123.3.585)\n");
    fprintf(stderr, "  compares medium- to low-frequency segregating sites.\n");
    fprintf(stderr, "Fay and Wu's H (Fay and Wu 2000; 10.1093/genetics/155.3.1405)\n");
    fprintf(stderr, "  compares medium- to high-frequency segregating sites.\n");
    fprintf(stderr, "Zeng's E (Zeng et al 2006; 10.1534/genetics.106.061432)\n");
    fprintf(stderr, "  compares low- to high-frequency segregating sites.\n");
    fprintf(stderr, "H and E are computed and scaled as recommended in 10.1534/genetics.106.061432.\n");
    fprintf(stderr, "Normalized versions of statistics are also provided:\n");
    fprintf(stderr, "D_Dmin is D divided by its minimum possible value, as proposed\n");
    fprintf(stderr, "  in Schaeffer (2002) (10.1017/s0016672302005955)\n");
    fprintf(stderr, "H_Hmin and E_Emin are computed in a similar way.\n");
    fprintf(stderr, "Dprime, Hprime, and Eprime are scaled by dividing negative values by\n");
    fprintf(stderr, "  the theoretical minimum value possible, and dividing positive\n");
    fprintf(stderr, "  values by the theoretical maximum possible value.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "If windows (BED format) are provided, then Tajima's D will be computed in\n");
    fprintf(stderr, "each window.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "The population file should be two-column, tab separated, where the first\n");
    fprintf(stderr, "column is individual names from the VCF, and the second is arbitrary text\n");
    fprintf(stderr, "population labels.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "===== OPTIONS =====\n");
    fprintf(stderr, "--vcf -v VCF/BCF file. REQUIRED.\n");
    fprintf(stderr, "--pops -p Population file (described above). REQUIRED.\n");
    fprintf(stderr, "--bed -b Compute Tajima's D in each BED interval (- for stdin)\n");
    fprintf(stderr, "--outgroup -o Population ID to use for outgroup to determine ancestral/derived\n");
    fprintf(stderr, "  alleles. This is required for computing Fay & Wu's H and Zeng's E, unless\n");
    fprintf(stderr, "  the ancestral allele (AA) field is populated in the VCF.\n");
    fprintf(stderr, "--ref_ancestral -r Assume the reference allele is always ancestral. Removes\n");
    fprintf(stderr, "  the need to supply --outgroup for Fay & Wu's H and Zeng's E.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "--help -h Display this message and exit.\n");
    exit(code);
}

void read_pops(string& popsfile, 
    map<string, vector<string> >& pops){
    ifstream inf(popsfile.c_str());
    string indvname;
    string popname;
    set<string> popuniq;
    while (inf >> indvname >> popname){
        if (pops.count(popname) == 0){
            vector<string> v;
            pops.insert(make_pair(popname, v));
        }
        pops[popname].push_back(indvname);
        popuniq.insert(popname);
    }
    if (popuniq.size() < 2){
        fprintf(stderr, "ERROR: fewer than two populations read from %s\n", popsfile.c_str());
        exit(1);
    }
}

void read_bed(string& bedfile, map<string, vector<pair<int, int> > >& bed){
    ifstream inf(bedfile);
    if (!inf.is_open()){
        fprintf(stderr, "ERROR opening file %s\n", bedfile.c_str());
        exit(1);
    }
    string line;
    while (getline(inf, line)){
        istringstream splitter(line);
        string field;
        int idx = 0;
        string chrom;
        int start;
        int end;
        while(getline(splitter, field, '\t')){
            if (idx == 0){
                chrom = field;
                if (bed.count(chrom) == 0){
                    vector<pair<int, int> > v;
                    bed.insert(make_pair(chrom, v));
                }
            }
            else if (idx == 1){
                start = atoi(field.c_str());
            }
            else if (idx == 2){
                end = atoi(field.c_str());
                bed[chrom].push_back(make_pair(start, end));
            }
            else{
                break;
            }
            ++idx;
        }
    }
    // Make sure everything is sorted.
    for (map<string, vector<pair<int, int> > >::iterator b = bed.begin(); b != bed.end(); ++b){
        sort(b->second.begin(), b->second.end());
    }
}

struct tajD_dat{
    double pi_sum;
    double pi_sum_min;
    double pi_sum_max;
    double pi_comps;
    double segsites;
    double totsites;
    double hsum;
    double hsum_min;
    double hsum_max;
    double n;

    tajD_dat(){
        pi_sum = 0;
        pi_sum_min = 0;
        pi_sum_max = 0;
        pi_comps = 0;
        segsites = 0;
        totsites = 0;
        hsum = 0;
        hsum_min = 0;
        hsum_max = 0;
        n = 0;
    }
};

void proc_record(map<string, tajD_dat>& dat,
    map<string, vector<int> >& pop_idx,
    bcf_hdr_t* bcf_header,
    bcf1_t* bcf_record,
    int num_samples,
    vector<int>& og_idx,
    bool ref_anc,
    bool has_AA){
    
    // Only consider biallelic SNPs
    if (bcf_record->n_allele == 2){
        bcf_unpack(bcf_record, BCF_UN_STR);
        // pass = biallelic, no indels, ref/alt both A/C/G/T
        bool pass = true;
        for (int i = 0; i < bcf_record->n_allele; ++i){
            if (strcmp(bcf_record->d.allele[i], "A") != 0 &&
                strcmp(bcf_record->d.allele[i], "C") != 0 &&
                strcmp(bcf_record->d.allele[i], "G") != 0 && 
                strcmp(bcf_record->d.allele[i], "T") != 0){
                pass = false;
                break;
            }
        }
        if (bcf_record->d.allele[0][0] == bcf_record->d.allele[1][0]){
            pass = false;
        }
        if (pass){

            // Get all available genotypes.
            int32_t* gts = NULL;
            int n_gts = 0;
            int nmiss = 0;
            int num_loaded = bcf_get_genotypes(bcf_header, bcf_record, &gts, &n_gts);
            if (num_loaded <= 0){
                fprintf(stderr, "ERROR loading genotypes at %ld %ld\n", 
                    (long int) bcf_record->rid, (long int) bcf_record->pos);
                exit(1);
            }
            int ploidy = n_gts / num_samples; 
            
            // Determine ancestral allele, if possible.
            int anc_allele = -1;
            if (ref_anc){
                anc_allele = 0;
            }
            else if (has_AA){
                // Look for AA tag.
                bcf_unpack(bcf_record, BCF_UN_INFO);
                char* aa_str = NULL;
                int aa_str_len = 0;
                int success = bcf_get_info_string(bcf_header, bcf_record, "AA", &aa_str, &aa_str_len);
                if (success && aa_str_len > 0 && aa_str != NULL){
                    for (int x = 0; x < aa_str_len; ++x){
                        if (aa_str[x] != 'A' && aa_str[x] != 'C' && aa_str[x] != 'G' && aa_str[x] != 'T'){
                            aa_str[x] = '\0';
                            aa_str_len = x;
                            break;
                        }
                    }
                    for (int aidx = 0; aidx < bcf_record->n_allele; ++aidx){
                        if (strcmp(bcf_record->d.allele[aidx], aa_str) == 0){
                            anc_allele = aidx;
                            break;
                        }
                    }
                }
                free(aa_str);
                if (anc_allele == -1){
                    //return;
                }
            }
            else if (og_idx.size() > 0){
                int c0 = 0;
                int c1 = 0;
                for (vector<int>::iterator a = og_idx.begin(); a != og_idx.end(); ++a){
                    int32_t* ptr = gts + (*a)*ploidy;
                    if (!bcf_gt_is_missing(ptr[0])){
                        for (int pi = 0; pi < ploidy; ++pi){
                            if (bcf_gt_allele(ptr[pi]) == 0){
                                c0++;
                            }
                            else{
                                c1++;
                            }
                        }
                    }
                }
                if (c0 > 0 && c1 == 0){
                    anc_allele = 0;
                }
                else if (c0 == 0 && c1 > 0){
                    anc_allele = 1;
                }
                else{
                    //return;
                }
            }
            
            for (map<string, vector<int> >::iterator p = pop_idx.begin(); 
                p != pop_idx.end(); ++p){
                
                int popcomps = 0;
                int popdiffs = 0;
                
                // Track derived allele freqs
                int daf = 0;

                int n_eff = 0;
                
                int n = 0;

                for (int i = 0; i < p->second.size(); ++i){
                    int i_x = p->second[i];
                    int32_t* ptr1 = gts + (i_x)*ploidy;
                    if (!bcf_gt_is_missing(ptr1[0])){
                        n += 2;
                        n_eff += ploidy;
                        // Account for self-self
                        for (int pi = 0; pi < ploidy; ++pi){
                            if (anc_allele != -1 && bcf_gt_allele(ptr1[pi]) != anc_allele){
                                daf++;
                            }
                            for (int pj = pi + 1; pj < ploidy; ++pj){
                                if (bcf_gt_allele(ptr1[pi]) != bcf_gt_allele(ptr1[pj])){
                                    popdiffs++;
                                }
                                popcomps++;
                            }
                        }
                        for (int j = i + 1; j < p->second.size(); ++j){
                            int j_x = p->second[j];
                            int32_t* ptr2 = gts + (j_x)*ploidy;
                            if (!bcf_gt_is_missing(ptr2[0])){
                                for (int pi = 0; pi < ploidy; ++pi){
                                    for (int pj = 0; pj < ploidy; ++pj){
                                        if (bcf_gt_allele(ptr1[pi]) != bcf_gt_allele(ptr2[pj])){
                                            popdiffs++;
                                        }
                                        popcomps++;
                                    }
                                }
                            }
                        }
                    }
                }
                
                if (dat.count(p->first) == 0){
                    dat.insert(make_pair(p->first, tajD_dat()));
                    dat[p->first].n = (double)p->second.size();
                }
                
                if (popdiffs > 0){
                    double nchoosek = (double)n*(double)(n-1) / 2.0;
                    dat[p->first].segsites++;
                    dat[p->first].pi_sum += (double)popdiffs / (double)popcomps;
                    dat[p->first].pi_sum_min += (1.0*((double)n - 1.0))/nchoosek;
                    dat[p->first].pi_sum_max += (pow((double)n/2.0, 2))/nchoosek;
                    //dat[p->first].pi_sum += popdiffs;
                    //dat[p->first].pi_comps += popcomps;
                }
                dat[p->first].totsites++;

                if (anc_allele != -1 && daf > 1 && daf < n){
                    // Incorporate derived allele freq data
                    dat[p->first].hsum += (1.0/(double)(n_eff-1)) * (double)daf;
                    dat[p->first].hsum_min += (1.0/(double)(n_eff-1)) * (double)1.0;
                    dat[p->first].hsum_max += (1.0/(double)(n_eff-1)) * (double)(n-1);
                }
            }
        }
    }
}

void print_record(const string& chrom,
    const pair<int, int>& interval,
    const string& pop,
    tajD_dat& dat,
    double a1,
    double e1,
    double e2,
    double hvar_fac1,
    double hvar_fac2,
    double bn,
    bool has_anc){

    if (dat.totsites == 0){
        return;
    }
    if (dat.segsites == 0){
        return;
    }

    if (interval.first == -1 && interval.second == -1){
        // No interval.
    }
    else{
        fprintf(stdout, "%s\t%d\t%d\t", chrom.c_str(), interval.first, interval.second);
    }
    if (pop != "all"){
        fprintf(stdout, "%s\t", pop.c_str());
    }

    if (dat.pi_comps == 0){
        dat.pi_comps = 1.0;
    }
    if (dat.totsites == 0){
        dat.totsites = 1.0;
    }

    // Watterson's estimator = low frequency (every site counts the same)
    // pi = mid frequency (mean variation)
    // theta_L = high frequency (weight sites by frequency)
    
    double nchoose2_inv = 2.0/(dat.n*(dat.n-1.0));
    //double pi = dat.pi_sum * nchoose2_inv;
    double pi = dat.pi_sum;
    //double pi = dat.pi_sum / dat.pi_comps;
    double wat = dat.segsites / a1;
    // Minimum pi happens when each seg site is MAF == 1
    double kmin = (double)dat.segsites * (1.0 - dat.n - 1.0) * nchoose2_inv;
    // Maximum pi happens when each seg site is MAF == 0.5
    double kmax = ((double)dat.segsites * ((dat.n/2.0)*(dat.n/2.0))) * nchoose2_inv;
    
    kmin = dat.pi_sum_min;
    kmax = dat.pi_sum_max;

    double v_D = sqrt(e1*dat.segsites + e2*dat.segsites*(dat.segsites - 1));
    double D = (pi - wat)/v_D;
    double Dmin = kmin - wat;
    double Dmax = kmax - wat;
    double Dprime1 = (pi - wat) / abs(Dmin);
    double Dprime = 0;
    if (D < 0){
        Dprime = (pi - wat) / abs(Dmin);
    }
    else if (D > 0){
        Dprime = (pi - wat) / abs(Dmax);
    }
    //double Dprime = (pi - wat) / abs(Dmin);
     
    fprintf(stdout, "%f\t%f\t%d\t%f\t%f\t%f", pi, wat, (int)round(dat.totsites),
        D, Dprime1, Dprime);

    if (has_anc){
        if (dat.hsum == 0.0){
            // No polarizable sites in this window
            fprintf(stdout, "\tNA\tNA\tNA\tNA\tNA\tNA\tNA");
        }
        else{
            // Normalized Fay & Wu's H
            double theta_L = dat.hsum;
            double t2 = ((double)dat.segsites*((double)dat.segsites - 1.0)) / (a1*a1 + bn);
            double h_var = hvar_fac1 * wat + hvar_fac2 * t2;
            double H = (pi - theta_L)/sqrt(h_var);
            // Normalize again in the spirit of D'
            double Hprime = 0;
            double Hmin = kmin - dat.hsum_max;
            double Hmax = kmax - dat.hsum_min;
            double Hprime1 = (pi - theta_L) / abs(Hmin);
            if (H < 0){
                Hprime = (pi - theta_L) / abs(Hmin);
            }
            else if (H > 0){
                Hprime = (pi - theta_L) / abs(Hmax);
            }
            //double Hprime = (pi - theta_L) / abs(kmin - theta_L);
            // Zeng's E
            double e_var = (dat.n/(2*(dat.n-1)) - 1.0/a1) * wat;
            double fac2 = bn/(a1*a1) + 2.0*pow(dat.n/(dat.n-1),2)*bn;
            fac2 -= (1*(dat.n*bn - dat.n + 1.0))/((dat.n-1.0)*a1);
            fac2 -= (3*dat.n + 1.0)/(dat.n - 1.0);
            e_var += fac2*t2;
            double E = (theta_L - wat)/sqrt(e_var);
            double Emin = (dat.hsum_min - wat);
            double Emax = (dat.hsum_max - wat);
            double Eprime1 = (theta_L - wat)/abs(Emin);
            double Eprime = 0.0;
            if (E < 0){
                Eprime = (theta_L - wat)/abs(Emin);
            }
            else if (E > 0){
                Eprime = (theta_L - wat)/abs(Emax);
            }
            fprintf(stdout, "\t%f\t%f\t%f\t%f\t%f\t%f\t%f", theta_L, H, Hprime1, Hprime, E, Eprime1, Eprime);
        } 
    }
    fprintf(stdout, "\n");
}

void read_vcf(htsFile* bcf_reader, 
    bcf_hdr_t* bcf_header,
    bcf1_t* bcf_record, 
    int num_samples,
    map<string, vector<int> >& pop_idx,
    map<string, vector<pair<int, int> > >& bed,
    map<string, double>& a1,
    map<string, double>& e1,
    map<string, double>& e2,
    map<string, double>& hvar_fac1,
    map<string, double>& hvar_fac2,
    map<string, double>& bn,
    bool ref_ancestral,
    vector<int>& og_idx,
    bool AA){
    
    bool has_anc = AA || ref_ancestral || og_idx.size() > 0;

    map<pair<int, int>, map<string, tajD_dat> > cur_dat;

    long int progress = 1000;
    long int nsnp = 0;
    
    int prevrid = -1;
    int prevbin = -1;

    vector<string> pop_names;
    map<string, int> pop_name2idx;

    for (map<string, vector<int> >::iterator pi = pop_idx.begin(); pi != 
        pop_idx.end(); ++pi){
        pop_name2idx.insert(make_pair(pi->first, pop_names.size()));
        pop_names.push_back(pi->first);
    }
    
    vector<pair<int, int> >::iterator cur_bed;
    string cur_chrom = "";
    bool has_chrom = false;
    
    pair<int, int> nullkey = make_pair(-1, -1);

    while(bcf_read(bcf_reader, bcf_header, bcf_record) == 0){
        
        bool use_snp = false;
        if (bed.size() == 0){
            if (cur_dat.count(nullkey) == 0){
                map<string, tajD_dat> m;
                cur_dat.insert(make_pair(nullkey, m));
            }
            proc_record(cur_dat[nullkey], pop_idx, bcf_header, bcf_record, num_samples,
                og_idx, ref_ancestral, AA);
        }
        else{
            string this_chrom = bcf_hdr_id2name(bcf_header, bcf_record->rid);
            if (this_chrom != cur_chrom){
                // Print data for prev chrom intervals
                if (cur_dat.size() > 0){
                    for (map<pair<int, int>, map<string, tajD_dat> >::iterator ival = cur_dat.begin();
                        ival != cur_dat.end(); ++ival){
                        for (map<string, tajD_dat>::iterator ival2 = ival->second.begin(); 
                            ival2 != ival->second.end(); ++ival2){
                            
                            print_record(cur_chrom, ival->first, ival2->first, ival2->second,
                                a1[ival2->first], e1[ival2->first], e2[ival2->first], 
                                hvar_fac1[ival2->first], hvar_fac2[ival2->first], 
                                bn[ival2->first], has_anc);
                        }
                    }
                    cur_dat.clear();
                }
                if (bed.count(this_chrom) > 0){
                    has_chrom = true;
                    cur_bed = bed[this_chrom].begin();
                }   
                else{
                    has_chrom = false;
                }         
                cur_chrom = this_chrom;
            }
            if (has_chrom){
                while (cur_bed != bed[cur_chrom].end() && cur_bed->second-1 < bcf_record->pos){
                    // Print out data from this interval
                    for (map<string, tajD_dat>::iterator ival = cur_dat[*cur_bed].begin(); 
                        ival != cur_dat[*cur_bed].end(); ++ival){
                        print_record(cur_chrom, *cur_bed, ival->first, ival->second,
                            a1[ival->first], e1[ival->first], e2[ival->first], 
                            hvar_fac1[ival->first], hvar_fac2[ival->first], 
                            bn[ival->first], has_anc);
                    }
                    cur_dat.erase(*cur_bed);
                    cur_bed = bed[cur_chrom].erase(cur_bed);
                }
                if (cur_bed != bed[cur_chrom].end() && cur_bed->first <= bcf_record->pos && 
                    cur_bed->second > bcf_record->pos){
                    vector<pair<int, int> >::iterator cur_bed2 = cur_bed;
                    while (cur_bed2 != bed[cur_chrom].end() && cur_bed2->first <= bcf_record->pos){
                        if (cur_bed2->first <= bcf_record->pos && cur_bed2->second > bcf_record->pos){
                            // Create entry if it doesn't exist.
                            if (cur_dat.count(*cur_bed2) == 0){
                                map<string, tajD_dat> m;
                                cur_dat.insert(make_pair(*cur_bed2, m));
                            }
                            proc_record(cur_dat[*cur_bed2], pop_idx, bcf_header, bcf_record, num_samples,
                                og_idx, ref_ancestral, AA);
                        }
                        ++cur_bed2;
                    }
                }       
            }
        }
        nsnp++;
        if (nsnp % progress == 0){
            fprintf(stderr, "Read %ld snps\r", nsnp);
        }
    }
    fprintf(stderr, "Read %ld snps\n", nsnp);
    
    // Print final window data
    for (map<pair<int, int>, map<string, tajD_dat> >::iterator dat = cur_dat.begin(); 
        dat != cur_dat.end(); ++dat){
        for (map<string, tajD_dat>::iterator dat2 = dat->second.begin(); dat2 != 
            dat->second.end(); ++dat2){
            print_record(cur_chrom, dat->first, dat2->first, dat2->second,
                a1[dat2->first], e1[dat2->first], e2[dat2->first], 
                hvar_fac1[dat2->first], hvar_fac2[dat2->first], 
                bn[dat2->first], has_anc);
        }
    }
}

/**
 * Checks ahead in a VCF file to see if the ancestral allele (AA) field is 
 * populated in the first n records. Returns true if so and false if not.
 */
bool check_AA(string& filename, int n_record){
    bcf_hdr_t* bcf_header;
    bcf1_t* bcf_record = bcf_init();
    htsFile* bcf_reader = bcf_open(filename.c_str(), "r");
    if (bcf_reader == NULL){
        fprintf(stderr, "ERROR interpreting %s as BCF format.\n", filename.c_str());
        exit(1);
    }
    bcf_header = bcf_hdr_read(bcf_reader);
    int num_samples = bcf_hdr_nsamples(bcf_header);
    
    int nrec = 0;
    while(bcf_read(bcf_reader, bcf_header, bcf_record) == 0){
        // Look for AA tag.
        bcf_unpack(bcf_record, BCF_UN_INFO);
        char* aa_str = NULL;
        int aa_str_len = 0;
        int success = bcf_get_info_string(bcf_header, bcf_record, "AA", &aa_str, &aa_str_len);
        if (success && aa_str_len > 0 && aa_str != NULL){
            free(aa_str);
            return true;
        }
        ++nrec;
        if (nrec >= n_record){
            free(aa_str);
            return false;
        }
        free(aa_str);
    }
    return false;
}


int main(int argc, char *argv[]) {    
    
    static struct option long_options[] = {
       {"vcf", required_argument, 0, 'v'},
       {"pops", required_argument, 0, 'p'},
       {"outgroup", required_argument, 0, 'o'},
       {"bed", required_argument, 0, 'b'},
       {"ref_ancestral", no_argument, 0, 'r'},
       {0, 0, 0, 0} 
    };
    
    // Set default values
    string vcf_file = "";
    string popsfile = "";
    string bedfile = "";
    string outgroup = "";
    bool ref_ancestral = false;

    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "v:p:b:o:rh", long_options, &option_index )) != -1){
        switch(ch){
            case 0:
                // This option set a flag. No need to do anything here.
                break;
            case 'h':
                help(0);
                break;
            case 'v':
                vcf_file = optarg;
                break;
            case 'p':
                popsfile = optarg;
                break;
            case 'b':
                bedfile = optarg;
                break;
            case 'o':
                outgroup = optarg;
                break;
            case 'r':
                ref_ancestral = true;
                break;
            default:
                help(0);
                break;
        }    
    }
    
    // Error check arguments.
    if (vcf_file == ""){
        fprintf(stderr, "ERROR: VCF/BCF file required\n");
        exit(1);
    }
     
    // Read in populations
    map<string, vector<string> > pops;
    if (popsfile != ""){
        read_pops(popsfile, pops);
        if (outgroup != "" && pops.count(outgroup) == 0){
            fprintf(stderr, "ERROR: outgroup population %s not found in population file.\n",
                outgroup.c_str());
            exit(1);
        }
    }
    else if (outgroup != ""){
        fprintf(stderr, "ERROR: outgroup provided without population file\n");
        exit(1);
    }
    if (ref_ancestral && outgroup != ""){
        fprintf(stderr, "ERROR: no need to provide ougroup if treating ref as ancestral alleles.\n");
        exit(1);
    }
   
    map<string, vector<pair<int, int> > > bed;
    if (bedfile != ""){
        read_bed(bedfile, bed);
    } 

    // Init BCF/VCF reader and take in header
    bcf_hdr_t* bcf_header;
    bcf1_t* bcf_record = bcf_init();
    htsFile* bcf_reader = bcf_open(vcf_file.c_str(), "r");
    if (bcf_reader == NULL){
        fprintf(stderr, "ERROR interpreting %s as BCF format.\n", vcf_file.c_str());
        exit(1);
    }
    bcf_header = bcf_hdr_read(bcf_reader);
    int num_samples = bcf_hdr_nsamples(bcf_header);
    map<string, int> samp2idx;
    vector<string> idx2samp;
    if (popsfile == ""){
        vector<string> v;
        pops.insert(make_pair("all", v));
    }
    for (int i = 0; i < num_samples; ++i){
        samp2idx.insert(make_pair(bcf_header->samples[i], i));
        idx2samp.push_back(bcf_header->samples[i]);
        if (popsfile == ""){
            pops["all"].push_back(bcf_header->samples[i]);
        }
    }
    
    // Set up population indices
    map<string, vector<int> > pop_idx;
    // Also pre-compute vals needed for normalization
    map<string, double> pop_a1;
    map<string, double> pop_e1;
    map<string, double> pop_e2;
    map<string, double> pop_hvar_fac1;
    map<string, double> pop_hvar_fac2;
    map<string, double> pop_bn; 

    for (map<string, vector<string> >::iterator p = pops.begin(); p != pops.end(); ++p){
        vector<int> v;
        pop_idx.insert(make_pair(p->first, v));
        for (vector<string>::iterator pi = p->second.begin();
            pi != p->second.end(); ++pi){
            if (samp2idx.count(*pi) == 0){
                fprintf(stderr, "ERROR: sample %s not found in VCF.\n", pi->c_str());
                exit(1);
            }
            pop_idx[p->first].push_back(samp2idx[*pi]);
        }
        double a1 = 0.0;
        double a2 = 0.0;
        int n = p->second.size()*2;
        for (int i = 1; i <= n - 1; ++i){
            a1 += 1.0/(double)(i);
            a2 += 1.0/pow((double)(i), 2);
        }
        double b1 = ((double)n+1)/(3.0*(double)n + 1);
        double b2 = (2*((double)n*(double)n + (double)n + 3))/(9*(double)n*((double)n-1));
        double c1 = b1 - 1.0/a1;
        double c2 = b2 - (double)(n+2)/(a1*(double)n) + a2/(a1*a1);
        double e1 = c1/a1;
        double e2 = c2/(a1*a1 + a2);
        
        pop_a1.insert(make_pair(p->first, a1));
        pop_e1.insert(make_pair(p->first, e1));
        pop_e2.insert(make_pair(p->first, e2));
        
        pop_hvar_fac1.insert(make_pair(p->first, (double)(n-2)/(double)(6*(n-1))));
        double bn_1 = 0.0;
        double bn = 0.0;
        for (int x = 1; x <= n; ++x){
            bn_1 += (1.0/(double)(x*x));
            if (x < n){
                bn += (1.0/(double)(x*x));
            }
        }
        double n_d = (double)n;
        double hvar_fac2 = (18.0*n_d*n_d*(3*n_d + 2.0)*bn_1 - (88.0*n_d*n_d*n_d + 9.0*n_d*n_d - 
            13.0*n_d + 6.0)) / (9*n_d*(n_d-1.0)*(n_d-1.0));
        pop_hvar_fac2.insert(make_pair(p->first, hvar_fac2));
        pop_bn.insert(make_pair(p->first, bn));
        
    }
   
    vector<int> og_idx;

    // Check for ancestral alleles
    bool has_anc = false;
    bool has_AA = false;
    if (outgroup == "" && !ref_ancestral){
        has_anc = check_AA(vcf_file, 5000);
        if (has_anc){
            has_AA = true;
        }
    }
    else{
        has_anc = true;
        if (outgroup != ""){
            og_idx = pop_idx[outgroup];
            pop_idx.erase(outgroup);
        }
    }

    if (bed.size() > 0){
        fprintf(stdout, "chrom\tstart\tend\t");
    } 
    if (popsfile != ""){
        fprintf(stdout, "pop\t");
    }
    fprintf(stdout, "theta_pi\ttheta_w\tnsites\tD\tD_dmin\tDprime");
    if (has_anc){
        fprintf(stdout, "\ttheta_H\tH\tH_hmin\tHprime\tE\tE_Emin\tEprime");
    }
    fprintf(stdout, "\n");
    
    // Ready to go
    read_vcf(bcf_reader, bcf_header, bcf_record, num_samples, pop_idx,
        bed, pop_a1, pop_e1, pop_e2, pop_hvar_fac1, pop_hvar_fac2,
        pop_bn, ref_ancestral, og_idx, has_AA);
}

