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
    fprintf(stderr, "dstat [OPTIONS]\n");
    fprintf(stderr, "Given a VCF and a file mapping individuals in the VCF to populations\n");
    fprintf(stderr, "(candidate admixed, candidate admixer, candidate un-admixed, and outgroup)\n");
    fprintf(stderr, "plus a window size (in bp), computes the D statistic for each combination\n");
    fprintf(stderr, "of individuals, along with weighted block jackknife standard errors and\n");
    fprintf(stderr, "p-values.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "The population file should be 2 columns, tab separated. The first column\n");
    fprintf(stderr, "is the name of an individual in the VCF and the second column is a population\n");
    fprintf(stderr, "identifier:\n");
    fprintf(stderr, "1 = Population 1 (in Green et al (2010), this would be Yoruba)\n");
    fprintf(stderr, "2 = Population 2 (in Green et al (2010), this would be French)\n");
    fprintf(stderr, "T = Test population (in Green et al (2010), this would be Neanderthal)\n");
    fprintf(stderr, "O = outgroup (in Green et al (2010), this would be chimpanzee)\n");
    fprintf(stderr, "D stats will be positive when T individuals are closer to 2, or negative\n");
    fprintf(stderr, "when T individuals are closer to 1.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "You can omit either 1 or 2, as long as the included population (1 or 2) has\n");
    fprintf(stderr, "more than one individual. In this case, each individual from the included\n");
    fprintf(stderr, "population will be tested as a representative of the other population.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "===== OPTIONS =====\n");
    fprintf(stderr, "--vcf -v VCF/BCF file. REQUIRED.\n");
    fprintf(stderr, "--pops -p Population file (described above). REQUIRED.\n");
    fprintf(stderr, "--window -w Window size for weighted block jackknife in Mb (default = 10)\n");
    fprintf(stderr, "  Set to 0 to disable calculating standard errors.\n");
    fprintf(stderr, "--f -f Also estimate admixture proportions using the f-hat statistic. This\n");
    fprintf(stderr, "  is computed, for each test individual, as D(1,2,T,O)/D(1,T1,T2,O) where T1\n");
    fprintf(stderr, "  and T2 are two individuals from the test (T) population. Requires at least\n");
    fprintf(stderr, "  two T individuals.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "--help -h Display this message and exit.\n");
    exit(code);
}
struct dcomp{
    int idx1;
    int idx2;
    int idxT;
    bool for_f;

    dcomp(){
        idx1 = -1;
        idx2 = -1;
        idxT = -1;
        for_f = false;
    }

    dcomp(int one, int two, int t){
        idx1 = one;
        idx2 = two;
        idxT = t;
        for_f = false;
    }
};

struct wbj_counts{
    double abba;
    double baba;
    int informative_sites;
    int chrom;
    int bin;

    wbj_counts(){
        abba = 0.0;
        baba = 0.0;
        informative_sites = 0;
        chrom = -1;
        bin = -1;
    }
};

void read_pops(string& popsfile, 
    vector<string>& pop1,
    vector<string>& pop2,
    vector<string>& popT,
    vector<string>& popO){
    ifstream inf(popsfile.c_str());
    string indvname;
    string popname;
    while (inf >> indvname >> popname){
        if (popname == "1"){
            pop1.push_back(indvname);
        }
        else if (popname == "2"){
            pop2.push_back(indvname);
        }
        else if (popname == "T" || popname == "t"){
            popT.push_back(indvname);
        }
        else if (popname == "O" || popname == "o"){
            popO.push_back(indvname);
        }
        else{
            fprintf(stderr, "ERROR: unknown population ID %s\n", popname.c_str());
            fprintf(stderr, "Pops should be 1,2,T, or O, where\n");
            fprintf(stderr, " 1 (P1) = Population 1\n");
            fprintf(stderr, " 2 (P2) = Population 2\n");
            fprintf(stderr, " T (P3) = Test - are these closer to pop 1 or pop 2?\n");
            fprintf(stderr, " O (P4) = Outgroup\n");
            fprintf(stderr, "\n");
            fprintf(stderr, " e.g. as used in Green (2010), pops would be:\n");
            fprintf(stderr, " 1 = Yoruba, 2 = French, T = Neanderthal, O = chimpanzee\n");
            exit(1);
        }
    }
    if (popO.size() == 0){
        fprintf(stderr, "ERROR: no outgroup individuals (pop = O) present.\n");
        exit(1);
    }
    if (popT.size() == 0){
        fprintf(stderr, "ERROR: no test individuals (pop = T) present.\n");
        exit(1);
    }
    if (pop1.size() == 0 && pop2.size() == 0){
        fprintf(stderr, "ERROR: no individuals from population 1 or 2 present.\n");
        exit(1);
    }
    else if (pop1.size() == 0){
        if (pop2.size() == 1){
            fprintf(stderr, "ERROR: no pop 1 individuals and only one pop 2 individual.\n");
            exit(1);
        }
        fprintf(stderr, "No pop 1 individuals. Treating pop 2 as both pop 1 and pop 2.\n"); 
    }
    else if (pop2.size() == 0){
        if (pop1.size() == 1){
            fprintf(stderr, "ERROR: no pop 2 individuals and only one pop 1 individual.\n");
            exit(1);
        }
        fprintf(stderr, "No pop 2 individuals. Treating pop 1 as both pop 1 and pop 2.\n");
    }
}

void read_vcf(htsFile* bcf_reader, 
    bcf_hdr_t* bcf_header,
    bcf1_t* bcf_record, 
    int num_samples,
    vector<dcomp>& dcomps,
    vector<int>& dcomp_inf_sites,
    set<int>& ingroup_all,
    vector<int>& outgroup_idx,
    vector<pair<double, double> >& abcounts,
    map<int, vector<wbj_counts> >& wbj_dat,
    int winsize){
    
    long int progress = 1000;
    long int nsnp = 0;
    
    int binprev = -1;
    int ridprev = -1;

    vector<double> derprobs(num_samples, -1.0);
    
    while(bcf_read(bcf_reader, bcf_header, bcf_record) == 0){
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
                //int ploidy = 2;
                int ploidy = n_gts / num_samples; 
                // Get current WBJ bin
                int binstart = (int)round(((double)bcf_record->pos / 
                    (double)winsize)) * (double)winsize;
                int chrom_idx = bcf_record->rid;
                bool new_bin = false;
                if (chrom_idx != ridprev || binstart != binprev){
                    new_bin = true;
                }
                ridprev = chrom_idx;
                binprev = binstart;
                 
                // Get ancestral allele, skip if not fixed in outgroup
                int anc_allele = -1;
                bool ogfixed = true;
                for (int i = 0; i < outgroup_idx.size(); ++i){
                    int32_t* optr = gts + outgroup_idx[i]*ploidy;
                    if (!bcf_gt_is_missing(optr[0])){
                        for (int pi = 0; pi < ploidy; ++pi){
                            if (anc_allele == -1){
                                anc_allele = bcf_gt_allele(optr[pi]);
                            }
                            else if (anc_allele != bcf_gt_allele(optr[pi])){
                                ogfixed = false;
                                break;
                            }
                        }
                    }
                }
                if (ogfixed){
                    // Get derived allele probs of all relevant individuals
                    //for (set<int>::iterator ig = ingroup_all.begin();
                    //    ig != ingroup_all.end(); ++ig){
                    //    int32_t* ptr = gts + (*ig)*ploidy;
                    for (int n = 0; n < num_samples; ++n){
                        int32_t* ptr = gts + n*ploidy;
                        if (!bcf_gt_is_missing(ptr[0])){
                            double gt = 0.0;
                            for (int pi = 0; pi < ploidy; ++pi){
                                if (bcf_gt_allele(ptr[pi]) != anc_allele){
                                    gt += 1.0/(double)ploidy;
                                }
                            }
                            //derprobs[*ig] = gt;
                            derprobs[n] = gt;
                        }
                        else{
                            //derprobs[*ig] = -1.0;
                            derprobs[n] = -1.0;
                        }
                    }
                    for (int i = 0; i < dcomps.size(); ++i){
                        if (derprobs[dcomps[i].idx1] >= 0 &&
                            derprobs[dcomps[i].idx2] >= 0 &&
                            derprobs[dcomps[i].idxT] >= 0){

                            dcomp* dc = &dcomps[i];
                            // ABBA = probability derived in T & 2
                            double aprob = derprobs[dc->idxT];
                            double abba = derprobs[dc->idx2] * aprob;

                            // BABA = probability derived in T & 1
                            double baba = derprobs[dc->idx1] * aprob;
                            abcounts[i].first += abba;
                            abcounts[i].second += baba;
                            dcomp_inf_sites[i]++;
                                
                            if (winsize > 0 && !dc->for_f){
                                wbj_counts* wc = NULL;
                                if (wbj_dat.count(i) == 0){
                                    vector<wbj_counts> v;
                                    wbj_dat.insert(make_pair(i, v));
                                }
                                if (new_bin){
                                    // Need a new bin no matter what
                                    wbj_counts counts;
                                    wbj_dat[i].push_back(counts);
                                    wc = &wbj_dat[i][wbj_dat[i].size()-1];
                                    wc->chrom = chrom_idx;
                                    wc->bin = binstart;
                                }
                                else{
                                    // Might need a new bin
                                    bool make_new = false;
                                    if (wbj_dat[i].size() == 0){
                                        make_new = true;
                                    }
                                    else if (wbj_dat[i][wbj_dat[i].size()-1].chrom != chrom_idx ||
                                        wbj_dat[i][wbj_dat[i].size()-1].bin != binstart){
                                        make_new = true;
                                    }
                                    if (make_new){
                                        wbj_counts counts;
                                        wbj_dat[i].push_back(counts);
                                        wc = &wbj_dat[i][wbj_dat[i].size()-1];
                                        wc->chrom = chrom_idx;
                                        wc->bin = binstart;
                                    }
                                    else{
                                        wc = &wbj_dat[i][wbj_dat[i].size()-1];
                                    }
                                }
                                // Store relevant weighted block jackknife info
                                wc->abba += abba;
                                wc->baba += baba;
                                wc->informative_sites++;
                            }
                        }
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
}

double normcdf(double x, double mu, double sigma){
    return 0.5 * erfc(-((x-mu)/sigma) * M_SQRT1_2);
}

int main(int argc, char *argv[]) {    
    
    static struct option long_options[] = {
       {"vcf", required_argument, 0, 'v'},
       {"pops", required_argument, 0, 'p'},
       {"window", required_argument, 0, 'w'},
       {"f", no_argument, 0, 'f'},
       {0, 0, 0, 0} 
    };
    
    // Set default values
    string vcf_file = "";
    string popsfile = "";
    int window = 10000000;
    bool calc_f = false;

    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "v:p:w:fh", long_options, &option_index )) != -1){
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
            case 'w':
                window = (int)round(atof(optarg) * 1000000.0);
                break;
            case 'f':
                calc_f = true;
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
    if (popsfile == ""){
        fprintf(stderr, "ERROR: pops file required.\n");
        exit(1);
    }
    if (window <= 0){
        fprintf(stderr, "Disabled calculation of standard errors/p-values\n");
    }
    
    // Read in populations
    vector<string> pops1;
    vector<string> pops2;
    vector<string> popsT;
    vector<string> popsO;
    read_pops(popsfile, pops1, pops2, popsT, popsO);
    
    if (calc_f && popsT.size() < 2){
        fprintf(stderr, "ERROR: cannot estimate admixture proportions (--f) with \
fewer than 2 individuals present in the admixing population.\n");
        exit(1);
    }

    // Sort everything alphabetically so it looks good later
    sort(pops1.begin(), pops1.end());
    sort(pops2.begin(), pops2.end());
    sort(popsT.begin(), popsT.end());
    sort(popsO.begin(), popsO.end());

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
    for (int i = 0; i < num_samples; ++i){
        samp2idx.insert(make_pair(bcf_header->samples[i], i));
        idx2samp.push_back(bcf_header->samples[i]);
    }
    
    // Set up ABBA/BABA indices
    vector<pair<double, double> > abcounts;
    vector<dcomp> dcomps;
    vector<int> outgroup_idx;
    set<int> ingroup_all;

    for (int i = 0; i < popsO.size(); ++i){
        outgroup_idx.push_back(samp2idx[popsO[i]]);
    }

    int n_i = pops1.size();
    int n_j = pops2.size();
    if (n_i == 0){
        n_i = pops2.size();
    }
    else if (n_j == 0){
        n_j = pops1.size();
    }
    for (int i = 0; i < n_i; ++i){
        int i_x;
        if (pops1.size() == 0){
            if (samp2idx.count(pops2[i]) == 0){
                fprintf(stderr, "ERROR: sample %s not found in VCF.\n", pops2[i].c_str());
                exit(1);
            }
            i_x = samp2idx[pops2[i]];
        }
        else{
            if (samp2idx.count(pops1[i]) == 0){
                fprintf(stderr, "ERROR: sample %s not found in VCF.\n", pops1[i].c_str());
                exit(1);
            }
            i_x = samp2idx[pops1[i]];
        }
        for (int j = 0; j < n_j; ++j){
            if ((pops2.size() > 0 && pops1.size() > 0) ||
                j > i){
                int j_x;
                if (pops2.size() == 0){
                    if (samp2idx.count(pops1[j]) == 0){
                        fprintf(stderr, "ERROR: sample %s not found in VCF.\n", pops1[j].c_str());
                        exit(1);
                    }
                    j_x = samp2idx[pops1[j]];
                }
                else{
                    if (samp2idx.count(pops2[j]) == 0){
                        fprintf(stderr, "ERROR: sample %s not found in VCF.\n", pops2[j].c_str());
                    }
                    j_x = samp2idx[pops2[j]];
                }
                for (int k = 0; k < popsT.size(); ++k){
                    if (samp2idx.count(popsT[k]) == 0){
                        fprintf(stderr, "ERROR: sample %s not found in VCF.\n", popsT[k].c_str());
                    }
                    int k_x = samp2idx[popsT[k]];
                    ingroup_all.insert(i_x);
                    ingroup_all.insert(j_x);
                    ingroup_all.insert(k_x);
                    dcomp d(i_x, j_x, k_x);
                    dcomps.push_back(d);
                    abcounts.push_back(make_pair(0,0));
                }
            }
        }
    }
    
    // Map pop1 and pop2 individuals to all dcomps indices that
    // are to be used to calculate f-hat with those individuals
    map<int, vector<int> > pop1tof_d;
    map<int, vector<int> > pop2tof_d;

    if (calc_f){
        // Add in additional D stat calculations required to calculate f-hat
        for (int t1 = 0; t1 < popsT.size()-1; ++t1){
            int t1_x = samp2idx[popsT[t1]];
            for (int t2 = t1 + 1; t2 < popsT.size(); ++t2){
                int t2_x = samp2idx[popsT[t2]];

                if (pops1.size() > 0){
                    for (int i = 0; i < pops1.size(); ++i){
                        int i_x = samp2idx[pops1[i]];
                        if (pop1tof_d.count(i_x) == 0){
                            vector<int> v;
                            pop1tof_d.insert(make_pair(i_x, v));
                        }
                        dcomp d1(i_x, t1_x, t2_x);
                        d1.for_f = true;
                        dcomp d2(i_x, t2_x, t1_x);
                        d2.for_f = true;
                        pop1tof_d[i_x].push_back(dcomps.size());
                        dcomps.push_back(d1);
                        pop1tof_d[i_x].push_back(dcomps.size());
                        dcomps.push_back(d2);
                        abcounts.push_back(make_pair(0.0, 0.0));
                        abcounts.push_back(make_pair(0.0, 0.0));
                    }
                }
                if (pops2.size() > 0){
                    for (int j = 0; j < pops2.size(); ++j){
                        int j_x = samp2idx[pops2[j]];
                        if (pop2tof_d.count(j_x) == 0){
                            vector<int> v;
                            pop2tof_d.insert(make_pair(j_x, v));
                        }
                        dcomp d1(j_x, t1_x, t2_x);
                        d1.for_f = true;
                        dcomp d2(j_x, t2_x, t1_x);
                        d2.for_f = true;
                        pop2tof_d[j_x].push_back(dcomps.size());
                        dcomps.push_back(d1);
                        pop2tof_d[j_x].push_back(dcomps.size());
                        dcomps.push_back(d2);
                        abcounts.push_back(make_pair(0.0, 0.0));
                        abcounts.push_back(make_pair(0.0, 0.0));
                    }
                }
            }
        }
    }


    // For weighted block jackknife
    // BCF stores chromosomes as numeric indices called "rid"
    // chromosome -> window -> dcomp idx -> abba, baba
    map<int, vector<wbj_counts> > wbj_dat;
    vector<int> dcomp_inf_sites;
    for (int i = 0; i < dcomps.size(); ++i){
        dcomp_inf_sites.push_back(0);
    }

    // Ready to go
    read_vcf(bcf_reader, bcf_header, bcf_record, num_samples, dcomps, 
        dcomp_inf_sites, ingroup_all, outgroup_idx, abcounts, wbj_dat, window);
    
    string ogstr = "";
    for (int i = 0; i < popsO.size(); ++i){
        if (i > 0){
            ogstr += ",";
        }
        ogstr += popsO[i];
    }
    
    fprintf(stdout, "P1\tP2\ttest\toutgroup\tD");

    if (window > 0){
        fprintf(stdout, "\tSE\tZ\tp");
    }
    
    if (calc_f){
        fprintf(stdout, "\tf\tf_SE");
    }
    
    fprintf(stdout, "\n");

    for (int i = 0; i < dcomps.size(); ++i){
        
        // Visit all WBJ windows
        double abba = abcounts[i].first;
        double baba = abcounts[i].second;
        int inf_sites = dcomp_inf_sites[i];
        
        if (abba + baba == 0){
            continue;
        }
        
        string P1 = idx2samp[dcomps[i].idx1];
        string P2 = idx2samp[dcomps[i].idx2];
        string t = idx2samp[dcomps[i].idxT];
        
        double d = (abba - baba)/(abba + baba);
        
        double f;
        double f_se;
        double f_t;
        double f_t_se;
        if (calc_f){
            vector<double> fs;
            double fmean = 0.0;
            for (vector<int>::iterator idx = pop1tof_d[dcomps[i].idx1].begin(); 
                idx != pop1tof_d[dcomps[i].idx1].end(); ++idx){
                double abba_other = abcounts[*idx].first;
                double baba_other = abcounts[*idx].second;
                double d_other = (abba_other - baba_other)/(abba_other + baba_other);
                double f_this = d/d_other;
                fmean += f_this;
                fs.push_back(f_this);
            }
            fmean /= (double)fs.size();
            double fvar = 0.0;
            for (vector<double>::iterator f_i = fs.begin(); f_i != fs.end(); ++f_i){
                fvar += pow(*f_i - fmean, 2);
            }
            fvar /= (double)(fs.size()-1);
            f_se = sqrt(fvar);
            f = fmean;

            vector<double> f_ts;
            double f_tmean = 0.0;
            for (vector<int>::iterator idx = pop2tof_d[dcomps[i].idx2].begin();
                idx != pop2tof_d[dcomps[i].idx2].end(); ++idx){
                double abba_other = abcounts[*idx].first;
                double baba_other = abcounts[*idx].second;
                double d_other = (abba_other - baba_other)/(abba_other + baba_other);
                double f_this = -d/d_other;
                f_tmean += f_this;
                f_ts.push_back(f_this);
            }
            f_tmean /= (double)f_ts.size();
            double f_tvar = 0.0;
            for (vector<double>::iterator f_i = f_ts.begin(); f_i != f_ts.end(); ++f_i){
                f_tvar += pow(*f_i - fmean, 2);
            }
            f_tvar /= (double)(f_ts.size()-1);
            f_t_se = sqrt(f_tvar);
            f_t = f_tmean;
        }

        if (window <= 0){
            fprintf(stdout, "%s\t%s\t%s\t%s\t%e", P1.c_str(), P2.c_str(),
                t.c_str(), ogstr.c_str(), d);    
            if (calc_f){
                fprintf(stdout, "\t%e\t%e", f, f_se);
            }
            fprintf(stdout, "\n");
            if (pops1.size() == 0 || pops2.size() == 0){
                // Really just one pop. Print both directions.
                fprintf(stdout, "%s\t%s\t%s\t%s\t%e", P2.c_str(), P1.c_str(),
                    t.c_str(), ogstr.c_str(), -d);
                if (calc_f){
                    fprintf(stdout, "\t%e\t%e", f_t, f_t_se);
                }
                fprintf(stdout, "\n");
            }
        }
        else{
            double d_mu = 0.0;
            vector<double> wbj_d;
            vector<double> wbj_w;
            for (int j = 0; j < wbj_dat[i].size(); ++j){
                wbj_counts* wc = &wbj_dat[i][j];
                double weight = (double)wc->informative_sites / (double)inf_sites;
                double abba_adj = abba - wc->abba;
                double baba_adj = baba - wc->baba;
                double d = (abba_adj - baba_adj) / (abba_adj + baba_adj);
                d_mu += d;
                wbj_d.push_back(d);
                wbj_w.push_back(weight);
            }
            d_mu /= (double)wbj_d.size();
            // Now compute WBJ variance
            double d_var = 0.0;
            for (int j = 0; j < wbj_d.size(); ++j){
                d_var += wbj_w[j] * pow(wbj_d[j] - d_mu, 2);
            }
            double se = sqrt(d_var) / sqrt((double)wbj_d.size());
            double z = d/se;
            double zprime = z;
            if (z > 0){
                zprime = -z;
            }
            double p = 2.0 * normcdf(zprime, 0, 1); 
            fprintf(stdout, "%s\t%s\t%s\t%s\t%e\t%e\t%e\t%e", 
                P1.c_str(), P2.c_str(), t.c_str(), ogstr.c_str(),
                d, se, z, p);
            if (calc_f){
                fprintf(stdout, "\t%e\t%e", f, f_se);
            }
            fprintf(stdout, "\n");
            if (pops1.size() == 0 || pops2.size() == 0){
                // Really just one pop. Compute both directions.
                fprintf(stdout, "%s\t%s\t%s\t%s\t%e\t%e\t%e\t%e",
                    P2.c_str(), P1.c_str(), t.c_str(), ogstr.c_str(),
                    -d, se, -z, p);
                if (calc_f){
                    fprintf(stdout, "\t%e\t%e", f_t, f_t_se);
                }
                fprintf(stdout, "\n");
            }
        }
    }
}
