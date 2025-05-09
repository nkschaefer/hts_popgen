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
    fprintf(stderr, "fst [OPTIONS]\n");
    fprintf(stderr, "Given a VCF and a file mapping individuals in the VCF to populations\n");
    fprintf(stderr, "computes Wright's fixation index (Fst) among pairs of populations, either\n");
    fprintf(stderr, "genome-wide or per SNP.\n");
    fprintf(stderr, "Uses the Hudson et al (1992) formulation after the recommendations made in\n");
    fprintf(stderr, "Bhatia et al (2013), https://doi.org/10.1101/gr.154831.113\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "The population file should be two-column, tab separated, where the first\n");
    fprintf(stderr, "column is individual names from the VCF, and the second is arbitrary text\n");
    fprintf(stderr, "population labels.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "===== OPTIONS =====\n");
    fprintf(stderr, "--vcf -v VCF/BCF file. REQUIRED.\n");
    fprintf(stderr, "--pops -p Population file (described above). REQUIRED.\n");
    fprintf(stderr, "--site -s Output statistic per site. Default = genome-wide.\n");
    fprintf(stderr, "--window -w Window size (in Mb) for weighted block jackknife. Default = 10\n");
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

double mean_fst(vector<pair<double, double> >& vec){
    double num_mean = 0.0;
    double denom_mean = 0.0;
    for (int i = 0; i < vec.size(); ++i){
        num_mean += vec[i].first;
        denom_mean += vec[i].second;
    }
    num_mean /= (double)vec.size();
    denom_mean /= (double)vec.size();
    return num_mean/denom_mean;
}

pair<double, double> mean_fst_separate(vector<pair<double, double> >& vec){
    double num_mean = 0.0;
    double denom_mean = 0.0;
    for (int i = 0; i < vec.size(); ++i){
        num_mean += vec[i].first;
        denom_mean += vec[i].second;
    }
    num_mean /= (double)vec.size();
    denom_mean /= (double)vec.size();
    return make_pair(num_mean, denom_mean);
}

void read_vcf(htsFile* bcf_reader, 
    bcf_hdr_t* bcf_header,
    bcf1_t* bcf_record, 
    int num_samples,
    map<string, vector<int> >& pop_idx,
    map<string, map<string, vector<pair<double, double> > > >& fsts,
    map<string, map<string, vector<double> > >& inf_sites,
    bool per_site,
    int window){
    
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
    
    map<int, map<int, vector<pair<double, double> > > > cur_block;
    map<int, map<int, pair<int, int> > > cur_loc;

    for (int i = 0; i < pop_names.size()-1; ++i){
        map<int, vector<pair<double, double> > > m1;
        cur_block.insert(make_pair(i, m1));
        map<int, pair<int, int> > m2;
        cur_loc.insert(make_pair(i, m2));
        for (int j = i + 1; j < pop_names.size(); ++j){
            vector<pair<double, double> > v;
            cur_block[i].insert(make_pair(j, v));
            cur_loc[i].insert(make_pair(j, make_pair(-1, -1))); 
        }
    }

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
                // Get current WBJ bin
                int binstart = (int)round(((double)bcf_record->pos / 
                    (double)window)) * (double)window;
                int chrom_idx = bcf_record->rid;
                bool new_bin = false;
                if (chrom_idx != prevrid || binstart != prevbin){
                    new_bin = true;
                }

                prevrid = chrom_idx;
                prevbin = binstart;

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
                
                // Store alt allele freq per population
                map<int, pair<double, double> > aaf;
                
                for (map<string, vector<int> >::iterator p = pop_idx.begin(); 
                    p != pop_idx.end(); ++p){
                    
                    double sum = 0.0;
                    double tot = 0.0;
                    for (vector<int>::iterator i = p->second.begin(); i != p->second.end();
                        ++i){
                        int32_t* ptr = gts + (*i)*ploidy;
                        if (!bcf_gt_is_missing(ptr[0])){
                            for (int pi = 0; pi < ploidy; ++pi){
                                if (bcf_gt_allele(ptr[pi]) != 0){
                                    sum++;
                                }
                                tot++;
                            }
                        }
                    }
                    if (tot > 0){
                        aaf.insert(make_pair(pop_name2idx[p->first], 
                            make_pair(sum, tot))); 
                    }
                }

                if (aaf.size() >= 2){

                    // Compute per-site F_ST for each population pair
                    for (int i = 0; i < pop_names.size()-1; ++i){
                        double f1 = aaf[i].first / aaf[i].second;
                        double n1 = aaf[i].second;
                        for (int j = i + 1; j < pop_names.size(); ++j){
                            if (aaf.count(i) > 0 && aaf.count(j) > 0){
                                double f2 = aaf[j].first / aaf[j].second;
                                if (f1 > 0 || f2 > 0){
                                    double n2 = aaf[j].second;
                                    double num = (f1-f2)*(f1-f2) - 
                                        f1*(1.0-f1)/(n1-1.0) - 
                                        f2*(1.0-f2)/(n2-1.0);
                                    double denom = f1 + f2 - 2.0*f1*f2;
                                    
                                    if (per_site & denom > 0){
                                        // Print data here and be done with it.
                                        string chrom = bcf_hdr_id2name(bcf_header, bcf_record->rid);
                                        // Convert to 1-based
                                        int pos = bcf_record->pos + 1;
                                        fprintf(stdout, "%s\t%d\t%s\t%s\t%f\n", chrom.c_str(),
                                            pos, pop_names[i].c_str(), pop_names[j].c_str(),
                                            num/denom);
                                    }
                                    else if (denom > 0){
                                        // Determine whether we're in a new bin.
                                        bool newbinthis = new_bin || cur_loc[i][j].first != chrom_idx ||
                                            cur_loc[i][j].second != binstart;
                                        if (newbinthis){
                                            // Get F_ST from current block
                                            if (cur_block[i][j].size() > 0){
                                                pair<double, double> mfst = mean_fst_separate(cur_block[i][j]);
                                                // Add to data structure
                                                if (fsts.count(pop_names[i]) == 0){
                                                    map<string, vector<pair<double, double> > > m;
                                                    fsts.insert(make_pair(pop_names[i], m));
                                                    map<string, vector<double> > m2;
                                                    inf_sites.insert(make_pair(pop_names[i], m2));
                                                }
                                                if (fsts[pop_names[i]].count(pop_names[j]) == 0){
                                                    vector<pair<double, double> > v;
                                                    fsts[pop_names[i]].insert(make_pair(pop_names[j],
                                                        v));
                                                    vector<double> v2;
                                                    inf_sites[pop_names[i]].insert(make_pair(pop_names[j],
                                                        v2));
                                                }
                                                fsts[pop_names[i]][pop_names[j]].push_back(mfst);
                                                inf_sites[pop_names[i]][pop_names[j]].push_back(
                                                    (double)cur_block[i][j].size());
                                            }
                                            // Make room for a new block
                                            cur_block[i][j].clear();
                                            cur_loc[i][j].first = chrom_idx;
                                            cur_loc[i][j].second = binstart;
                                        }
                                        // Add new data to current block
                                        cur_block[i][j].push_back(make_pair(num, denom));
                                    }
                                }
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

    // Add final blocks to data structure
    for (int i = 0; i < pop_names.size()-1; ++i){
        for (int j = i + 1; j < pop_names.size(); ++j){
            if (cur_block[i][j].size() > 0){
                pair<double, double> mfst = mean_fst_separate(cur_block[i][j]);
                // Add to data structure
                if (fsts.count(pop_names[i]) == 0){
                    map<string, vector<pair<double, double> > > m;
                    fsts.insert(make_pair(pop_names[i], m));
                    map<string, vector<double> > m2;
                    inf_sites.insert(make_pair(pop_names[i], m2));
                }
                if (fsts[pop_names[i]].count(pop_names[j]) == 0){
                    vector<pair<double, double> > v;
                    fsts[pop_names[i]].insert(make_pair(pop_names[j],
                        v));
                    vector<double> v2;
                    inf_sites[pop_names[i]].insert(make_pair(pop_names[j],
                        v2));
                }
                fsts[pop_names[i]][pop_names[j]].push_back(mfst);
                inf_sites[pop_names[i]][pop_names[j]].push_back(
                    (double)cur_block[i][j].size());
            }
        }
    }
}

double normcdf(double x, double mu, double sigma){
    return 0.5 * erfc(-((x-mu)/sigma) * M_SQRT1_2);
}



int main(int argc, char *argv[]) {    
    
    static struct option long_options[] = {
       {"vcf", required_argument, 0, 'v'},
       {"pops", required_argument, 0, 'p'},
       {"site", no_argument, 0, 'w'},
       {"window", required_argument, 0, 'w'},
       {0, 0, 0, 0} 
    };
    
    // Set default values
    string vcf_file = "";
    string popsfile = "";
    bool per_site = false;
    int window = 10000000;

    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "v:p:b:w:sh", long_options, &option_index )) != -1){
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
            case 's':
                per_site = true;
                break;
            case 'w':
                window = (int)round(atof(optarg)*1000000.0);
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
    
    srand(time(NULL));

    // Read in populations
    map<string, vector<string> > pops;
    read_pops(popsfile, pops);
    
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
    
    // Set up population indices
    map<string, vector<int> > pop_idx;
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
    } 
    
    // Create data structures for averaging values genome-wide
    map<string, map<string, vector<pair<double, double> > > > fsts;
    // Track the number of informative sites that went into each F_ST calculation
    // for each block of SNPs
    map<string, map<string, vector<double> > > inf_sites;
    
    if (per_site){
        fprintf(stdout, "chrom\tpos\tpop1\tpop2\tfst\n");
    }

    // Ready to go
    read_vcf(bcf_reader, bcf_header, bcf_record, num_samples, pop_idx,
        fsts, inf_sites, per_site, window);

    if (!per_site){
        bool hdr_printed = false;

        for (map<string, map<string, vector<pair<double, double> > > >::iterator x = 
            fsts.begin(); x != fsts.end(); ++x){
            for (map<string, vector<pair<double, double> > >::iterator y = x->second.begin();
                y != x->second.end(); ++y){
                
                if (y->second.size() == 1){
                    if (!hdr_printed){
                        fprintf(stdout, "pop1\tpop2\tfst\n");
                        hdr_printed = true;
                    }
                    fprintf(stdout, "%s\t%s\t%f\n", x->first.c_str(), y->first.c_str(),
                        y->second[0].first/y->second[0].second);
                }
                else{
                    if (!hdr_printed){
                        fprintf(stdout, "pop1\tpop2\tfst\tSE\tlower95CI\tupper95CI\n");
                        hdr_printed = true;
                    }
                    // Convert informative site counts to weights & compute overall mean
                    double fst_num = 0.0;
                    double fst_denom = 0.0;
                    double inf_site_tot = 0.0;
                    for (int i = 0; i < y->second.size(); ++i){
                        inf_site_tot += inf_sites[x->first][y->first][i];
                    }
                    vector<double> nums;
                    vector<double> denoms;
                    vector<double> ws;
                    for (int i = 0; i < y->second.size(); ++i){
                        double w = inf_sites[x->first][y->first][i]/inf_site_tot;
                        fst_num += y->second[i].first * w;
                        fst_denom += y->second[i].second * w;
                    }
                    double fst_mean = fst_num/fst_denom;
                    double fst_var = 0.0;
                    for (int i = 0; i < y->second.size(); ++i){
                        // Remove this window's contribution to the mean
                        // How much did it contribute?
                        double w = inf_sites[x->first][y->first][i]/inf_site_tot;
                        double fst_num_without = fst_num - w*y->second[i].first;
                        double fst_denom_without = fst_denom - w*y->second[i].second;
                        fst_num_without /= (1.0 - w);
                        fst_denom_without /= (1.0 - w);
                        double fst_jack = fst_num_without/fst_denom_without;
                        fst_var += w*pow(fst_jack - fst_mean, 2);
                    }
                    double fst_se = sqrt(fst_var);
                    fprintf(stdout, "%s\t%s\t%f\t%f\t%f\t%f\n", x->first.c_str(),
                        y->first.c_str(), fst_mean, fst_se, fst_mean-1.96*fst_se, 
                        fst_mean+1.96*fst_se);
                    /*
                    double fst_mean = mean_fst(y->second); 
                    
                    if (bootstraps > 0){
                        double bootmean = 0.0;
                        vector<double> fst_boots;
                        fprintf(stderr, "Bootstrap %s %s...\n",
                            x->first.c_str(), y->first.c_str());
                        for (int bs = 0; bs < bootstraps; ++bs){
                            fprintf(stderr, "Rep %d\r", bs);
                            vector<pair<double, double> > fst_boot;
                            for (int rep = 0; rep < y->second.size(); ++rep){
                                int ri = rand() % y->second.size();                    
                                fst_boot.push_back(y->second[ri]);
                            }
                            double fst_boot_mean = mean_fst(fst_boot);
                            fst_boots.push_back(fst_boot_mean);
                            bootmean += fst_boot_mean * (1.0/(double)bootstraps);
                        }
                        // Get boostrap variance
                        double bootvar = 0.0;
                        for (int bs = 0; bs < bootstraps; ++bs){
                            bootvar += pow(fst_boots[bs] - bootmean, 2);
                        }
                        bootvar /= (double)(bootstraps-1.0);
                        double se = sqrt(bootvar); 
                        fprintf(stderr, "\n");
                        fprintf(stdout, "%s\t%s\t%f\t%f\n", x->first.c_str(),
                            y->first.c_str(), fst_mean, se);
                    }
                    else{
                        fprintf(stdout, "%s\t%s\t%f\n", x->first.c_str(),
                            y->first.c_str(), fst_mean); 
                    }
                    */
                }
            }
        }
    }
}
