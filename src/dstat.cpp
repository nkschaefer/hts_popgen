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
#include <random>
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
    fprintf(stderr, "E (OPTIONAL): exclude sites that are derived in any individual from this\n");
    fprintf(stderr, "population\n");
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
    fprintf(stderr, "--site_test -s Instead of computing D, stores positions of ABBA sites in the\n");
    fprintf(stderr, "  genome for each comparison and then tests for clustering of sites (this can\n");
    fprintf(stderr, "  be used to further suggest admixture). Argument: a chromosome lengths file in\n");
    fprintf(stderr, "  the same format as required by BEDTools shuffle (chrom<tab>length)\n");
    fprintf(stderr, "--print_abba -P Print ABBA sites in BED format. This will be a lot of output\n");
    fprintf(stderr, "  especially if you are doing many comparisons.\n");
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

struct site_pos{
    int pos;
    double prob;
    short type;

    site_pos(int p, double P, short t){
        this->pos = p;
        this->prob = P;
        this->type = t;
    }
};

bool operator<(const site_pos& p1, const site_pos& p2){
    return p1.pos < p2.pos;
}

void read_pops(string& popsfile, 
    vector<string>& pop1,
    vector<string>& pop2,
    vector<string>& popT,
    vector<string>& popO,
    vector<string>& popE){

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
        else if (popname == "E" || popname == "e"){
            popE.push_back(indvname);
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
        fprintf(stderr, "WARNING: no outgroup individuals (pop = O) present.\n");
        //exit(1);
        fprintf(stderr, "This is required unless the ancestral allele (AA) field is populated\n");
        fprintf(stderr, "in the VCF. If this field is missing, you will get no results.\n");
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
    set<int>& excl_idx,
    vector<pair<double, double> >& abcounts,
    map<int, vector<wbj_counts> >& wbj_dat,
    int winsize,
    bool per_site,
    vector<map<int, vector<site_pos> > >& sitedat){
    
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
                if (outgroup_idx.size() == 0){
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
                }
                else{
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
                }
                if (ogfixed && anc_allele >= 0){
                    bool excl = false;
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
                            if (gt > 0 && excl_idx.find(n) != excl_idx.end()){
                                excl = true;
                                break;
                            }
                            //derprobs[*ig] = gt;
                            derprobs[n] = gt;
                        }
                        else{
                            //derprobs[*ig] = -1.0;
                            derprobs[n] = -1.0;
                        }
                    }
                    if (!excl){
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
                                
                                if (per_site && abba > 0 && baba == 0){
                                    // Store this site.
                                    if (sitedat[i].count(chrom_idx) == 0){
                                        vector<site_pos> v;
                                        sitedat[i].insert(make_pair(chrom_idx, v));
                                    }
                                    sitedat[i][chrom_idx].push_back(site_pos(bcf_record->pos, abba, 0));
                                }
                                else if (per_site && baba > 0 && abba == 0){
                                    // Store this site.
                                    if (sitedat[i].count(chrom_idx) == 0){
                                        vector<site_pos> v;
                                        sitedat[i].insert(make_pair(chrom_idx, v));
                                    }
                                    sitedat[i][chrom_idx].push_back(site_pos(bcf_record->pos, baba, 1));
                                }

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

void read_chromlens(string& genomefile, map<string, long int>& chromlens){
    ifstream inf(genomefile);
    string chrom;
    long int len;
    while (inf >> chrom >> len){
        chromlens.insert(make_pair(chrom, len));
    }
}

double score_dists(map<int, vector<site_pos> >& sitedat){
    
    // Sum of 1/(dist to closest variant)
    // Variants alone on chromosomes contribute 0
    
    vector<double> abbaruns;
    vector<double> babaruns;
    set<double> runs_all;
    short cur = 2;
    double prob = -1;
    
    double abbatot = 0;
    double babatot = 0;

    for (map<int, vector<site_pos> >::iterator x = sitedat.begin(); x != sitedat.end(); ++x){
        if (x->second.size() > 1){
            
            short cur = 2;
            double prob = -1;

            // Ensure it's sorted
            sort(x->second.begin(), x->second.end());
            for (int i = 0; i < x->second.size(); ++i){
                if (x->second[i].type != cur){
                    // New run.
                    if (prob > 0){
                        runs_all.insert(prob);
                        if (cur == 0){
                            abbatot++;
                            abbaruns.push_back(prob);
                        }
                        else if (cur == 1){
                            babatot++;
                            babaruns.push_back(prob);
                        }
                    }
                    cur = x->second[i].type;
                    prob = x->second[i].prob/0.5;
                }
                else{
                    prob += x->second[i].prob/0.5;
                }
            }
            
            if (prob != -1){
                runs_all.insert(prob);
                if (cur == 0){
                    abbatot++;
                    abbaruns.push_back(prob);
                }
                else if (cur == 1){
                    babatot++;
                    babaruns.push_back(prob);
                }
            }
        }
    }
    
    // Do a KS-test
    double abbamean = 0.0;
    double babamean = 0.0;

    map<double, double> abbacdf;
    map<double, double> babacdf;
    sort(abbaruns.begin(), abbaruns.end());
    double abba_cumulative = 0.0;
    set<double>::iterator cur_run = runs_all.begin();
    for (int i = 0; i < abbaruns.size(); ++i){
        abbamean += (1.0/abbatot)*abbaruns[i];
        abba_cumulative++;
        if (i == abbaruns.size()-1 || abbaruns[i] != abbaruns[i+1]){
            abbacdf.insert(make_pair(abbaruns[i], abba_cumulative/abbatot));
            while (cur_run != runs_all.end() && *cur_run < abbaruns[i]){
                abbacdf.insert(make_pair(*cur_run, abba_cumulative/abbatot));
                ++cur_run;
            }
            if (*cur_run == abbaruns[i]){
                ++cur_run;
            }
        }
    }
    while (cur_run != runs_all.end()){
        abbacdf.insert(make_pair(*cur_run, 1.0));
        ++cur_run;
    }
    sort(babaruns.begin(), babaruns.end());
    double baba_cumulative = 0.0;
    cur_run = runs_all.begin();
    for (int i = 0; i < babaruns.size(); ++i){
        babamean += (1.0/babatot)*babaruns[i];
        baba_cumulative++;
        if (i == babaruns.size()-1 || babaruns[i] != babaruns[i+1]){
            babacdf.insert(make_pair(babaruns[i], baba_cumulative/babatot));
            while (cur_run != runs_all.end() && *cur_run < babaruns[i]){
                babacdf.insert(make_pair(*cur_run, baba_cumulative/babatot));
                ++cur_run;
            }
            if (*cur_run == babaruns[i]){
                ++cur_run;
            }
        }
    }
    while (cur_run != runs_all.end()){
        babacdf.insert(make_pair(*cur_run, 1.0));
        ++cur_run;
    }
    double D = 0.0;
    for (set<double>::iterator r = runs_all.begin(); r != runs_all.end(); ++r){
        double diff = abs(abbacdf[*r] - babacdf[*r]);
        if (diff > D){
            D = diff;
        }
    }
    // Account for sample sizes
    double n = (double)abbaruns.size();
    double m = (double)babaruns.size();
    D = D / sqrt((n+m)/(n*m));
    double p = 2*exp(-2*D*D);
    double logp = log(2.0) - 2*D*D;
    fprintf(stderr, "MEAN %f %f\n", abbamean, babamean);
    return logp;
}

/*
double score_bootstrap(map<int, vector<int> >& abbasites,
    vector<string>& chromnames,
    map<string, long int>& chromlens){
    
    // First, randomly re-sample ABBA sites
    map<int, vector<int> > abbasites_sample;
    
    // Lay out genome largest -> smallest
    map<string, int> chrom2name;
    vector<pair<long int, string> > gsort;
    long int g = 0;

    for (int i = 0; i < chromnames.size(); ++i){
        chrom2name.insert(make_pair(chromnames[i], i));
        gsort.push_back(make_pair(-chromlens[chromnames[i]], chromnames[i]));
        g += chromlens[chromnames[i]];
    }
    sort(gsort.begin(), gsort.end());
    
    // Set up random sampler
    std::random_device rd;  // non-deterministic seed
    std::mt19937 gen(rd()); // Mersenne Twister RNG
    std::uniform_int_distribution<long int> dist(0, g - 1);
    
    for (map<int, vector<int> >::iterator x = abbasites.begin(); x != abbasites.end(); ++x){
        for (vector<int>::iterator y = x->second.begin(); y != x->second.end(); ++y){
            // Generate random base from the entire genome
            long int randbase = dist(gen);

            // Figure out where this base belongs in the genome.
            long int runtot = 0;
            bool found = false;
            for (int i = 0; i < gsort.size(); ++i){
                if (randbase >= runtot && randbase < runtot + -gsort[i].first){
                    // On this chrom.
                    int chrom_idx = chrom2name[gsort[i].second];
                    int chrom_pos = (int) (randbase - runtot);
                    if (abbasites_sample.count(chrom_idx) == 0){
                        vector<int> v;
                        abbasites_sample.insert(make_pair(chrom_idx, v));
                    }
                    abbasites_sample[chrom_idx].push_back(chrom_pos);
                    found = true;
                    break;
                }
                runtot += -gsort[i].first;
            }
            if (!found){
                // Shouldn't happen.
                fprintf(stderr, "ERROR: couldn't find random sample position %ld\n", randbase);
                exit(1);
            }
        }
    }
    
    return score_dists(abbasites_sample);
}
*/

void cluster_test(vector<map<int, vector<site_pos> > >& sitedat,
    vector<dcomp>& dcomps, 
    int reps,
    map<string, long int>& chromlens,
    vector<string>& chromnames,
    vector<string>& idx2samp,
    string& ogstr){
    
    // Perform independent test for each D comparison
    for (int i = 0; i < dcomps.size(); ++i){
        double s = score_dists(sitedat[i]);
        fprintf(stdout, "%s\t%s\t%s\t%s\t%f\n", idx2samp[dcomps[i].idx1].c_str(),
            idx2samp[dcomps[i].idx2].c_str(), idx2samp[dcomps[i].idxT].c_str(), ogstr.c_str(), s);
        continue;
        /*
        vector<double> score_boots;
        double scorebootmean = 0.0;
        for (int j = 0; j < reps; ++j){
            double score_samp = score_bootstrap(abbasites[i], chromnames, chromlens);
            score_boots.push_back(score_samp);
            scorebootmean += (1.0/(double)reps) * score_samp;
        }
        double var = 0.0;
        for (int j = 0; j < reps; ++j){
            var += (1.0/(double)(reps-1)) * pow(score_boots[j] - scorebootmean, 2);
        }
        double sd = sqrt(var);
        double p = 1.0 - normcdf(score_true, scorebootmean, sd);

        fprintf(stderr, "%s\t%s\t%s\t%s\t%f\t%f\t%f\t%f\n", idx2samp[dcomps[i].idx1].c_str(),
            idx2samp[dcomps[i].idx2].c_str(), idx2samp[dcomps[i].idxT].c_str(), ogstr.c_str(),
            score_true, scorebootmean, sd, p);
        */
    }
}

int main(int argc, char *argv[]) {    
    
    static struct option long_options[] = {
       {"vcf", required_argument, 0, 'v'},
       {"pops", required_argument, 0, 'p'},
       {"window", required_argument, 0, 'w'},
       {"f", no_argument, 0, 'f'},
       {"site_test", required_argument, 0, 's'},
       {"print_abba", no_argument, 0, 'P'},
       {0, 0, 0, 0} 
    };
    
    // Set default values
    string vcf_file = "";
    string popsfile = "";
    int window = 10000000;
    bool calc_f = false;
    bool per_site = false;
    bool print_abba = false;
    string genomefile = "";

    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "v:p:w:s:Pfh", long_options, &option_index )) != -1){
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
            case 's':
                per_site = true;
                genomefile = optarg;
                break;
            case 'P':
                print_abba = true;
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
    else if (calc_f){
        window = -1;
        fprintf(stderr, "Weighted block jackknife disabled for f calculations.\n");
    }
    if (per_site && calc_f){
        fprintf(stderr, "ERROR: only one of -s/-f is allowed.\n");
        exit(1);
    }
    if (per_site && print_abba){
        fprintf(stderr, "ERROR: only one of -s/-P is allowed.\n");
        exit(1);
    }
    if (calc_f && print_abba){
        fprintf(stderr, "ERROR: only one of -f/-P is allowed.\n");
        exit(1);
    }
    
    map<string, long int> chromlens;

    if (per_site){
        if (genomefile == ""){
            fprintf(stderr, "ERROR: genome (chromome lengths) file required for --site_test.\n");
            exit(1);
        }
        read_chromlens(genomefile, chromlens);
        window = -1;
        fprintf(stderr, "Weighted block jackknife disabled for site test.\n");
    }
    else if (print_abba){
        // Store individual ABBA site counts.
        per_site = true;
    }

    // Read in populations
    vector<string> pops1;
    vector<string> pops2;
    vector<string> popsT;
    vector<string> popsO;
    vector<string> popsE;
    read_pops(popsfile, pops1, pops2, popsT, popsO, popsE);
    
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
    set<int> excl_idx;

    for (int i = 0; i < popsO.size(); ++i){
        outgroup_idx.push_back(samp2idx[popsO[i]]);
    }
    for (int i = 0; i < popsE.size(); ++i){
        excl_idx.insert(samp2idx[popsE[i]]);
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
                    d.for_f = false;
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
                        
                        if (pops2.size() == 0){
                            if (pop2tof_d.count(i_x) == 0){
                                vector<int> v;
                                pop2tof_d.insert(make_pair(i_x, v));
                            }
                            pop2tof_d[i_x].push_back(dcomps.size()-2);
                            pop2tof_d[i_x].push_back(dcomps.size()-1);
                        }
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

                        if (pops1.size() == 0){
                            if (pop1tof_d.count(j_x) == 0){
                                vector<int> v;
                                pop1tof_d.insert(make_pair(j_x, v));
                            }
                            pop1tof_d[j_x].push_back(dcomps.size()-2);
                            pop1tof_d[j_x].push_back(dcomps.size()-1);
                        }
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
    
    vector<map<int, vector<site_pos> > > sitedat;
    if (per_site){
        for (int i = 0; i < dcomps.size(); ++i){
            map<int, vector<site_pos> > m;
            sitedat.push_back(m);
        }
    }

    // Ready to go
    read_vcf(bcf_reader, bcf_header, bcf_record, num_samples, dcomps, 
        dcomp_inf_sites, ingroup_all, outgroup_idx, excl_idx, abcounts, 
        wbj_dat, window, per_site, sitedat);
    
    string ogstr = "";
    for (int i = 0; i < popsO.size(); ++i){
        if (i > 0){
            ogstr += ",";
        }
        ogstr += popsO[i];
    }
    if (popsO.size() == 0){
        ogstr = "AA";
    }
    
    if (per_site){
        
        fprintf(stderr, "DCOMPS %ld\n", dcomps.size());

        vector<string> chromnames;
        int n_chroms = bcf_header->n[BCF_DT_CTG];
        for (int i = 0; i < n_chroms; ++i){
            const char *chrom_buf = bcf_hdr_id2name(bcf_header, i);
            string chrom_str = chrom_buf;
            chromnames.push_back(chrom_str);
        }
        if (print_abba){
            /*
            for (int i = 0; i < dcomps.size(); ++i){
                for (map<int, vector<pair<int, double> > >::iterator x = abbasites[i].begin(); x != abbasites[i].end(); ++x){
                    for (vector<pair<int, double> >::iterator y = x->second.begin(); y != x->second.end(); ++y){
                        fprintf(stdout, "%s\t%d\t%d\tABBA\t%f\tD(%s,%s,%s,(%s))\n", chromnames[x->first].c_str(),
                            y->first, y->first+1, y->second, idx2samp[dcomps[i].idx1].c_str(), 
                            idx2samp[dcomps[i].idx2].c_str(),
                            idx2samp[dcomps[i].idxT].c_str(), ogstr.c_str());
                    }
                }
                for (map<int, vector<pair<int, double> > >::iterator x = babasites[i].begin(); x != babasites[i].end(); ++x){
                    for (vector<pair<int, double> >::iterator y = x->second.begin(); y != x->second.end(); ++y){
                        fprintf(stdout, "%s\t%d\t%d\tBABA\t%f\tD(%s,%s,%s,(%s))\n", chromnames[x->first].c_str(),
                            y->first, y->first+1, y->second, idx2samp[dcomps[i].idx1].c_str(), 
                            idx2samp[dcomps[i].idx2].c_str(),
                            idx2samp[dcomps[i].idxT].c_str(), ogstr.c_str());
                    }
                }
                exit(0);
            }
            */
        }
        else{
            fprintf(stdout, "P1\tP2\ttest\toutgroup\tlogp\n");

            // Do a test per site.
            cluster_test(sitedat, dcomps, 100, chromlens, chromnames, idx2samp, ogstr);
        }
        // Done
        return 0;
    }

    fprintf(stdout, "P1\tP2\ttest\toutgroup\tD");

    if (window > 0 && !calc_f){
        fprintf(stdout, "\tSE\tZ\tp");
    }
    
    if (calc_f){
        fprintf(stdout, "\ttestA\ttestB\tf");
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
        if (calc_f && dcomps[i].for_f){
            continue;
        }

        string P1 = idx2samp[dcomps[i].idx1];
        string P2 = idx2samp[dcomps[i].idx2];
        string t = idx2samp[dcomps[i].idxT];
        
        double d = (abba - baba)/(abba + baba);
        
        if (calc_f){
            vector<double> fs;
            double fmean = 0.0;
            for (vector<int>::iterator idx = pop1tof_d[dcomps[i].idx1].begin(); 
                idx != pop1tof_d[dcomps[i].idx1].end(); ++idx){
                
                string ta = idx2samp[dcomps[*idx].idx2];
                string tb = idx2samp[dcomps[*idx].idxT];
                
                double abba_other = abcounts[*idx].first;
                double baba_other = abcounts[*idx].second;
                double d_other = (abba_other - baba_other)/(abba_other + baba_other);
                double f_this = d/d_other;
                
                fprintf(stdout, "%s\t%s\t%s\t%s\t%e\t%s\t%s\t%e\n", P1.c_str(), P2.c_str(),
                    t.c_str(), ogstr.c_str(), d, ta.c_str(), tb.c_str(), f_this);

            }
            for (vector<int>::iterator idx = pop2tof_d[dcomps[i].idx2].begin();
                idx != pop2tof_d[dcomps[i].idx2].end(); ++idx){
                
                string ta = idx2samp[dcomps[*idx].idx2];
                string tb = idx2samp[dcomps[*idx].idxT];

                double abba_other = abcounts[*idx].first;
                double baba_other = abcounts[*idx].second;
                double d_other = (abba_other - baba_other)/(abba_other + baba_other);
                //double f_this = -d/d_other;
                double f_this = d/d_other;
                
                fprintf(stdout, "%s\t%s\t%s\t%s\t%e\t%s\t%s\t%e\n", P1.c_str(), P2.c_str(),
                    t.c_str(), ogstr.c_str(), d, ta.c_str(), tb.c_str(), f_this);

                //fprintf(stdout, "%s\t%s\t%s\t%s\t%e\t%s\t%s\t%e\n", P2.c_str(), P1.c_str(),
                //    t.c_str(), ogstr.c_str(), -d, ta.c_str(), tb.c_str(), f_this);

            }
        }
        else{
            if (window <= 0){
                fprintf(stdout, "%s\t%s\t%s\t%s\t%e", P1.c_str(), P2.c_str(),
                    t.c_str(), ogstr.c_str(), d);    
                //if (calc_f){
                //    fprintf(stdout, "\t%e\t%e", f, f_se);
                //}
                fprintf(stdout, "\n");
                if (pops1.size() == 0 || pops2.size() == 0){
                    // Really just one pop. Print both directions.
                    fprintf(stdout, "%s\t%s\t%s\t%s\t%e", P2.c_str(), P1.c_str(),
                        t.c_str(), ogstr.c_str(), -d);
                    //if (calc_f){
                    //    fprintf(stdout, "\t%e\t%e", f_t, f_t_se);
                    //}
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
                //if (calc_f){
                //    fprintf(stdout, "\t%e\t%e", f, f_se);
                //}
                fprintf(stdout, "\n");
                if (pops1.size() == 0 || pops2.size() == 0){
                    // Really just one pop. Compute both directions.
                    fprintf(stdout, "%s\t%s\t%s\t%s\t%e\t%e\t%e\t%e",
                        P2.c_str(), P1.c_str(), t.c_str(), ogstr.c_str(),
                        -d, se, -z, p);
                    //if (calc_f){
                    //    fprintf(stdout, "\t%e\t%e", f_t, f_t_se);
                    //}
                    fprintf(stdout, "\n");
                }
            }
        }
    }
    return 0;
}
