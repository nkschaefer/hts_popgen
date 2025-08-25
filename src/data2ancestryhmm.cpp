#include <getopt.h>
#include <string>
#include <algorithm>
#include <vector>
#include <iterator>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <cstdio>
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
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/kstring.h>

using std::cout;
using std::endl;
using namespace std;

/**
 * Print a help message to the terminal and exit.
 */
void help(int code){
    fprintf(stderr, "data2ancestryhmm [OPTIONS]\n");
    fprintf(stderr, "Produces an input file for AncestryHMM (Corbett-Detig et al 2017)\n");
    fprintf(stderr, "given BAM files, or a VCF file plus and population identifiers for\n"); 
    fprintf(stderr, "reference populations and admixed individuals.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "===== OPTIONS =====\n");
    fprintf(stderr, "--out -o Output file prefix. Will create [out].ahmm and [out].sample.\n");
    fprintf(stderr, "----- Using BAM files as input -----\n");
    fprintf(stderr, "----- DO NOT USE THIS METHOD YET -----\n");
    fprintf(stderr, "This is liable to output way too many sites and is not suitable unless\n");
    fprintf(stderr, "you cannot call variants. This needs to be tested more before becoming\n");
    fprintf(stderr, "useful.\n");
    fprintf(stderr, "--refbams -b A 2 column tab separated text file, where column 1 is path to\n");
    fprintf(stderr, "    a BAM file and column 2 is the name of a reference population the\n");
    fprintf(stderr, "    individual belongs to. REQUIRED.\n");
    fprintf(stderr, "--admixbams -A A file listing text (admixed) bam files. Should be two columns,\n");
    fprintf(stderr, "    tab-separated. First column is file path; second column is individual name.\n");
    fprintf(stderr, "--mapqual -q Minimum map quality (default 30)\n");
    fprintf(stderr, "--basequal -Q Minimum base quality (default 30)\n");
    fprintf(stderr, "----- Using VCF file as input -----\n");
    fprintf(stderr, "--vcf -v VCF/BCF file. REQUIRED.\n");
    fprintf(stderr, "--refpops -p File mapping indv to pop. This should include all reference.\n");
    fprintf(stderr, "    admixer populations desired when running AncestryHMM. REQUIRED.\n");
    fprintf(stderr, "--admixed -a File listing test (admixed) individual names, one per line.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "With either input method, reference populations are listed in the output\n");
    fprintf(stderr, "    data file in alphabetical order (e.g. if providing CEU, JPT, and YRI, then\n");
    fprintf(stderr, "    CEU = 0, JPT = 1, and YRI = 2).\n");
    fprintf(stderr, "--help -h Display this message and exit.\n");
    exit(code);
}

void read_vcf(htsFile* bcf_reader, 
    bcf_hdr_t* bcf_header,
    bcf1_t* bcf_record,
    int num_samples,
    FILE* outf,
    map<int, int>& idx2pop,
    set<int>& idxadm){
    
    long int progress = 1000;
    long int nsnp = 0;
    
    int prevrid = -1;
    string curchrom = "";
    double prevr = 0;

    char genobuf[num_samples+1];
    genobuf[num_samples] = '\0';

    int snp_idx = 0;
    
    set<int> idxall;
    for (map<int, int>::iterator p = idx2pop.begin(); p != idx2pop.end(); ++p){
        idxall.insert(p->first);
    }
    for (set<int>::iterator a = idxadm.begin(); a != idxadm.end(); ++a){
        idxall.insert(*a);
    }

    while(bcf_read(bcf_reader, bcf_header, bcf_record) == 0){
        
        if (bcf_record->n_allele == 2){
            bcf_unpack(bcf_record, BCF_UN_STR);
            bool snp = true;
            for (int i = 0; i < bcf_record->n_allele; ++i){
                if (strcmp(bcf_record->d.allele[i], "A") != 0 &&
                    strcmp(bcf_record->d.allele[i], "C") != 0 &&
                    strcmp(bcf_record->d.allele[i], "G") != 0 && 
                    strcmp(bcf_record->d.allele[i], "T") != 0){
                    snp = false;
                    break;
                }
            }
            if (snp){
                // Get GTs            
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
                
                // Make sure it's biallelic in the sample.
                int nref = 0;
                int nalt = 0;
                for (set<int>::iterator idx = idxall.begin(); idx != idxall.end(); ++idx){
                    int32_t* ptr = gts + (*idx)*ploidy;
                    if (!bcf_gt_is_missing(ptr[0])){
                        for (int pi = 0; pi < ploidy; ++pi){
                            if (bcf_gt_allele(ptr[pi]) == 0){
                                nref++;
                            }
                            else{
                                nalt++;
                            }
                        }
                        if (nref > 0 && nalt > 0){
                            // It's biallelic
                            break;
                        }
                    }
                }
                if (nref > 0 && nalt > 0){
                    // Use this SNP.
                    string chrom = bcf_hdr_id2name(bcf_header, bcf_record->rid);
                    if (chrom != curchrom){
                        prevr = 0.0;
                        curchrom = chrom;
                    }               
                     
                    //int ploidy = 2;
                    int ploidy = n_gts / num_samples; 
                    // Assume 1cM/Mb recombination rate
                    double r = (double)(bcf_record->pos + 1) / (double)100000000;

                    fprintf(outf, "%s\t%ld", chrom.c_str(), bcf_record->pos + 1);
                    map<int, pair<int, int> > pcounts;

                    for (map<int, int>::iterator p = idx2pop.begin(); p != idx2pop.end(); ++p){
                        if (pcounts.count(p->second) == 0){
                            pcounts.insert(make_pair(p->second, make_pair(0,0)));
                        }
                        int32_t* ptr = gts + (p->first)*ploidy;
                        if (!bcf_gt_is_missing(ptr[0])){
                            for (int pi = 0; pi < ploidy; ++pi){
                                if (bcf_gt_allele(ptr[pi]) == 0){
                                    pcounts[p->second].first++;       
                                }
                                else{
                                    pcounts[p->second].second++;
                                }
                            }
                        }
                    }
                    
                    for (map<int, pair<int, int> >::iterator pc = pcounts.begin(); pc != pcounts.end();
                        ++pc){
                        fprintf(outf, "\t%d\t%d", pc->second.first, pc->second.second);    
                    }
                    
                    fprintf(outf, "\t%.8f", r - prevr);
                    prevr = r;

                    for (set<int>::iterator a = idxadm.begin(); a != idxadm.end(); ++a){
                        int refcount = 0;
                        int altcount = 0;
                        int32_t* ptr = gts + (*a)*ploidy;
                        if (!bcf_gt_is_missing(ptr[0])){
                            for (int pi = 0; pi < ploidy; ++pi){
                                if (bcf_gt_allele(ptr[pi]) == 0){
                                    refcount++;
                                }
                                else{
                                    altcount++;
                                }
                            }
                        }
                        fprintf(outf, "\t%d\t%d", refcount, altcount);
                    }
                    fprintf(outf, "\n");
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

void parse_pops(string& popfile, map<string, string>& indv2pop){
    ifstream inf(popfile);
    string indv;
    string pop;
    while (inf >> indv >> pop){
        indv2pop.insert(make_pair(indv, pop));
    }
}

void parse_adm(string& admfile, set<string>& adm_indv){
    ifstream inf(admfile);
    string indv;
    while (inf >> indv){
        adm_indv.insert(indv);
    }
}

// Info to pass to function (below) that will be called by pileup function
struct infile_t {

    const char* fname;
    samFile* fp;
    sam_hdr_t* fp_hdr;
    hts_itr_t* itr;
    hts_idx_t* idx;
    
    void init(const char* fname){
        this->fname = fname;
        
        this->fp = sam_open(fname, "r");
        this->fp_hdr = sam_hdr_read(this->fp);   
        
        // Load/create BAM index
        this->idx = hts_idx_load(fname, HTS_FMT_BAI);
        
        // Set iterator to use index 
        this->itr = sam_itr_queryi(this->idx, HTS_IDX_START, 0, INT32_MAX);
    };

    // Constructor
    infile_t(const char* fname){
        this->init(fname);  
    };

    infile_t(){
        fp = NULL;
        fp_hdr = NULL;
        itr = NULL;
        idx = NULL;
    }
    
    // Destructor
    ~infile_t(){
        sam_hdr_destroy(this->fp_hdr);
        sam_close(this->fp);
        hts_itr_destroy(this->itr);
        hts_idx_destroy(this->idx);
    }; 
};

// Iterator function for pileup
static int readaln(void *data, bam1_t *b){
    // Retrieve data in usable format
    infile_t *g = (infile_t*)data;
    int ret;
    while (1){
        ret = sam_itr_next(g->fp, g->itr, b); 
        if (ret < 0) break;
        if ( b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP) ) continue;
        break;
    }
    return ret;
}

/**
 * Log PDF of binomial distribution wrt n, k, p
 */
double logbinom(double n, double k, double p){
    double ll = k * log(p) + (n-k)*log(1.0-p);
    // Compute log binomial coefficient
    if (k < n && k != 0){
        // Use Stirling's approximation
        double logn = log(n);
        double logk = log(k);
        double logn_k = log(n-k);
        ll += n*logn - k*logk - (n-k)*logn_k + 0.5*(logn - logk - logn_k - log(2*M_PI));
    }
    return ll;
}

bool is_allele(double tot,
    double base_tot,
    vector<double>& base_files,
    vector<double>& tot_files,
    double p_error){
    
    // Require at least 10 reads for minor allele
    if (base_tot < 10){
        return false;
    }
    double ll_error = logbinom(tot, base_tot, p_error);
    double ll_var = logbinom(tot, base_tot, 0.5*(1.0/(double)tot_files.size()));
    if (ll_var > ll_error){
        for (int i = 0; i < base_files.size(); ++i){
            double p_private = 0.5;
            double ll_private = logbinom(tot_files[i], base_files[i], p_private);
            if (ll_private > logbinom(tot_files[i], base_files[i], p_error)){
                return true;
            }
        }    
    }
    return false;
    /*
    
    return false;
    */
}

void read_bams(vector<string>& bamfiles, 
    vector<int>& bamgroups, 
    FILE* outf,
    double p_error,
    int mapq,
    int baseq){
    
    int nfiles = bamfiles.size();
    
    bam_mplp_t mplp;
    
    //infile_t** data = (infile_t**)malloc(nfiles*sizeof(infile_t*));
    
    //vector<infile_t*> infiles;
    vector<void*> infiles;
    for (int i = 0; i < nfiles; ++i){
        infile_t* inf = new infile_t();
        inf->init(bamfiles[i].c_str());
        infiles.push_back((void*)inf);
        //data[i]->init(bamfiles[i].c_str());
    }
    
    //mplp = bam_mplp_init(nfiles, readaln, (void**)data);
    mplp = bam_mplp_init(nfiles, readaln, infiles.data());
    if (!mplp){
        fprintf(stderr, "ERROR initializing multi-pileup\n");
        exit(1);
    }
    int ovsuccess = bam_mplp_init_overlaps(mplp);
    if (ovsuccess < 0){
        fprintf(stderr, "ERROR ignoring mate pair overlaps\n");
        exit(1);
    }
    // Chromosome index - names in infile.fp_hdr->target_name array
    int tid = 0;
    // 0-based chromosome position
    int pos = 0;
    // Number of reads at position (per file)
    int n_plp[nfiles];
    // Pileup data (per file)
    const bam_pileup1_t** plp = (const bam_pileup1_t**)malloc(nfiles*sizeof(bam_pileup1_t*));

    int ret; 
    vector<double> acount_files(nfiles);
    vector<double> ccount_files(nfiles);
    vector<double> gcount_files(nfiles);
    vector<double> tcount_files(nfiles);
    
    vector<double> tot_files(nfiles);

    double prevr = 0.0;
    
    long int progress = 1000;
    long int nsnp = 0;
    
    int nrefpop = 0;
    vector<double> refpoptots;
    int bgprev = -1;
    for (int i = 0; i < bamgroups.size(); ++i){
        if (bamgroups[i] >= 0){
            if (bamgroups[i] != bgprev){
                ++nrefpop;
                refpoptots.push_back(0);
                bgprev = bamgroups[i];
            }
        }
        else{
            break;
        }
    } 
    int n = 0;
    while ((ret = bam_mplp_auto(mplp, &tid, &pos, &n_plp[0], plp)) > 0){
        if (tid < 0){
            // Unmapped reads
            break;
        }
        // Check base counts within files.
        double atot = 0.0;
        double ctot = 0.0;
        double gtot = 0.0;
        double ttot = 0.0;

        for (int i = 0; i < nrefpop; ++i){
            refpoptots[i] = 0.0;
        }
        
        for (int i = 0; i < nfiles; ++i) {
            
            acount_files[i] = 0.0;
            ccount_files[i] = 0.0;
            gcount_files[i] = 0.0;
            tcount_files[i] = 0.0;
            
            tot_files[i] = 0.0;

            for (int j = 0; j < n_plp[i]; ++j) {
                const bam_pileup1_t *p = &plp[i][j];

                uint8_t* seq = bam_get_seq(p->b);
                uint8_t* qual = bam_get_qual(p->b);
                uint16_t flag = p->b->core.flag;
                int del_len, is_rev = bam_is_rev(p->b);
                
                int mq = p->b->core.qual;
                int bq = qual[p->qpos] + 33;

                if (!p->is_del && p->indel == 0 && !p->is_refskip && mq >= mapq && bq >= baseq){
                    
                    // check map quality
                    double prob_map = 1.0 - pow(10, (double)p->b->core.qual / -10.0);
                    
                    // check base quality
                    unsigned char qualchr = qual[p->qpos] + 33;
                    double prob_base = 1.0 - pow(10, (double)qualchr / -10.0);
                    
                    //double count = prob_map*prob_base;
                    double count = 1.0;

                    // get base
                    unsigned char c = seq_nt16_str[bam_seqi(seq, p->qpos)];
                    switch(c){
                        case 'a':
                        case 'A':
                            acount_files[i] += count;
                            atot += count;
                            tot_files[i] += count;
                            if (bamgroups[i] >= 0){
                                refpoptots[bamgroups[i]] += count;
                            }
                            break;
                        case 'c':
                        case 'C':
                            ccount_files[i] += count;
                            ctot += count;
                            tot_files[i] += count;
                            if (bamgroups[i] >= 0){
                                refpoptots[bamgroups[i]] += count;
                            }
                            break;
                        case 'g':
                        case 'G':
                            gcount_files[i] += count;
                            gtot += count;
                            tot_files[i] += count;
                            if (bamgroups[i] >= 0){
                                refpoptots[bamgroups[i]] += count;
                            }
                            break;
                        case 't':
                        case 'T':
                            tcount_files[i] += count;
                            ttot += count;
                            tot_files[i] += count;
                            if (bamgroups[i] >= 0){
                                refpoptots[bamgroups[i]] += count;
                            }
                            break;
                    } 
                }        
            }
        }
        
        bool skipsite = false; 
        for (int i = 0; i < nrefpop; ++i){
            if (refpoptots[i] < 10.0){
                // Skip any variant that is low coverage in any reference population
                skipsite = true;
                break;
            }
        }
        if (skipsite){
            continue;
        }
        double tot = atot + ctot + gtot + ttot;
        
        int n_alleles = 0;
        bool a_is_allele = is_allele(tot, atot, acount_files, tot_files, p_error);
        if (a_is_allele){
            n_alleles++;
        }
        bool c_is_allele = is_allele(tot, ctot, ccount_files, tot_files, p_error);
        if (c_is_allele){
            n_alleles++;
        }
        bool g_is_allele = is_allele(tot, gtot, gcount_files, tot_files, p_error);
        if (g_is_allele){
            n_alleles++;
        }
        bool t_is_allele = is_allele(tot, ttot, tcount_files, tot_files, p_error);
        if (t_is_allele){
            n_alleles++;
        }
        
        // Only accept biallelic
        if (n_alleles == 2){
            
            // Go ahead and print data about the site.
            infile_t* ptr = (infile_t*)infiles[0];
            fprintf(outf, "%s\t%d", ptr->fp_hdr->target_name[tid], pos+1);

            // Find top two bases
            vector<pair<double, char> > sorts = {make_pair(-atot, 'A'), make_pair(-ctot, 'C'), 
                make_pair(-gtot, 'G'), make_pair(-ttot, 'T')};
            sort(sorts.begin(), sorts.end());
            
            // Process ref pops
            int idx = 0;
            int prevgrp = 9999;
            double count1 = 0.0;
            double count2 = 0.0;
            while (idx < nfiles){
                if (bamgroups[idx] < 0){
                    // Hitting the admixed individuals.
                    break;
                }
                else{
                    if (prevgrp != bamgroups[idx]){
                        if (prevgrp != 9999){
                            fprintf(outf, "\t%d\t%d", (int)round(count1), (int)round(count2));
                        }
                        count1 = 0.0;
                        count2 = 0.0;
                        prevgrp = bamgroups[idx];
                    }
                    if (sorts[0].second == 'A'){
                        count1 += acount_files[idx];
                    }
                    else if (sorts[0].second == 'C'){
                        count1 += ccount_files[idx];
                    }
                    else if (sorts[0].second == 'G'){
                        count1 += gcount_files[idx];
                    }
                    else if (sorts[0].second == 'T'){
                        count1 += tcount_files[idx];
                    }
                    if (sorts[1].second == 'A'){
                        count2 += acount_files[idx];
                    }
                    else if (sorts[1].second == 'C'){
                        count2 += ccount_files[idx];
                    }
                    else if (sorts[1].second == 'G'){
                        count2 += gcount_files[idx];
                    }
                    else if (sorts[1].second == 'T'){
                        count2 += tcount_files[idx];
                    }
                    ++idx;
                }
            }
            if (count1 + count2 > 0.0 && prevgrp != 9999){
                fprintf(outf, "\t%d\t%d", (int)round(count1), (int)round(count2));
            }

            // Print genetic distance
            // Assume 1 cM/Mb
            double r = (double)(pos + 1)/(double)(100000000);
            fprintf(outf, "\t%.8f", r - prevr);
            prevr = r;

            // Handle admixed individuals
            // For these, we don't need to sum reads across multiple individuals - each
            // one should be a single individual. 
            while (idx < nfiles){
                int count1;
                int count2;
                if (sorts[0].second == 'A'){
                    count1 = (int)round(acount_files[idx]);
                }
                else if (sorts[0].second == 'C'){
                    count1 = (int)round(ccount_files[idx]);
                }
                else if (sorts[0].second == 'G'){
                    count1 = (int)round(gcount_files[idx]);
                }
                else if (sorts[0].second == 'T'){
                    count1 = (int)round(tcount_files[idx]);
                }
                if (sorts[1].second == 'A'){
                    count2 = (int)round(acount_files[idx]);
                }
                else if (sorts[1].second == 'C'){
                    count2 = (int)round(ccount_files[idx]);
                }
                else if (sorts[1].second == 'G'){
                    count2 = (int)round(gcount_files[idx]);
                }
                else if (sorts[1].second == 'T'){
                    count2 = (int)round(tcount_files[idx]);
                }
                fprintf(outf, "\t%d\t%d", count1, count2);
                ++idx;
            }
            fprintf(outf, "\n");
        }

        nsnp++;
        if (nsnp % progress == 0){
            fprintf(stderr, "Read %ld sites\r", nsnp);
        }
    }
    
    fprintf(stderr, "Read %ld sites\n", nsnp);

    for (int i = 0; i < infiles.size(); ++i){
        infile_t* ptr = (infile_t*)infiles[i];
        delete ptr;

    }
    infiles.clear();

    free(plp);
    //free(data);

    /*
    vector<bam_hdr_t*> headers(nfiles);
    vector<hts_idx_t*> indices(nfiles);
    vector<void*> data(nfiles);
    vector<bam_plp_t> plps(nfiles);

    // Open BAMs and read headers
    for (int i = 0; i < nfiles; ++i) {

        puwrapper* puw = new puwrapper;
        data[i] = (void*)puw;

        // Set up single pileup for each BAM
        puw->fp = sam_open(bamfiles[i].c_str(), "r");
        headers[i] = sam_hdr_read(puw->fp);
        indices[i] = sam_index_load(puw->fp, bamfiles[i].c_str());
        puw->itr = sam_itr_queryi(indices[i], HTS_IDX_START, 0, INT32_MAX);
        
        plps[i] = bam_plp_init(fetch_func, data[i]);
    }
    
    // Multi-pileup
    bam_mplp_t mplp = bam_mplp_init(nfiles, [](void* d, bam1_t* b) {
        return fetch_func(d, b);
    }, data.data());

    int tid, pos, n_plp;
    const bam_pileup1_t** plp = (const bam_pileup1_t**) calloc(nfiles, sizeof(bam_pileup1_t*));
       while (bam_mplp_next(mplp, &tid, &pos, n_plp, plp) > 0) {
        
        
    }
    
    bam_mplp_destroy(mplp);
    for (int i = 0; i < nfiles; ++i) {
        bam_hdr_destroy(headers[i]);
        bam_plp_destroy(plps[i]);
        //hts_idx_destroy(indices[i]);
    }
    for (int i = 0; i < nfiles; ++i){
        delete (puwrapper*)data[i];
        hts_idx_destroy(plps[i]);
    }
    free(plp);
    //free(n_plp);
    */
}

int main(int argc, char *argv[]) {    
    
    static struct option long_options[] = {
       {"vcf", required_argument, 0, 'v'},
       {"refpop", required_argument, 0, 'p'},
       {"admixed", required_argument, 0, 'a'},
       {"refbams", required_argument, 0, 'b'},
       {"admixbams", required_argument, 0, 'A'},
       {"out", required_argument, 0, 'o'},
       {"mapqual", required_argument, 0, 'q'},
       {"basequal", required_argument, 0, 'Q'},
       {0, 0, 0, 0} 
    };
    
    // Set default values
    string vcf_file = "";
    string outprefix = "";
    string popfile = "";
    string admfile = "";
    string bampopfile = "";
    string bamadmfile = "";
    int mapqual = 30;
    int basequal = 30;

    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "v:o:p:a:b:A:q:Q:h", long_options, &option_index )) != -1){
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
                popfile = optarg;
                break;
            case 'o':
                outprefix = optarg;
                break;
            case 'a':
                admfile = optarg;
                break;
            case 'b':
                bampopfile = optarg;
                break;
            case 'A':
                bamadmfile = optarg;
                break;
            case 'q':
                mapqual = atoi(optarg);
                break;
            case 'Q':
                basequal = atoi(optarg);
                break;
            default:
                help(0);
                break;
        }    
    }
    
    // Error check arguments.
    bool bam_mode = false;
    bool vcf_mode = false;
    if (vcf_file != "" || popfile != "" || admfile != ""){
        vcf_mode = true;
    }
    if (bampopfile != "" || bamadmfile != ""){
        bam_mode = true;
    }
    if (bam_mode && vcf_mode){
        fprintf(stderr, "ERROR: input can be either BAM files or VCF, but not both.\n");
        exit(1);
    }
    if (!bam_mode && !vcf_mode){
        fprintf(stderr, "ERROR: you must either provide BAM or VCF input data.\n");
        exit(1);
    }
    if (vcf_mode){
        if (vcf_file == ""){
            fprintf(stderr, "ERROR: VCF/BCF file required\n");
            exit(1);
        }
        if (popfile == ""){
            fprintf(stderr, "ERROR: no reference populations given.\n");
            exit(1);
        }
        if (admfile == ""){
            fprintf(stderr, "ERROR: no admixed individuals given.\n");
            exit(1);
        }
    }
    if (bam_mode){
        if (bampopfile == ""){
            fprintf(stderr, "ERROR: --refbams/-b is required.\n");
            exit(1);
        }
        if (bamadmfile == ""){
            fprintf(stderr, "ERROR: --admixbams/-A is required.\n");
            exit(1);
        }
    }
    if (outprefix == ""){
        fprintf(stderr, "ERROR: output prefix is required.\n");
        exit(1);
    }
    
    string out_ahmm = outprefix + ".ahmm";
    FILE* outf = fopen(out_ahmm.c_str(), "w");
    if (outf == NULL){
        fprintf(stderr, "ERROR: cannot open file %s for writing.\n", out_ahmm.c_str());
        exit(1);
    }
    
    string out_sample = outprefix + ".sample";
    FILE* outs = fopen(out_sample.c_str(), "w");
    if (!outs){
        fprintf(stderr, "ERROR opening %s for writing.\n", out_sample.c_str());
        exit(1);
    }
    
    if (bam_mode){
        map<string, string> bam2pop;
        parse_pops(bampopfile, bam2pop);
        map<string, string> bam2adm;
        parse_pops(bamadmfile, bam2adm);
        
        // Write out sample file
        map<string, vector<string> > adm2bam;
        for (map<string, string>::iterator ba = bam2adm.begin(); ba != bam2adm.end(); ++ba){
            if (adm2bam.count(ba->second) == 0){
                vector<string> v;
                adm2bam.insert(make_pair(ba->second, v));
            }
            adm2bam[ba->second].push_back(ba->first);
        }
        for (map<string, vector<string> >::iterator ab = adm2bam.begin(); ab != adm2bam.end(); 
            ++ab){
            // Assume ploidy = 2
            fprintf(outs, "%s\t2\n", ab->first.c_str());
        }
        
        // Put files in order needed in ancestry_hmm output file
        vector<string> bamfiles;
        vector<int> bamgroups;

        map<string, vector<string> > pop2bam;
        for (map<string, string>::iterator p = bam2pop.begin(); p != bam2pop.end(); ++p){
            if (pop2bam.count(p->second) == 0){
                vector<string> v;
                pop2bam.insert(make_pair(p->second, v));
            }
            pop2bam[p->second].push_back(p->first);
        }
        int pi = 0;
        for (map<string, vector<string> >::iterator pb = pop2bam.begin(); pb != pop2bam.end(); ++pb){
            for (vector<string>::iterator b = pb->second.begin(); b != pb->second.end(); ++b){
                bamgroups.push_back(pi);
                bamfiles.push_back(*b);
            }
            ++pi;
        }
        int ai = -1;
        for (map<string, vector<string> >::iterator ab = adm2bam.begin(); ab != adm2bam.end(); ++ab){
            for (vector<string>::iterator b = ab->second.begin(); b != ab->second.end(); ++b){
                bamgroups.push_back(ai);
                bamfiles.push_back(*b);
            }
            --ai;
        }
        read_bams(bamfiles, bamgroups, outf, 0.01, mapqual, basequal);
    }
    else if (vcf_mode){
        map<string, string> indv2pop;
        parse_pops(popfile, indv2pop);
        set<string> adm_indv;
        parse_adm(admfile, adm_indv);

        // Translate to numeric.
        set<string> pops;
        for (map<string, string>::iterator p = indv2pop.begin(); p != indv2pop.end(); ++p){
            pops.insert(p->second);
        }
        int pi = 0;
        map<string, int> pop2idx;
        for (set<string>::iterator p = pops.begin(); p != pops.end(); ++p){
            pop2idx.insert(make_pair(*p, pi));
            pi++;
        }
        
        // Numeric indv idx mapped to pop idx    
        map<int, int> idx2pop_num;
        // Numeric indv idx of admixers
        set<int> adm_num;

        // Init BCF/VCF reader and take in header
        bcf_hdr_t* bcf_header;
        bcf1_t* bcf_record = bcf_init();
        htsFile* bcf_reader = bcf_open(vcf_file.c_str(), "r");
        if (bcf_reader == NULL){
            fprintf(stderr, "ERROR interpreting %s as BCF format.\n", vcf_file.c_str());
            exit(1);
        }
        bcf_header = bcf_hdr_read(bcf_reader);
        
        // Map samples to idx
        int nsamples = bcf_hdr_nsamples(bcf_header);
        for (int i = 0; i < nsamples; ++i) {
            string samp = bcf_header->samples[i];
            if (adm_indv.find(samp) != adm_indv.end()){
                if (indv2pop.count(samp) > 0){
                    fprintf(stderr, "ERROR: the following individual was listed both as candidate\n");
                    fprintf(stderr, "admixed and belonging to ref population %s:\n", indv2pop[samp].c_str());
                    fprintf(stderr, "%s\n", samp.c_str());
                    exit(1);
                }
                adm_num.insert(i);
                // Assume ploidy = 2
                fprintf(outs, "%s\t2\n", samp.c_str());
            }
            else if (indv2pop.count(samp) > 0){
                idx2pop_num.insert(make_pair(i, pop2idx[indv2pop[samp]]));
            }
        }
        
        if (adm_num.size() == 0){
            fprintf(stderr, "ERROR: no candidate admixers found in VCF.\n");
            exit(1);
        } 
        if (idx2pop_num.size() == 0){
            fprintf(stderr, "ERROR: no reference pop individuals found in VCF.\n");
            exit(1);
        }
        set<int> puniq;
        for (map<int, int>::iterator ix = idx2pop_num.begin(); ix != idx2pop_num.end(); ++ix){
            puniq.insert(ix->second);
        }
        if (puniq.size() < 2){
            fprintf(stderr, "ERROR: only %ld reference populations were specified that have\n", puniq.size());
            fprintf(stderr, "  individuals in the VCF.\n");
            exit(1);
        }
        
        // Ready to go
        read_vcf(bcf_reader, bcf_header, bcf_record, nsamples, outf, idx2pop_num, adm_num);
    }
    
    fclose(outs);
    fclose(outf);
}
