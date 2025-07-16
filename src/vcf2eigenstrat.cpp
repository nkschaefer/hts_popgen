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
    fprintf(stderr, "vcf2eigenstrat [OPTIONS]\n");
    fprintf(stderr, "Produces EIGENSTRAT format data (e.g. for Admixtools or SMARTPCA)\n");
    fprintf(stderr, "given a VCF and individual-population mapping.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "===== OPTIONS =====\n");
    fprintf(stderr, "--vcf -v VCF/BCF file. REQUIRED.\n");
    fprintf(stderr, "--pop -p File mapping indv to pop. STRONGLY RECOMMENDED.\n");
    fprintf(stderr, "    If you omit --pop, a unique \"population\" label will be added to each individual.\n");
    fprintf(stderr, "--out -o Output file REQUIRED.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "--help -h Display this message and exit.\n");
    exit(code);
}

void read_vcf(htsFile* bcf_reader, 
    bcf_hdr_t* bcf_header,
    bcf1_t* bcf_record,
    int num_samples,
    FILE* outf_geno,
    FILE* outf_snp,
    vector<int>& indv_idx){

    long int progress = 1000;
    long int nsnp = 0;
    
    int prevrid = -1;
    string curchrom = "";

    char genobuf[num_samples+1];
    genobuf[num_samples] = '\0';

    int snp_idx = 0;

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
                for (vector<int>::iterator idx = indv_idx.begin(); idx != indv_idx.end(); ++idx){
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
                    //int ploidy = 2;
                    int ploidy = n_gts / num_samples; 
                    
                    // Assume 1cM/Mb recombination rate
                    double r = (double)(bcf_record->pos + 1) / (double)100000000;

                    // Use 1-based chromosome index as standin for chrom name

                    fprintf(outf_snp, "snp%d\t%d\t%.8f\t%ld\t%c\t%c\n", snp_idx, bcf_record->rid + 1,
                        r, bcf_record->pos + 1, bcf_record->d.allele[0][0], bcf_record->d.allele[1][0]);
                    
                    ++snp_idx;
                    
                    int i = 0;
                    for (vector<int>::iterator idx = indv_idx.begin(); idx != indv_idx.end(); ++idx){
                        int32_t* ptr = gts + (*idx)*ploidy;
                        if (bcf_gt_is_missing(ptr[0])){
                            genobuf[i] = '9';
                        }
                        else{
                            int ac = 0;
                            for (int pi = 0; pi < ploidy; ++pi){
                                if (bcf_gt_allele(ptr[pi]) != 0){
                                    ac++;
                                }
                            }
                            sprintf(&genobuf[i], "%d", ac);
                        }
                        ++i;
                    }
                    
                    fprintf(outf_geno, "%s\n", genobuf);
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

int main(int argc, char *argv[]) {    
    
    static struct option long_options[] = {
       {"vcf", required_argument, 0, 'v'},
       {"pop", required_argument, 0, 'p'},
       {"out", required_argument, 0, 'o'},
       {0, 0, 0, 0} 
    };
    
    // Set default values
    string vcf_file = "";
    string outprefix = "";
    string popfile = "";

    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "v:o:p:h", long_options, &option_index )) != -1){
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
    if (outprefix == ""){
        fprintf(stderr, "ERROR: output file required\n");
        exit(1);
    }
    map<string, string> indv2pop;
    if (popfile == ""){
        fprintf(stderr, "WARNING: no pops given. Adding a unique population ID for each individual.\n");
    }
    else{
        parse_pops(popfile, indv2pop);
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
    
    vector<int> keep_indvs;

    string ind_name = outprefix + ".ind";
    FILE* outf = fopen(ind_name.c_str(), "w");
    // Map samples to idx
    int nsamples = bcf_hdr_nsamples(bcf_header);
    for (int i = 0; i < nsamples; ++i) {
        if (indv2pop.size() > 0){
            string samp = bcf_header->samples[i];
            if (indv2pop.count(samp) > 0){
                fprintf(outf, "%s\tU\t%s\n", bcf_header->samples[i], indv2pop[samp].c_str());
                keep_indvs.push_back(i);
            }
        }
        else{
            fprintf(outf, "%s\tU\t%s\n", bcf_header->samples[i], bcf_header->samples[i]);
            keep_indvs.push_back(i);
        }
    }
    fclose(outf);
    
    string geno_name = outprefix + ".geno";
    string snp_name = outprefix + ".snp";
    FILE* outf_geno = fopen(geno_name.c_str(), "w");
    FILE* outf_snp = fopen(snp_name.c_str(), "w");
    
    // Ready to go
    read_vcf(bcf_reader, bcf_header, bcf_record, nsamples, outf_geno, outf_snp, keep_indvs);
    
    fclose(outf_geno);
    fclose(outf_snp);

}
