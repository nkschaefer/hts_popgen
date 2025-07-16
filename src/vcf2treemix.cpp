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
    fprintf(stderr, "vcf2treemix [OPTIONS]\n");
    fprintf(stderr, "Produces a TreeMix input file given a VCF and individual-population mapping.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "===== OPTIONS =====\n");
    fprintf(stderr, "--vcf -v VCF/BCF file. REQUIRED.\n");
    fprintf(stderr, "--pop -p File mapping indv to pop. REQUIRED.\n");
    fprintf(stderr, "--out -o Output file (REQUIRED)\n");
    fprintf(stderr, "--allow_missing -M allow missing data (default = no missing data)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "--help -h Display this message and exit.\n");
    exit(code);
}

void read_vcf(htsFile* bcf_reader, 
    bcf_hdr_t* bcf_header,
    bcf1_t* bcf_record,
    map<int, int>& indv2pop_num,
    int num_samples,
    bool allow_missing,
    gzFile& outf){

    long int progress = 1000;
    long int nsnp = 0;
    
    int prevrid = -1;
    string curchrom = "";

    char linebuf[64];

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
                //int ploidy = 2;
                int ploidy = n_gts / num_samples; 
                
                // Store allele freqs per pop
                map<int, pair<int, int> > afs;
                bool any_miss = false;
                
                for (map<int, int>::iterator ip = indv2pop_num.begin(); ip != indv2pop_num.end(); ++ip){ 
                    int miss = 0;
                    
                    if (afs.count(ip->second) == 0){
                        afs.insert(make_pair(ip->second, make_pair(0,0)));
                    }

                    int32_t* ptr = gts + (ip->first)*ploidy;
                    if (!bcf_gt_is_missing(ptr[0])){
                        for (int pi = 0; pi < ploidy; ++pi){
                            if (bcf_gt_allele(ptr[pi]) == 0){
                                afs[ip->second].first++;
                            }
                            else{
                                afs[ip->second].second++;
                            }
                        }
                    }
                    else{
                        miss++;
                        any_miss = true;
                    }
                }
                bool first = true;
                if (allow_missing || !any_miss){
                    for (map<int, pair<int, int> >::iterator dat = afs.begin(); dat != afs.end(); ++dat){
                        if (!first){
                            gzwrite(outf, " ", 1);
                        }
                        sprintf(&linebuf[0], "%d,%d", dat->second.first, dat->second.second);
                        gzwrite(outf, &linebuf[0], strlen(linebuf));
                        first = false;
                    }
                    gzwrite(outf, "\n", 1);
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
       {"allow_missing", no_argument, 0, 'M'},
       {0, 0, 0, 0} 
    };
    
    // Set default values
    string vcf_file = "";
    string popfile = "";
    string outfile = "";
    bool allow_missing = false;

    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "v:p:o:Mh", long_options, &option_index )) != -1){
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
                outfile = optarg;
                break;
            case 'M':
                allow_missing = true;
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
    if (popfile == ""){
        fprintf(stderr, "ERROR: pop file required.\n");
        exit(1);
    }
    if (outfile == ""){
        fprintf(stderr, "ERROR: output file required\n");
        exit(1);
    }
    if (outfile.length() < 3 || outfile.substr(outfile.length()-3, 3) != ".gz"){
        outfile += ".gz";
    }

    map<string, string> indv2pop;
    parse_pops(popfile, indv2pop);
    map<string, int> pop2idx;
    set<string> pop;

    for (map<string, string>::iterator p = indv2pop.begin(); p != indv2pop.end(); ++p){
        pop.insert(p->second);
    }
    
    gzFile outf = gzopen(outfile.c_str(), "w"); 
    if (!outf){
        fprintf(stderr, "ERROR opening file %s for writing.\n", outfile.c_str());
        exit(1);
    }

    // Print header
    bool first = true;
    int pi = 0;
    for (set<string>::iterator p = pop.begin(); p != pop.end(); ++p){
        if (!first){
            gzwrite(outf, " ", 1);
        }
        gzwrite(outf, p->c_str(), p->length());
        first = false;
        pop2idx.insert(make_pair(*p, pi));
        pi++;
    }
    gzwrite(outf, "\n", 1);

    // Init BCF/VCF reader and take in header
    bcf_hdr_t* bcf_header;
    bcf1_t* bcf_record = bcf_init();
    htsFile* bcf_reader = bcf_open(vcf_file.c_str(), "r");
    if (bcf_reader == NULL){
        fprintf(stderr, "ERROR interpreting %s as BCF format.\n", vcf_file.c_str());
        exit(1);
    }
    bcf_header = bcf_hdr_read(bcf_reader);
    map<int, int> indv2pop_num;

    // Map samples to idx
    map<int, int> popcounts;
    int nsamples = bcf_hdr_nsamples(bcf_header);
    for (int i = 0; i < nsamples; ++i) {
        string samp = bcf_header->samples[i];
        if (indv2pop.count(samp) > 0){
            int popnum = pop2idx[indv2pop[samp]];
            indv2pop_num.insert(make_pair(i, popnum));
            if (popcounts.count(popnum) == 0){
                popcounts.insert(make_pair(popnum, 0));
            }
            popcounts[popnum]++;
        }
    }
    for (map<string, int>::iterator p = pop2idx.begin(); p != pop2idx.end(); ++p){
        if (popcounts.count(p->second) == 0 || popcounts[p->second] == 0){
            fprintf(stderr, "ERROR: no samples found for pop %s\n", p->first.c_str());
            exit(1);
        }
    } 

    // Ready to go
    read_vcf(bcf_reader, bcf_header, bcf_record, indv2pop_num, nsamples, allow_missing, outf);
    
    gzclose(outf);

}
