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
    fprintf(stderr, "mk [OPTIONS]\n");
    fprintf(stderr, "Given a VCF you have annotated to mark codon positions using ann_codon_pos,\n");
    fprintf(stderr, "plus a file mapping individuals in the VCF to populations, looks for fixed\n");
    fprintf(stderr, "differences between populations and polymorphisms within populations, separated\n");
    fprintf(stderr, "into 1st/2nd codon position mutations (stand-in for non-silent mutations) and\n");
    fprintf(stderr, "3rd codon position (stand-in for silent mutations).\n");
    fprintf(stderr, "Outputs McDonald-Kreitman test-relevant statistics per transcript.\n");
    fprintf(stderr, "Performs a McDonald-Kreitman test for each available transcript.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "The population file should be two-column, tab separated, where the first\n");
    fprintf(stderr, "column is individual names from the VCF, and the second is arbitrary text\n");
    fprintf(stderr, "population labels.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "===== OPTIONS =====\n");
    fprintf(stderr, "--vcf -v VCF/BCF file. REQUIRED.\n");
    fprintf(stderr, "--pops -p Population file (described above). REQUIRED.\n");
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

/**
 * Log binomial coefficient (n choose k)
 * Uses Stirling's approximation
 */
double binom_coef_log(double n, double k){
    if (k > n){
        return 0;
    }
    if (n == 0){
        return 0;
    }
    if (n == k){
        return 0;
    }
    if (k == 0){
        return 0; // 1 way to choose
    }
    double nf = (double)n;
    double kf = (double)k;

    // https://math.stackexchange.com/questions/64716/approximating-the-logarithm-of-the-binomial-coefficient
    // Use Stirling's approximation
    return nf*log(nf) - kf*log(kf) -
        (nf-kf)*log(nf-kf) + 0.5*(log(nf) - log(kf) -
        log(nf - kf) - log(2*M_PI));
}

double fisher_exact(double a1, double a2, double b1, double b2){
    // a = a1 c = a2 b = b1 d = b2
    double numer = binom_coef_log(a1 + a2, a1) + 
        binom_coef_log(b1 + b2, b1);
    double denom = binom_coef_log(a1 + a2 + b1 + b2, a1 + b1);
    return numer - denom;
}

void print_gene_aux(const string& txid,
    const string& gene,
    const string& pop1,
    const string& pop2,
    int dn,
    int ds,
    int pn,
    int ps,
    int pop1pn,
    int pop1ps,
    int pop2pn,
    int pop2ps,
    FILE* outf){
    
    // Eyre-Walker alpha
    double alpha = 1.0 - ((double)ds*(double)pn) / 
        ((double)dn*(double)ps);
    
    // Neutrality index
    double ni = ((double)pn / (double)ps) / 
        ((double)dn/(double)ds);
    
    double denom1 = (double)(dn + ds);
    if (denom1 == 0){
        denom1 = 1.0;
    }
    double denom2 = (double)(pn + ps);
    if (denom2 == 0){
        denom2 = 1.0;
    } 
    double dos = (double)dn/denom1 - 
       (double)pn/denom2;

    // Fisher's exact test p-value
    // order = dN, dS, pN, pS
    double logprob = fisher_exact(dn, ds, pn, ps);
    
    fprintf(outf, "%s\t%s\t%s\t%s", txid.c_str(),
        gene.c_str(), pop1.c_str(), 
        pop2.c_str());
    fprintf(outf, "\t%d\t%d\t%d\t%d",
        dn, ds, pn, ps);
    if (dn > 0 && ds > 0 && ps > 0){
        fprintf(outf, "\t%f", ni);
    }
    else{
        fprintf(outf, "\tNA");
    }
    fprintf(outf, "\t%f", logprob);
    if (true){
    //if (dn + ds > 0 && pn + ps > 0){
        fprintf(outf, "\t%f", dos);
    }
    else{
        fprintf(outf, "\tNA");
    }
    if (dn > 0 && ps > 0){
        fprintf(outf, "\t%f", alpha);
    }
    else{
        fprintf(outf, "\tNA");
    }
    fprintf(outf, "\t%d\t%d\t%d\t%d\n", pop1pn, pop1ps,
        pop2pn, pop2ps);
}

void print_genes(map<string, string>& tx2gene,
    vector<string>& popnames,
    map<string, map<int, pair<int, int> > >& tx_pop_poly,
    map<string, map<pair<int, int>, pair<int, int> > >& tx_2pop_poly,
    map<string, map<pair<int, int>, pair<int, int> > >& tx_pop_fixed,
    FILE* out_countsf){
    
    fprintf(out_countsf, "transcript\tgene\tpop1\tpop2\tdN\tdS\tpN\tpS\tNI\tlogp\tDoS\talpha\tpop1pN\tpop1pS\tpop2pN\tpop2pS\n");
    
    for (map<string, string>::iterator tg = tx2gene.begin(); tg != tx2gene.end(); ++tg){
        if (tx_pop_fixed.count(tg->first) > 0 || tx_pop_poly.count(tg->first) > 0){
            for (int i = 0; i < popnames.size()-1; ++i){
                int poly_i_12 = tx_pop_poly[tg->first][i].first;
                int poly_i_3 = tx_pop_poly[tg->first][i].second;

                for (int j = i + 1; j < popnames.size(); ++j){
                    int poly_j_12 = tx_pop_poly[tg->first][j].first;
                    int poly_j_3 = tx_pop_poly[tg->first][j].second;
                    
                    pair<int, int> key = make_pair(i,j);
                    
                    int poly_12 = tx_2pop_poly[tg->first][key].first;
                    int poly_3 = tx_2pop_poly[tg->first][key].second;

                    int fixed_12 = tx_pop_fixed[tg->first][key].first;
                    int fixed_3 = tx_pop_fixed[tg->first][key].second;
                    
                    print_gene_aux(tg->first,
                        tg->second,
                        popnames[i],
                        popnames[j],
                        fixed_12,
                        fixed_3,
                        poly_12,
                        poly_3,
                        poly_i_12,
                        poly_i_3,
                        poly_j_12,
                        poly_j_3,
                        out_countsf);  
                }
            }
            tx_pop_poly.erase(tg->first);
            tx_2pop_poly.erase(tg->first);
            tx_pop_fixed.erase(tg->first);
        }
    }
}


void read_vcf(htsFile* bcf_reader, 
    bcf_hdr_t* bcf_header,
    bcf1_t* bcf_record, 
    int num_samples,
    map<string, vector<int> >& pop_idx,
    map<string, string>& tx2gene,
    map<string, map<int, pair<int, int> > >& tx_pop_poly,
    map<string, map<pair<int, int>, pair<int, int> > >& tx_2pop_poly,
    map<string, map<pair<int, int>, pair<int, int> > >& tx_pop_fixed){

    long int progress = 1000;
    long int nsnp = 0;
    
    int prevrid = -1;
    string curchrom = "";
    
    string genomekey = "genome";
    {
        map<int, pair<int, int> > m;
        tx_pop_poly.insert(make_pair(genomekey, m));
        map<pair<int, int>, pair<int, int> > m2;
        tx_2pop_poly.insert(make_pair(genomekey, m2));
        tx_pop_fixed.insert(make_pair(genomekey, m2));
    }
    tx2gene.insert(make_pair(genomekey, genomekey));
    
    //char* tagbuf = (char*)malloc(1024*sizeof(char));
    char* tagbuf = NULL;
    char chunkbuf[1024];
    
    int tag_count = 0;

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
                
                int nchar = 0;

                // Try to get codon position info
                vector<string> txids;
                vector<int> cps;
                
                int success = bcf_get_info_string(bcf_header, bcf_record, "CP", &tagbuf, &nchar);
                if (success){
                    tag_count++;
                    int cbidx = 0;
                    int eltidx = 0;
                    for (int i = 0; i < nchar; ++i){
                        if (tagbuf[i] == ','){
                            chunkbuf[cbidx] = '\0';
                            // Last chunk should be codon pos
                            int cp = atoi(&chunkbuf[0]);
                            cps.push_back(cp);
                            cbidx = 0;    
                            eltidx = 0;
                        }
                        else if (tagbuf[i] == '|'){
                            chunkbuf[cbidx] = '\0';
                            if (eltidx == 0){
                                txids.emplace_back(&chunkbuf[0]);
                            }
                            else if (eltidx == 1){
                                if (cbidx > 0){
                                    string gid = &chunkbuf[0];
                                    tx2gene.insert(make_pair(txids[txids.size()-1], gid)); 
                                }
                                else{
                                    tx2gene.insert(make_pair(txids[txids.size()-1], ""));
                                }
                            }
                            cbidx = 0;
                            eltidx++;
                        }
                        else{
                            chunkbuf[cbidx] = tagbuf[i];
                            cbidx++;
                        }
                    }
                    // Last item should be codon pos
                    chunkbuf[cbidx] = '\0';
                    int cp = atoi(&chunkbuf[0]);
                    cps.push_back(cp);

                    for (int i = 0; i < txids.size(); ++i){
                    
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
                        
                        // Skip site if it's fixed in all individuals
                        int refcount = 0;
                        int altcount = 0;

                        vector<bool> pop_poly;
                        vector<short> pop_fixed_allele;

                        // First, check for whether the site is polymorphic within populations.
                        for (map<string, vector<int> >::iterator p = pop_idx.begin(); 
                            p != pop_idx.end(); ++p){
                            int pop_ref = 0;
                            int pop_alt = 0;
                            for (vector<int>::iterator i = p->second.begin(); i != p->second.end();
                                ++i){
                                int32_t* ptr = gts + (*i)*ploidy;
                                if (!bcf_gt_is_missing(ptr[0])){
                                    for (int pi = 0; pi < ploidy; ++pi){
                                        if (bcf_gt_allele(ptr[pi]) != 0){
                                            pop_alt++;
                                            altcount++;
                                        }
                                        else{
                                            pop_ref++;
                                            refcount++;
                                        }
                                    }
                                }
                            }
                            if (pop_ref + pop_alt == 0){
                                pop_poly.push_back(false);
                                pop_fixed_allele.push_back(-1);
                            }
                            else if (pop_ref > 0 && pop_alt > 0){
                                pop_poly.push_back(true);
                                pop_fixed_allele.push_back(-1);
                            }
                            else{
                               pop_poly.push_back(false);
                               if (pop_ref == 0){
                                   pop_fixed_allele.push_back(1);
                               }
                               else{
                                   pop_fixed_allele.push_back(0);
                               }
                            }
                        }
                        if (refcount > 0 && altcount > 0){
                            // The site is not fixed across all pops
                            // Make space in data structures
                            if (tx_pop_poly.count(txids[i]) == 0){
                                map<int, pair<int, int> > m1;
                                tx_pop_poly.insert(make_pair(txids[i], m1));
                                map<pair<int, int>, pair<int, int> > m2;
                                tx_pop_fixed.insert(make_pair(txids[i], m2));
                                tx_2pop_poly.insert(make_pair(txids[i], m2));
                            }
                            for (int j = 0; j < pop_fixed_allele.size(); ++j){
                                if (pop_poly[j]){
                                    if (tx_pop_poly[txids[i]].count(i) == 0){
                                        tx_pop_poly[txids[i]].insert(make_pair(j, make_pair(0,0)));
                                    }
                                    if (cps[i] == 2){
                                        tx_pop_poly[txids[i]][j].second++;
                                        tx_pop_poly[genomekey][j].second++;
                                    }
                                    else{
                                        tx_pop_poly[txids[i]][j].first++;
                                        tx_pop_poly[genomekey][j].first++;
                                    }
                                }

                                if (pop_fixed_allele[j] != -1){
                                    for (int k = j + 1; k < pop_fixed_allele.size(); ++k){
                                        if (pop_poly[j] || pop_poly[k]){
                                            pair<int, int> key = make_pair(j, k);
                                            if (tx_2pop_poly[txids[i]].count(key) == 0){
                                                tx_2pop_poly[txids[i]].insert(make_pair(
                                                    key, make_pair(0,0)));
                                            }
                                            if (cps[i] == 2){
                                                tx_2pop_poly[txids[i]][key].second++;
                                                tx_2pop_poly[genomekey][key].second++;
                                            }
                                            else{
                                                tx_2pop_poly[txids[i]][key].first++;
                                                tx_2pop_poly[genomekey][key].first++;
                                            }
                                        }
                                        else if (!pop_poly[k] && pop_fixed_allele[k] != -1 &&
                                            pop_fixed_allele[k] != pop_fixed_allele[j]){
                                            // Fixed difference.
                                            pair<int, int> key = make_pair(j, k);
                                            if (tx_pop_fixed[txids[i]].count(key) == 0){
                                                tx_pop_fixed[txids[i]].insert(make_pair(
                                                    key, make_pair(0,0)));
                                            }
                                            if (cps[i] == 2){
                                                tx_pop_fixed[txids[i]][key].second++;
                                                tx_pop_fixed[genomekey][key].second++;
                                            }
                                            else{
                                                tx_pop_fixed[txids[i]][key].first++;
                                                tx_pop_fixed[genomekey][key].first++;
                                            }
                                        }
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
    free(tagbuf);
    if (tag_count == 0){
        fprintf(stderr, "ERROR: no CP tags found. Did you run ann_codon_pos on this VCF?\n");
        exit(1);
    }
}

int main(int argc, char *argv[]) {    
    
    static struct option long_options[] = {
       {"vcf", required_argument, 0, 'v'},
       {"pops", required_argument, 0, 'p'},
       {0, 0, 0, 0} 
    };
    
    // Set default values
    string vcf_file = "";
    string popsfile = "";

    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "v:p:h", long_options, &option_index )) != -1){
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
    
    map<string, string> tx2gene;

    // Read in populations
    map<string, vector<string> > pops;
    if (popsfile != ""){
        read_pops(popsfile, pops);
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
    for (int i = 0; i < num_samples; ++i){
        samp2idx.insert(make_pair(bcf_header->samples[i], i));
        idx2samp.push_back(bcf_header->samples[i]);
    }
    
    // Set up population indices
    map<string, vector<int> > pop_idx;
    vector<string> idx2pop;
    for (map<string, vector<string> >::iterator p = pops.begin(); p != pops.end(); ++p){
        idx2pop.push_back(p->first);
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
    
    map<string, map<int, pair<int, int> > > tx_pop_poly;
    map<string, map<pair<int, int>, pair<int, int> > > tx_2pop_poly;
    map<string, map<pair<int, int>, pair<int, int> > > tx_pop_fixed;
    
    // Ready to go
    read_vcf(bcf_reader, bcf_header, bcf_record, num_samples, pop_idx,
        tx2gene, tx_pop_poly, tx_2pop_poly, tx_pop_fixed);
    
    vector<string> pop_names;
    map<string, int> pop_name2idx;

    for (map<string, vector<int> >::iterator pi = pop_idx.begin(); pi != 
        pop_idx.end(); ++pi){
        pop_name2idx.insert(make_pair(pi->first, pop_names.size()));
        pop_names.push_back(pi->first);
    }

    // Print data
    print_genes(tx2gene, pop_names, tx_pop_poly, tx_2pop_poly, tx_pop_fixed, stdout);   

}
