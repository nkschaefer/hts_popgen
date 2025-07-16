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

struct codon_pos{
    string txid;
    short frame;
    bool rev;
    int codon_idx;
    int txlen;
    codon_pos(const string& t, short f, bool r, int i, int l){
        txid = t;
        frame = f;
        rev = r;
        codon_idx = i;
        txlen = l;
    };
};


/**
 * Print a help message to the terminal and exit.
 */
void help(int code){
    fprintf(stderr, "ann_codon_pos [OPTIONS]\n");
    fprintf(stderr, "Given a VCF and a GFF/GFF3 annotation, adds codon position information to\n");
    fprintf(stderr, "each variant, where available, and outputs a new annotated VCF.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "===== OPTIONS =====\n");
    fprintf(stderr, "--vcf -v VCF/BCF file. REQUIRED.\n");
    fprintf(stderr, "--gff -g Annotation in GFF/GFF3 format. REQUIRED.\n");
    fprintf(stderr, "--out -o Output file (- for stdout). OPTIONAL (default = stdout)\n");
    fprintf(stderr, "--outfmt -O Output file format: v = vcf (default), z = gz compressed VCF, b = BCF\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "--help -h Display this message and exit.\n");
    exit(code);
}

void read_vcf(htsFile* bcf_reader, 
    bcf_hdr_t* bcf_header,
    bcf1_t* bcf_record,
    htsFile* outf,
    map<string, multimap<int, codon_pos> >& chrom_pos_tx_codon,
    map<string, string>& tx2gene){

    long int progress = 1000;
    long int nsnp = 0;
    
    int prevrid = -1;
    string curchrom = "";

    char framebuf[96];
    string framestr;

    while(bcf_read(bcf_reader, bcf_header, bcf_record) == 0){
        
        string chrom = bcf_hdr_id2name(bcf_header, bcf_record->rid);
        if (chrom_pos_tx_codon.count(chrom) > 0 && 
            chrom_pos_tx_codon[chrom].count(bcf_record->pos) > 0){
            
            // Add codon info to record.
            
            // codon info only really makes sense for SNPs
            if (bcf_record->n_allele > 1){
                
                bcf_unpack(bcf_record, BCF_UN_INFO);

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
                    // Good to add
                    string infostr = "";
                    pair<multimap<int, codon_pos>::iterator, multimap<int, codon_pos>::iterator> range = 
                        chrom_pos_tx_codon[chrom].equal_range(bcf_record->pos);
                    bool first = true;
                    for (multimap<int, codon_pos>::iterator it = range.first; it != range.second; ++it){
                        string infostr_elt = it->second.txid + "|";
                        if (tx2gene.count(it->second.txid) > 0){
                            infostr_elt += tx2gene[it->second.txid] + "|";
                        }
                        else{
                            infostr_elt += "|";
                        }
                        sprintf(&framebuf[0], "%d|%d|%d", it->second.txlen, it->second.codon_idx, it->second.frame);
                        framestr = framebuf;
                        infostr_elt += framestr;
                        if (!first){
                            infostr += ",";
                        }
                        infostr += infostr_elt;
                        first = false;
                    }
                    int success = bcf_update_info_string(bcf_header, bcf_record, "CP", infostr.c_str());
                    if (success < 0){
                        fprintf(stderr, "ERROR: unable to update codon pos tag\n");
                        exit(1);
                    }
                }
            }
        }
        // Write record.
        int success = bcf_write1(outf, bcf_header, bcf_record);
        
        nsnp++;
        if (nsnp % progress == 0){
            fprintf(stderr, "Read %ld snps\r", nsnp);
        }
    }
    fprintf(stderr, "Read %ld snps\n", nsnp);
}

/**
 * Get attributes from a feature in either GFF3 or GTF format
 */
void parse_attrs(string& attrs, map<string, string>& parsed){
    string field;
    istringstream splitter(attrs);
    while (getline(splitter, field, ';')){
        size_t eqpos = field.find("=");
        if (eqpos != string::npos){
            string key = field.substr(0, eqpos);
            string val = field.substr(eqpos + 1, field.length()-eqpos-1);
            if (val[0] == '"'){
                val = val.substr(1, val.length()-1);
            }
            if (val[val.length()-1] == '"'){
                val = val.substr(0, val.length()-1);
            }
            parsed.insert(make_pair(key, val));
        }
    }
}

void parse_gff(string& gff_file,
    map<string, multimap<int, codon_pos > >& dat,
    map<string, string>& tx2gene){
    
    string chrom;
    string source;
    string type;
    int start;
    int end;
    string score;
    string strand;
    string phase;
    string attrs;
    
    map<string, int> tx_cds_len; 
    map<string, string> tx2chrom;
    map<string, bool> txrev;
    map<string, vector<pair<int, int> > > tx2cds;
    map<string, string> attrs_parsed;

    ifstream inf(gff_file);
    string line;
    while (getline(inf, line)){
        string field;
        int idx = 0;
        istringstream splitter(line);
        while (getline(splitter, field, '\t')){
            switch(idx){
                case 0:
                    chrom = field;
                    break;
                case 1:
                    source = field;
                    break;
                case 2:
                    type = field;
                    break;
                case 3:
                    start = atoi(field.c_str());
                    break;
                case 4:
                    end = atoi(field.c_str());
                    break;
                case 5:
                    score = field;
                    break;
                case 6:
                    strand = field;
                    break;
                case 7:
                    phase = field;
                    break;
                case 8:
                    attrs = field;
                    break;
            }
            idx++;
        }
        attrs_parsed.clear();
        if (type == "CDS"){
            parse_attrs(attrs, attrs_parsed);
            string txid = "";
            if (attrs_parsed.count("Parent") > 0){
                txid = attrs_parsed["Parent"];
            }
            else if (attrs_parsed.count("parent") > 0){
                txid = attrs_parsed["parent"];
            }
            else if (attrs_parsed.count("transcript_id") > 0){
                txid = attrs_parsed["transcript_id"];
            }
            if (txid == ""){
                fprintf(stderr, "ERROR: unable to extract txid from %s\n", attrs.c_str());
                fprintf(stderr, "offending line:\n");
                fprintf(stderr, "%s\n", line.c_str());
                exit(1);
            }
            if (tx2cds.count(txid) == 0){
                vector<pair<int, int> > v;
                tx2cds.insert(make_pair(txid, v));
                tx_cds_len.insert(make_pair(txid, 0));
            }
            // Convert to BED-like coords (start is 0-based, end is 0-based+1)
            start -= 1;
            // Assume phase == 0 if missing (e.g. single-exon genes)
            int phase_int = 0;
            if (phase != "."){
                phase_int = atoi(phase.c_str());
            }
            //start += phase_int;
            tx2cds[txid].push_back(make_pair(start, end));
            tx_cds_len[txid] += (end-start);
        }
        else if (type == "transcript"){
            parse_attrs(attrs, attrs_parsed);
            string txid = "";
            string gid = "";
            if (attrs_parsed.count("ID") > 0){
                txid = attrs_parsed["ID"];
            }
            else if (attrs_parsed.count("transcript_id") > 0){
                txid = attrs_parsed["transcript_id"];
            }
            else if (attrs_parsed.count("id") > 0){
                txid = attrs_parsed["id"];
            }
            if (attrs_parsed.count("Parent") > 0){
                gid = attrs_parsed["Parent"];
            }
            else if (attrs_parsed.count("parent") > 0){
                gid = attrs_parsed["parent"];
            }
            else if (attrs_parsed.count("gene_id") > 0){
                gid = attrs_parsed["gene_id"];
            }
            if (txid == "" || gid == ""){
                fprintf(stderr, "ERROR: unable to parse transcript attrs: %s\n", attrs.c_str());
                fprintf(stderr, "offending line:\n");
                fprintf(stderr, "%s\n", line.c_str());
                exit(1);
            }
            tx2gene.insert(make_pair(txid, gid));
            tx2chrom.insert(make_pair(txid, chrom));
            if (strand == "+" || strand == "-"){
                bool rev = strand == "-";
                txrev.insert(make_pair(txid, rev));
            }
            else{
                fprintf(stderr, "ERROR: invalid strand info %s for transcript %s\n", strand.c_str(), txid.c_str());
                fprintf(stderr, "offending line:\n");
                fprintf(stderr, "%s\n", line.c_str());
                exit(1);
            }
        }
    }

    // Reformat
    for (map<string, vector<pair<int, int> > >::iterator x = tx2cds.begin(); x != 
        tx2cds.end(); ++x){
        string chrom = tx2chrom[x->first];
        if (dat.count(chrom) == 0){
            multimap<int, codon_pos> m;
            dat.insert(make_pair(chrom, m));
        }
        if (txrev[x->first]){
            int pos_cds = 0;
            for (int i = x->second.size()-1; i >= 0; i--){
                for (int pos_genome = x->second[i].second-1;
                    pos_genome >= x->second[i].first; --pos_genome){
                    short frame = pos_cds % 3;
                    dat[chrom].insert(make_pair(pos_genome, codon_pos(x->first, frame, true, pos_cds, 
                        tx_cds_len[x->first]/3)));
                    pos_cds++;            
                }
            }
        }
        else{
            int pos_cds = 0;
            for (vector<pair<int, int> >::iterator v = x->second.begin(); v != x->second.end();
                ++v){
                for (int pos_genome = v->first; pos_genome < v->second; ++pos_genome){
                    short frame = pos_cds % 3;
                    dat[chrom].insert(make_pair(pos_genome, codon_pos(x->first, frame, false, pos_cds,
                        tx_cds_len[x->first]/3)));
                    pos_cds++;
                }
            }
            
        }
    }
}

int main(int argc, char *argv[]) {    
    
    static struct option long_options[] = {
       {"vcf", required_argument, 0, 'v'},
       {"gff", required_argument, 0, 'g'},
       {"out", required_argument, 0, 'o'},
       {"outfmt", required_argument, 0, 'O'},
       {0, 0, 0, 0} 
    };
    
    // Set default values
    string vcf_file = "";
    string gff_file = "";
    string outfile = "-";
    string outfmt = "v";

    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "v:g:o:O:h", long_options, &option_index )) != -1){
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
            case 'g':
                gff_file = optarg;
                break;
            case 'o':
                outfile = optarg;
                break;
            case 'O':
                outfmt = optarg;
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
    if (gff_file == ""){
        fprintf(stderr, "ERROR: annotation file required.\n");
        exit(1);
    }
    if (outfmt != "v" && outfmt != "z" && outfmt != "b"){
        fprintf(stderr, "ERROR: unrecognized output format %s\n", outfmt.c_str());
        exit(1);
    }

    map<string, multimap<int, codon_pos > > chrom_pos_tx_codon;
    map<string, string> tx2gene;
    parse_gff(gff_file, chrom_pos_tx_codon, tx2gene);
    
    // Init BCF/VCF reader and take in header
    bcf_hdr_t* bcf_header;
    bcf1_t* bcf_record = bcf_init();
    htsFile* bcf_reader = bcf_open(vcf_file.c_str(), "r");
    if (bcf_reader == NULL){
        fprintf(stderr, "ERROR interpreting %s as BCF format.\n", vcf_file.c_str());
        exit(1);
    }
    bcf_header = bcf_hdr_read(bcf_reader);
    
    // Add description of new field to header.
    string cp_line = "##INFO=<ID=CP,Number=.,Type=String,Description=\"Codon position: 'Transcript ID | Gene ID | Codon position (1,2,3)'\">";
    int success = bcf_hdr_append(bcf_header, cp_line.c_str());
    if (success != 0){
        fprintf(stderr, "ERROR appending INFO field to header\n");
        exit(1);
    }
    success = bcf_hdr_sync(bcf_header);
    
    // Open output VCF (stdout)
    // Write gz-compressed by default
    htsFile* outf;
    if (outfmt == "v"){
        outf = hts_open(outfile.c_str(), "w");
    }
    else if (outfmt == "z"){
        outf = hts_open(outfile.c_str(), "wz");
    }
    else if (outfmt == "b"){
        outf = hts_open(outfile.c_str(), "wb");
    }
    else{
        fprintf(stderr, "ERROR: unrecognized output format %s\n", outfmt.c_str());
        exit(1);
    } 
    int write_success = bcf_hdr_write(outf, bcf_header);

    // Ready to go
    read_vcf(bcf_reader, bcf_header, bcf_record, outf,
        chrom_pos_tx_codon, tx2gene);
    
    hts_close(outf);

}
