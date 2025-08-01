SHELL = bash
COMP = g++
CCOMP = gcc
PREFIX ?= /usr/local
CXXIFLAGS = -I$(PREFIX)/include -Iinclude
CIFLAGS = -I$(PREFIX)/include -Iinclude
LFLAGS = -L$(PREFIX)/lib -Llib
DEPS = -lz -lhts

all: ann_codon_pos data2ancestryhmm dstat fst mk sfs vcf2eigenstrat vcf2treemix vcf_filter

ann_codon_pos: src/ann_codon_pos.cpp
	$(COMP) src/ann_codon_pos.cpp -o ann_codon_pos $(DEPS)

data2ancestryhmm: src/data2ancestryhmm.cpp
	$(COMP) -O3 src/data2ancestryhmm.cpp -o data2ancestryhmm $(DEPS)

dstat: src/dstat.cpp
	$(COMP) -O3 src/dstat.cpp -o dstat $(DEPS)

fst: src/fst.cpp
	$(COMP) -O3 src/fst.cpp -o fst $(DEPS)

mk: src/mk.cpp
	$(COMP) src/mk.cpp -o mk $(DEPS)

sfs: src/sfs.cpp
	$(COMP) src/sfs.cpp -o sfs $(DEPS)

vcf2eigenstrat: src/vcf2eigenstrat.cpp
	$(COMP) -O3 src/vcf2eigenstrat.cpp -o vcf2eigenstrat $(DEPS)

vcf2treemix: src/vcf2treemix.cpp
	$(COMP) -O3 src/vcf2treemix.cpp -o vcf2treemix $(DEPS)

vcf_filter: src/vcf_filter.cpp
	$(COMP) -O3 src/vcf_filter.cpp -o vcf_filter $(DEPS)

clean:
	rm ann_codon_pos
	rm data2ancestryhmm
	rm dstat
	rm fst
	rm mk
	rm sfs
	rm vcf2eigenstrat
	rm vcf2treemix
	rm vcf_filter
