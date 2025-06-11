SHELL = bash
COMP = g++
CCOMP = gcc
PREFIX ?= /usr/local
CXXIFLAGS = -I$(PREFIX)/include -Iinclude
CIFLAGS = -I$(PREFIX)/include -Iinclude
LFLAGS = -L$(PREFIX)/lib -Llib
DEPS = -lz -lhts

all: ann_codon_pos dstat fst mk vcf_filter

ann_codon_pos: src/ann_codon_pos.cpp
	$(COMP) src/ann_codon_pos.cpp -o ann_codon_pos $(DEPS)

dstat: src/dstat.cpp
	$(COMP) -O3 src/dstat.cpp -o dstat $(DEPS)

fst: src/fst.cpp
	$(COMP) -O3 src/fst.cpp -o fst $(DEPS)

mk: src/mk.cpp
	$(COMP) src/mk.cpp -o mk $(DEPS)

vcf_filter: src/vcf_filter.cpp
	$(COMP) -O3 src/vcf_filter.cpp -o vcf_filter $(DEPS)

clean:
	rm ann_codon_pos
	rm dstat
	rm fst
	rm mk
	rm vcf_filter
