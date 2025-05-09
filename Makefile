SHELL = bash
COMP = g++
CCOMP = gcc
PREFIX ?= /usr/local
CXXIFLAGS = -I$(PREFIX)/include -Iinclude
CIFLAGS = -I$(PREFIX)/include -Iinclude
LFLAGS = -L$(PREFIX)/lib -Llib
DEPS = -lz -lhts

all: dstat fst vcf_filter

dstat: src/dstat.cpp
	$(COMP) -O3 src/dstat.cpp -o dstat $(DEPS)

fst: src/fst.cpp
	$(COMP) -O3 src/fst.cpp -o fst $(DEPS)

vcf_filter: src/vcf_filter.cpp
	$(COMP) -O3 src/vcf_filter.cpp -o vcf_filter $(DEPS)

clean:
	rm dstat
	rm fst
	rm vcf_filter
