# hts_popgen
Filter VCFs &amp; calculate pop gen-relevant stuff from VCFs using HTSlib

These programs are meant to be fast and to require minimal effort from users to convert/create input files.

## Requirements
* [htslib](https://github.com/samtools/htslib)
* gcc
* To plot VCF filtering results: [R](https://www.r-project.org/) and [ggplot2](https://ggplot2.tidyverse.org/)

## Build
* Clone this repository
* `cd hts_popgen && make`

## Programs

All programs will show options if run with `-h` or with no arguments.

* [`ann_codon_pos`](#ann_codon_pos)
* [`data2ancestryhmm`](#data2ancestryhmm)
* [`dstat`](#dstat)
* [`fst`](#fst)
* [`mk`](#mk)
* [`sfs`](#sfs)
* [`vcf2eigenstrat`](#vcf2eigenstrat)
* [`vcf2treemix`](#vcf2treemix)
* [`vcf_filter`](#vcf_filter)
  
### ann_codon_pos
Adds codon position information to variants in a [VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf) file, where applicable. It requires a VCF file and an annotation (GFF3/GFF format). Wherever a variant affects the coding sequence of one or more transcripts, it adds a tag called "CP" (for codon position). The field is a comma-separated text string, where commas separate information for individual transcripts. For each transcript, information is pipe (|) separated, with the following fields:
1. Transcript ID
2. Gene ID (or blank if none)
3. Total number of codons in transcript CDS
4. 0-based index of codon in transcript CDS
5. 0-based index of SNP within codon (0, 1, or 2)

Output is a new VCF with CP tags added, which can then be summarized using the `mk` program, if desired.

### data2ancestryhmm
Prepares input data for [AncestryHMM](https://github.com/russcd/Ancestry_HMM), a tool that maps locus-specific ancestry across the genomes of hybrid individuals, also inferring the time of admixture. The tool is described [Corbett-Detig and Nielsen (2017)](https://pmc.ncbi.nlm.nih.gov/articles/PMC5242547/).

It can theoretically take BAM files as input, but this method is difficult, untested, and prone to outputting way too many variant sites. This is meant for cases where input data are scarce, but more testing is needed for this option to work.

With [VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf) file input, you can provide one file listing the names of all candidate admixed individuals (one per line), and another mapping individuals in the VCF to candidate "admixer" populations. This file should be tab-separated.

The output will include a file with extension `.ahmm`, which is the main input file to `Ancestry_HMM`. In this file, candidate admixer populations will be listed in alphabetical order. In other words, if you have provided an admixer population file listing the populations `CEU`, `JPT`, and `YRI`, then population 0 in the `.ahmm` file is `CEU`, population 1 is `JPT`, and population 2 is `YRI`, regardless of the order in which you listed these populations in the input file.

The output will also include a file with extension `.sample`, which just lists the candidate admixed populations and their ploidy (always assumed by this program to be 2. You can edit these files manually to alter ploidy.

To run `Ancestry_HMM`, you will need to specify parameters specific to your use case, but the `.ahmm` file will be the `-i` argument, the `.sample` file will be the `-s` argument, and I include the `--output_ancestry` flag. Results will be stored in `[name].posterior`, where `[name]` is an individual name from the `.sample` file. 

### dstat
* Calculate [D-statistics](https://avianhybrids.wordpress.com/2019/11/09/d-statistics-for-dummies-a-simple-test-for-introgression/) for a set of individuals. Requires a [VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf) of variant data plus a tab-separated text file mapping individuals (column 1) to populations (column 2).
* For more information on D-statistics, $\hat{f}$, and the weighted block jackknife, see the section in [this paper](https://pubmed.ncbi.nlm.nih.gov/26826668/).
* D measures whether individuals from a test population are closer to reference population 1 (negative) or reference population 2 (positive) by counting derived alleles shared with each of the two groups.
  * Populations are:
    * 1: Reference population 1, **A**BBA-**B**ABA (e.g. Yoruba in D(Yoruba, French, Neanderthal, Chimpanzee), which is positive and shows excess allele sharing between Neanderthal/French)
    * 2: Reference population 2, A**B**BA-B**A**BA (e.g. French in D(Yoruba, French, Neanderthal, Chimpanzee))
    * T: Test population, AB**B**A-BA**B**A (e.g. Neanderthal in D(Yoruba, French, Neanderthal, Chimpanzee))
    * O: Outgroup population, ABB**A**-BAB**A** (e.g. Chimpanzee in D(Yoruba, French, Neanderthal, Chimpanzee))
* All possible combinations of individuals will be tested, except for outgroup individuals -- D will only consider sites fixed in the outgroup for the sake of polarizing ancestral/derived alleles. In other words, if population 1 has 2 individuals a and b, population 2 has 2 individuals c and d, population 3 has 2 individuals e and f, and the outgroup has 2 individuals g and h, then it will calculate:
  * D(a,c,e,g+h)
  * D(b,c,e,g+h)
  * D(a,d,e,g+h)
  * D(b,d,e,g+h)
  * D(a,c,f,g+h)
  * D(b,c,f,g+h)
  * D(a,d,f,g+h)
  * D(b,d,f,g+h)
* If you omit either population **1** or **2** from the population file, then each individual from the included population will be used to stand in for the missing population. In other words, if you wanted to know if any humans appeared to match Neanderthal more closely than any other humans, you could provide a file mapping every human to population 2, and each configuration of D will be computed, with each individual used as population 1 in turn.
* The weighted block jackknife can be used to test significance - this is especially useful if you do not have many individuals to test. This should be set to a reasonably high number to overcome local LD. This allows for the variance and significance of estimates to be computed. If you have many individuals, this will slow things down and might be unnecessary (e.g. significance could be determined at the population level instead).
* You can also compute $\hat{f}$, an estimate of admixture proportion, if you have multiple **T** population individuals. This disables the weighted block jackknife computation.
  
### fst
* Computes [Wright's fixation index](https://en.wikipedia.org/wiki/Fixation_index) between populations using a VCF. Follows the [Hudson *et al* (1992)](https://pubmed.ncbi.nlm.nih.gov/1427045/) formulation after recommendations made by [Bhatia *et al* (2013)](https://doi.org/10.1101/gr.154831.113).
* Requires a [VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf) of genotype data and a tab-separated text file mapping individuals in the VCF to populations.
* If given a window size, computes significance using the weighted block jackknife, which is described in the D-statistic section in [this paper](https://pubmed.ncbi.nlm.nih.gov/26826668/).
* This program also has the option to compute $F_{ST}$ per site. Note that these per-SNP values can be very noisy and some can also be negative.

### mk
Given a [VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf) file annotated with codon position information (by `ann_codon_pos`) and a file mapping individuals to populations, tabulates statistics relevant to the [McDonald-Kreitman test](https://en.wikipedia.org/wiki/McDonald%25E2%2580%2593Kreitman_test) for each transcript. 
* The optional `--hybrid/-H` flag can be used if one or more populations is an F1 hybrid population. This will treat variants that are heterozygous in all individuals as fixed. This is pretty niche and irrelevant to most use cases.
* All statistics output approximate nonsynonymous mutations as those in codon positions 1 & 2 and synonymous as those in codon position 3. This may change in the future.
  
Output is a text-format table with the following fields:

* `transcript`: transcript ID
* `gene`: gene ID
* `pop1`: population 1 label
* `pop2`: population 2 label
* `Ntot`: Total sites in transcript CDS where nonsynonymous (codon position 1 & 2) mutations are possible
* `Stot`: Total sites in transcript CDS where synonymous (codon position 3) mutations are possible
* `dN`: Total nonsynonymous (codon position 1 & 2) mutations fixed between species (raw count, not normalized)
* `dS`: Total synonymous (codon position 3) mutations fixed between species (raw count, not normalized)
* `pN`: Total polymorphic mutations in either population (union) at codon positions 1 & 2 (raw count, not normalized)
* `pS`: Total polymorphic mutations in either population (union) at codon position 3 (raw count, not normalized)
* `NI`: McDonald-Kreitman test [neutrality index](https://academic.oup.com/mbe/article/28/1/63/986360)
* `logp`: Log-transformed Fisher's exact test p-value for comparing (dN, dS) to (pN, pS)
* `DoS`: [Direction of selection](https://academic.oup.com/mbe/article/28/1/63/986360) statistic for loci/populations with little data -- handles zeroes better than the neutrality index
* `alpha`: [Estimation of the proportion of mutations fixed by natural selection](https://www.nature.com/articles/4151022a)
* `pop1pN`: Polymorphisms at codon position 1 & 2 in population 1 (raw count, not normalized)
* `pop1pS`: Polymorphisms at codon position 3 in population 1 (raw count, not normalized)
* `pop2pN`: Polymorphisms at codon position 1 & 2 in population 2 (raw count, not normalized)
* `pop2pS`: Polymorphisms at codon position 3 in population 2 (raw count, not normalized)

### sfs
Computes site frequency spectrum (SFS)-based statistics that can be used to infer selection from sequence polymorphism without considering coding sequence. Takes a [VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf) file, [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) file, a file mapping individuals to populations (if desired), and, optionally, a way to infer ancestral versus derived alleles. Ancestral alleles can be inferred in the following ways:
* Specify an outgroup population. Sites fixed in this outgroup will be considered, and alleles that differ from the outgroup will be considered derived.
* `--ref_ancestral/-r` tells the program to assume that all reference alleles are ancestral (this could be appropriate, e.g. if you are using an imputed ancestral genome from [Cactus](https://github.com/ComparativeGenomicsToolkit/cactus) as the reference genome)
* Ancestral alleles are already encoded in the VCF via the `AA` tag, as in the [1000 Genomes Project data](https://www.internationalgenome.org/faq/how-do-you-calculate-ancestral-alleles/)
* If none of these options are available, only Tajima's D (which does not need to know about ancestral versus derived alleles) will be calculated.

Each statistic measures the population genetic parameter $\theta$ in two different ways and reports their difference. Calculation of $\theta$ can be distorted by demographic scenarios (like a population bottleneck followed by expansion) as well as natural selection. 

`sfs` outputs, per BED interval, a set of statistics:
* `chrom`: chromosome name
* `start`: 0-based start coordinate of feature interval
* `end`: 1-based end coordinate of feature interval
* `pop`: Population identifier (if given)
* `theta_pi`: Sequence polymorphism-based estimator of $\theta$, reflects mid-frequency polymorphisms
* `theta_w`: [Watterson estimator](https://en.wikipedia.org/wiki/Watterson_estimator) of $\theta$ based on number of segregating sites. Reflects low-frequency polymorphisms
* `D`: [Tajima's D](https://en.wikipedia.org/wiki/Tajima%27s_D), which compares `theta_pi` - `theta_w` (mid-to-low-frequency polymorphisms) and can be negative in cases of recent population expansion after a bottleneck, or selection.
* `D_dmin`: Tajima's D normalized by dividing by its minimum possible value, as proposed in [Schaeffer (2002)](https://www.cambridge.org/core/journals/genetics-research/article/molecular-population-genetics-of-sequence-length-diversity-in-the-adh-region-of-drosophila-pseudoobscura/C26FE8697447AC9ECBB2EC7FD08D683B)
* `Dprime`: Tajima's D normalized by dividing by its minimum possible value (equivalent to `D_dmin`) when negative, or its maximum possible value when positive
* `theta_h`: Estimator of $\theta$ that weights derived allele frequency; reflects high-frequency polymorphisms and proposed in [Fay and Wu (2000)](https://pmc.ncbi.nlm.nih.gov/articles/instance/1461156/pdf/10880498.pdf)
* `H`: [Fay and Wu's H](https://pmc.ncbi.nlm.nih.gov/articles/instance/1461156/pdf/10880498.pdf), which compares `theta_h` - `theta_pi` (high-to-mid-frequency polymorphisms) and can be positive in cases of positive selection
* `H_hmin`: Fay and Wu's H normalized by dividing by its minimum possible value, similar to `D_dmin`
* `Hprime`: Fay and Wu's H normalized by dividing by its minimum (when negative) or maximum (when positive) possible value
* `E` [Zeng _et al_'s E](https://pmc.ncbi.nlm.nih.gov/articles/PMC1667063/), which compares `theta_h` - `theta_w` (high-to-low-frequency polymorphisms) and can be positive in cases of positive selection
* `E_emin`: Zeng _et al_'s E normalized by its minimum possible value, like `D_dmin`
* `Eprime`: Zeng _et al_'s E normalized by its minimum (when negative) or maximum (when positive) possible value

**NOTE**: all statistics (`D`, `H`, and `E`) are normalized as suggested in [Zeng _et al_ (2006)](https://pmc.ncbi.nlm.nih.gov/articles/PMC1667063/)

### vcf2eigenstrat
Converts a [VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf) file into [EIGENSTRAT](https://reich.hms.harvard.edu/software/InputFileFormats) format, used as input for multiple Reich lab programs, including [admixtools2](https://uqrmaie1.github.io/admixtools/), Eigensoft SmartPCA, [and others](https://reich.hms.harvard.edu/software)

### vcf2treemix
Converts a [VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf) file into format expected by [TreeMix](https://bitbucket.org/nygcresearch/treemix/wiki/Home)

### vcf_filter
* Filters [VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf) files using some common strategies (variant quality, genotype quality, probability of Hardy-Weinberg equilibrium given population identifiers), but also does per-sample depth filtering. Outputs a filtered VCF and a file of QC-relevant information.
* Per-sample sequencing depth histograms often have a very long right tail. Because of this, the program computes the histogram of depth per sample and finds the modal value for each sample. It then models depth as an exponential distribution with mean = $\frac{1}{modal value}$ and chooses the specified percentile of this exponential distribution as an upper cutoff for that sample (default = 0.995).
* For a lower-depth cutoff per sample, it uses a percentile of the genome-wide distribution (default = 0.005).
* Hard cutoffs for depth floor and ceiling are also used -- typically, the floor value will be used as the true cutoff (default = 5 reads) and the percentile-based cutoff described above will serve as the true upper cutoff (default upper ceiling = 1000 reads).
* Samples falling outside of their acceptable depth range as determined here will be set to missing.
* You can also provide a cutoff on the fraction of missing genotypes allowed per variant before the variant is discarded.
* A plotting script is included (`plot_vcf_filter.R`) that will plot information contained in the `.hist` file output from this program.
