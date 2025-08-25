# hts_popgen
Filter VCFs &amp; calculate pop gen-relevant stuff from VCFs using HTSlib

## Requirements
* [htslib](https://github.com/samtools/htslib)
* gcc
* To plot VCF filtering results: [R](https://www.r-project.org/) and [ggplot2](https://ggplot2.tidyverse.org/)

## Build
* Clone this repository
* `cd hts_popgen && make`

## Programs

All programs will show options if run with `-h` or with no arguments.

### dstat
* Calculate [D-statistics](https://avianhybrids.wordpress.com/2019/11/09/d-statistics-for-dummies-a-simple-test-for-introgression/) for a set of individuals. Requires a VCF of variant data plus a tab-separated text file mapping individuals (column 1) to populations (column 2).
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
* Requires a VCF of genotype data and a tab-separated text file mapping individuals in the VCF to populations.
* If given a window size, computes significance using the weighted block jackknife, which is described in the D-statistic section in [this paper](https://pubmed.ncbi.nlm.nih.gov/26826668/).
* This program also has the option to compute $F_{ST}$ per site. Note that these per-SNP values can be very noisy and some can also be negative.

### vcf2eigenstrat
Converts a [VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf) file into [EIGENSTRAT](https://reich.hms.harvard.edu/software/InputFileFormats) format, used as input for multiple Reich lab programs, including [admixtools2](https://uqrmaie1.github.io/admixtools/), Eigensoft SmartPCA, [and others](https://reich.hms.harvard.edu/software)

## vcf2treemix
Converts a VCF file into format expected by [TreeMix](https://bitbucket.org/nygcresearch/treemix/wiki/Home)

### vcf_filter
* Filters VCF files using some common strategies (variant quality, genotype quality, probability of Hardy-Weinberg equilibrium given population identifiers), but also does per-sample depth filtering. Outputs a filtered VCF and a file of QC-relevant information.
* Per-sample sequencing depth histograms often have a very long right tail. Because of this, the program computes the histogram of depth per sample and finds the modal value for each sample. It then models depth as an exponential distribution with mean = $\frac{1}{modal value}$ and chooses the specified percentile of this exponential distribution as an upper cutoff for that sample (default = 0.995).
* For a lower-depth cutoff per sample, it uses a percentile of the genome-wide distribution (default = 0.005).
* Hard cutoffs for depth floor and ceiling are also used -- typically, the floor value will be used as the true cutoff (default = 5 reads) and the percentile-based cutoff described above will serve as the true upper cutoff (default upper ceiling = 1000 reads).
* Samples falling outside of their acceptable depth range as determined here will be set to missing.
* You can also provide a cutoff on the fraction of missing genotypes allowed per variant before the variant is discarded.
* A plotting script is included (`plot_vcf_filter.R`) that will plot information contained in the `.hist` file output from this program.
