Estimation of allele frequencies, SNP calling, and genotype calling from NGS data
=================================================================================

Make directories and set environmental variables for this session

```bash
mkdir ngs_analysis
mkdir ~/ngs_analysis/output
DIR=~/ngs_analysis
DATDIR=/ricco/data/tyler
BAMLIST=$DATDIR/cichlid_bams.list
CICHREF=$DATDIR/ref/GCA_900246225.3_fAstCal1.2_genomic_chromnames_mt.fa
BAMDIR=$DATDIR/bams
SCRIPTS=$DATDIR/scripts
ANGSD=/ricco/data/tyler/prog/bin/angsd
```

ANGSD overview

	$ANGSD --help


