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

Today will be performing the first steps of most any NGS data analysis with [ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD)

ANGSD overview

	$ANGSD --help

Many of the calculation you do with ANGSD with be based on genotype likelihoods (GL)s,
which will improve inference accuracy over hard calline genotypes (avoids snowballing errors).
<br>
Let's calculate genotype likelihoods for the cichlids across the same 50 kb region that you filtered yesterday.
We'll also apply some quality control filters when calculating GLs from the bam files. You can see a full list 
options by running

	$ANGSD -bam

<details>

<summary> click here for bam parshing options </summary>

```bash
$ $ANGSD -bam
	-> angsd version: 0.935-48-gff0c042 (htslib: 1.13-3-gd16bed5) build(Jul 23 2021 22:06:16)
	-> /ricco/data/tyler/prog/bin/angsd -bam 
	-> Analysis helpbox/synopsis information:
	-> Fri Jul 23 22:25:05 2021
---------------
parseArgs_bambi.cpp: bam reader:
	-bam/-b		(null)	(list of BAM/CRAM files)
	-i		(null)	(Single BAM/CRAM file)
	-r		(null)	Supply a single region in commandline (see examples below)
	-rf		(null)	Supply multiple regions in a file (see examples below)
	-remove_bads	1	Discard 'bad' reads, (flag & 512) 
	-uniqueOnly	0	Discards reads that doesn't map uniquely
	-show		0	Mimic 'samtools mpileup' also supply -ref fasta for printing reference column
	-minMapQ	0	Discard reads with mapping quality below
	-minQ		13	Discard bases with base quality below
	-trim		0	Number of based to discard at both ends of the reads
	-trim		0	Number of based to discard at 5' ends of the reads
	-trim		0	Number of based to discard at 3' ends of the reads
	-only_proper_pairs 1	Only use reads where the mate could be mapped
	-C		0	adjust mapQ for excessive mismatches (as SAMtools), supply -ref
	-baq		0	adjust qscores around indels (1=normal baq 2= extended(as SAMtools)), supply -ref
	-redo-baq		0 (recompute baq, instead of using BQ tag)
	-setQscore	-1	Set qscore to this value, relevant for missing qscores
	-checkBamHeaders 1	Exit if difference in BAM headers
	-doCheck	1	Keep going even if datafile is not suffixed with .bam/.cram
	-downSample	0.000	Downsample to the fraction of original data
	-nReads		50	Number of reads to pop from each BAM/CRAMs
	-minChunkSize	250	Minimum size of chunk sent to analyses
	--ignore-RG	1	(dev only)
	+RG		(null)	Readgroups to include in analysis(can be filename)
	+LB		(null)	Libraries to include in analysis(can be filename)

Examples for region specification:
		chr:		Use entire chromosome: chr
		chr:start-	Use region from start to end of chr
		chr:-stop	Use region from beginning of chromosome: chr to stop
		chr:start-stop	Use region from start to stop from chromosome: chr
		chr:site	Use single site on chromosome: chr
```
</details>

<br>
Now lets calculate the genotype likelihoods. ANGSD provides the option to run it's analyses over a region
(see `$ANGSD -bam` for region specifications) and/or a subset of sites with `sites`. So, let's say we want to
only calculate genotype likelihoods for the region spanning sites 20,000 to 40,000 of chr7 at only the quality-controlled
sites that we generated yesterday (if you don't have the sites file a copy is at /ricco/data/tyler/output/qc_sites.pos).
Note that some options are on by default (e.g. -remove_bads), but we'll specify them to be explicit. We'll also specify some
commonly-used filters that are redundant with the QC we performed to make the `-sites` file, since they will demonstrate
how to these things if you want to be more stringent with your total site depth for example. Note that `-setMinDepth` and
`-setMaxDepth` require reads to be counted, i.e. `-doCounts 1`. You should look through each of the arguments used a make 
sure that you understand them. Information for ANGSD filters is [here](http://www.popgen.dk/angsd/index.php/Filters)
<br>
The GLs we calculate will use the SAMtools likelihood model `-GL 1` (see `$ANGSD -GL`) and output them as a text file
that we can easily examine.

```bash
$ANGSD -b $BAMLIST -ref $CICHREF -r chr7:20000-40000 -sites ~/ngs_intro/output/qc_sites.pos \
-remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -minQ 20 -minMapQ 20 -baq 1 -C 50 \
-setMinDepth 40 -setMaxDepth 1000 -doCounts 1 -GL 1 -doGlf 4 -out $DIR/output/calmas_region
```

The base alignment quality (BAQ) adjustment will recalibrate the base quality scores around INDELS. The `-C` option
adjusts mapping quality for reads with excessive mismatches. The [samtools](http://www.htslib.org/doc/samtools-mpileup.html)
documentation setting this to 50 for reads mapped with BWA. `-baq` and `-C` require that you supply the reference genome with `-ref`
<br>
Another useful option can be `-minInd X` in conjunction with `-minIndDepth Y`, which will require at least **X** individuals to be
covered by at least **Y** (default 1) reads to keep a site. We already applied this type of filter yesterday in a bit more of a sophisted way:
We required at least 5 individuals from *each* ecomorph to be covered by at least 2 reads in order for the site to pass QC. Note that in
ANGSD version 0.935-48-gff0c042 `-minIndDepth` is not recognized, but `-minInd` can still be used (the default required individual coverage of
1 will be used).

Looking at the standard output what percentage of the sites provided to ANGSD were actually retained?

<details>

<summary> click here for help </summary>

```bash
-> Total number of sites analyzed: 18433
-> Number of sites retained after filtering: 15754 
```
So ~85% of the sites were kept.

</details>

ANGSD always dumps a log file with information on how it was run. Check it out:

	less $DIR/output/calmas_region.arg

TEST
