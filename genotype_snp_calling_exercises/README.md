Estimation of allele frequencies, SNP calling, and genotype calling from NGS data
=================================================================================

Make directories and set environmental variables for this session

```bash
# set up directories

mkdir ~/ngs_analysis
mkdir ~/ngs_analysis/output

# set environment variables

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

<summary> click here for bam parsing options </summary>

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
only calculate genotype likelihoods for the region spanning sites 1 to 600,000 of chr7 at only the quality-controlled
sites that we generated yesterday (if you don't have the sites file a copy is at /ricco/data/tyler/output/qc_sites.pos).
Note that some options are on by default (e.g. -remove_bads), but we'll specify them to be explicit. We'll also specify some
commonly-used filters that are redundant with the QC we performed to make the `-sites` file, since they will demonstrate
how to these things if you want to be more stringent with your total site depth for example. Note that `-setMinDepth` and
`-setMaxDepth` require reads to be counted, i.e. `-doCounts 1`. You should look through each of the arguments used a make 
sure that you understand them. Information for ANGSD filters is [here](http://www.popgen.dk/angsd/index.php/Filters)
<br>
The GLs we calculate will use the SAMtools likelihood model `-GL 1` (see `$ANGSD -GL`) and output them as a text file
that we can easily examine. This will take ~2.5 minutes.

```bash
$ANGSD -b $BAMLIST -ref $CICHREF -r chr7:1-600000 -sites ~/ngs_intro/output/qc_sites.pos \
-remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -minQ 20 -minMapQ 20 -baq 1 -C 50 \
-setMinDepth 40 -setMaxDepth 700 -doCounts 1 -GL 1 -doGlf 4 -out $DIR/output/calmas_region
```

The base alignment quality (BAQ) adjustment will recalibrate the base quality scores around INDELS. The `-C` option
adjusts mapping quality for reads with excessive mismatches. The [samtools](http://www.htslib.org/doc/samtools-mpileup.html)
documentation setting this to 50 for reads mapped with BWA. `-baq` and `-C` require that you supply the reference genome with `-ref`
<br>
Another useful option can be `-minInd X` in conjunction with `-minIndDepth Y`, which will require at least **X** individuals to be
covered by at least **Y** (default 1) reads to keep a site. We already applied this type of filter yesterday in a bit more of a sophisted way:
We required at least 5 individuals from *each* ecomorph to be covered by at least 2 reads in order for the site to pass QC. Note that in
ANGSD version 0.935-48-gff0c042 `-minIndDepth` is not recognized, but `-minInd` can still be used (the default required individual depth of
1 will be used).

Looking at the standard output what percentage of the sites provided to ANGSD were actually retained?

<details>

<summary> click here for help </summary>

```bash
-> Total number of sites analyzed: 593192
-> Number of sites retained after filtering: 538150
```
So ~90% of the sites were kept.

</details>

ANGSD always dumps a log file with information on how it was run. Check it out:

	less $DIR/output/calmas_region.arg

Have a look at the GLs. The first two columns are the reference sequence (e.g. chromososome) and position. Then you have 10 likelihoods
for all possible genotypes in the order AA,AC,AG,AT,CC,CG,CT,GG,GT,TT. This set of 10 likelihoods is repeated sequentially starting from
the left of the file for each individual in the row order of individuals in the BAM file. The values are log-scaled likelihood ratios,
all scaled by the most likely genotype.

	less -S $DIR/output/calmas_region.glf.gz

We are analyzing 40 individuals, so we should have 402 fields in the glf file. You should confirm this and try to print the likelihoods
for the individual called CMASS6608007 at position chr7:1005 (each bam file is named by the individual, i.e. <idividual ID>.bam). What is 
their most likely genotype? If you need help you can click below.

<details>

<summary> click for help extracting GL info </summary>

```bash
# Count number of columns and subtract 2 (chromosome and position fields) to get the number of likelihood values

echo "$(($(zcat $DIR/output/calmas_region.glf.gz | head -n1 | wc -w)-2))"

# You should see that indeed there are 400 likelihood values.
# figure out what line CMASS6608007 is in the bam list.

INDNUM=$(grep -n "CMASS6608007.bam" $BAMLIST | cut -f1 -d':')
echo "$INDNUM"

# So this individual is at row 25 in the bam list. Now we can extract their likelihoods.

zcat $DIR/output/calmas_region.glf.gz | grep -m 1 $'^chr7\t10005\t' | cut -f 3- | perl -se '$start=($n-1)*10; @arr = split(/\t/,<>); print "@arr[$start .. $start+9]\n"' -- -n=$INDNUM
```
Since the likelihoods have been scaled to the most likely and log-transformed, the most likely genotype will
have a value of 0. For CMASS6608007 the 8th likelihood is zero, corresponding to the genotype 'GG'.

You could confirm this by looking at this site in the bcf file that you generated yesterday.

```bash
bcftools view -H -s CMASS6608007 -r chr7:10005 ~/ngs_intro/output/calmas_allsites.bcf.gz
```

Yep, it looks like this site is monomorphic with the reference a 'G' and this individual has data (your output should look like):

	chr7	10005	.	G	.	5018.12	.	DP=173;AD=173;SCR=21;MQSBZ=0.472176;FS=0;MQ0F=0;AN=2;DP4=96,77,0,0;MQ=57;NS=39	GT:DP:AD:SCR:QS	0/0:2:2:0:64

</details>

## allele frequency estimation

Now will estimate allele frequencies using the GLs we just calculated as input. Note that you can use the bams as input again,
but you'd have to recalculate the likelihoods (with `-GL` as before), which is redundant. If we do supply GLs as input we also
need to provide the number of individuals in the GL file, `-nInd`, and the reference index file, `-fai`.

We can get some information about how to estimate allele frequencies with `$ANGSD -doMaf`.

	-doMaf	0 (Calculate persite frequencies '.mafs.gz')
		1: Frequency (fixed major and minor)
		2: Frequency (fixed major unknown minor)
		4: Frequency from genotype probabilities
		8: AlleleCounts based method (known major minor)
		NB. Filedumping is supressed if value is negative

It's also useful to know how ANGSD can identify major and minor alleles, `$ANGSD -doMajorMinor`

	-doMajorMinor	0
	1: Infer major and minor from GL
	2: Infer major and minor from allele counts
	3: use major and minor from a file (requires -sites file.txt)
	4: Use reference allele as major (requires -ref)
	5: Use ancestral allele as major (requires -anc)
	6: Use EBD for major and minor (requires read data)
	7: Use EBD for minor (Major is ref) (requires -ref and read data, (similar to SAMtools))
	-rmTrans: remove transitions 0
	-skipTriallelic	0


We'll use `-doMajorMinor 1` to have ANGSD figure out what the major and minor alleles are from the GLs and then, based on these 
identified, alleles calculate their allele frequencies with `-doMaf 1`. If we wanted to account for more uncertainty in the 
identificaton of what the minor allele is we could use `-doMaf 2`. This latter approach would take a bit longer since there would
be three minor alleles to consider instead of just one. We'll also skip any sites that appear to have more than 2 alleles with
`-skipTriallelic 1` (though we should have filtered these out already with our qc_sites.pos file).

```bash
$ANGSD -glf10_text $DIR/output/calmas_region.glf.gz -nInd 40 -fai $CICHREF.fai \
-doMajorMinor 1 -doMaf 1 -skipTriallelic 1 -out $DIR/output/calmas_region_af
```

Check out the output: `less $DIR/output/calmas_region_af.mafs.gz` The columns are<br>
(1) chromosome<br>
(2) position<br>
(3) major allele<br>
(4) minor allele<br>
(5) maximum likelihood estimate of the minor allele frequency<br>
(6) the number of individuals with data<br>

Here's the first 5 sites:

	chromo	position	major	minor	knownEM	nInd
	chr7	115	A	C	0.000002	33
	chr7	116	A	C	0.000002	33
	chr7	117	A	C	0.000002	33
	chr7	118	A	C	0.000003	34
	chr7	119	G	A	0.000001	35

Now that you know how to extract allele frequencies, one interesting thing we could do is estimate the absolute divergence 
between the two ecomorphs of Astatotilapia calliptera in your data. Specifcally, we can use the allele frequencies
in the respective ecomorphs to calculate Dxy, which is the average number of pairwise nucleotide differences between them.
It's important for this calculation that for each site we estimate the allele frequency for the *same* allele in both ecomorphs.
In order to do this we can set which allele is the major allele using `-doMajorMinor 4`, which will assume that the reference allele
is major (which may not be true, but that's okay because we are just wanting to differentiate between alleles). Then for a biallelic 
site, the "minor" (or other) allele will be the same in both ecomorphs. Note that the minor allele is inferred from the GLs. We need 
to do this because the *actual* major, i.e. the most frequent allele, in the ecomorphs could be different. So, let's get started...
<br>
We need a maf file for each ecomorph, so let's start by generating this for the littoral morphs. We could simply split the glf file
containing all individuals into a file with only littoral GLs and another with only benthic GLs. Alternatively, we can calculate
the allele frequencies using two different bam lists as input, which is what we'll do here. The bam list for littoral individuals is
/ricco/data/tyler/littoral_bams.list and the bam list for benthics is /ricco/data/tyler/benthic_bams.list.

```bash
# estimate littoral morph allele frequencies (-ref is required for -doMajorMinor 4)

$ANGSD -b $DATDIR/littoral_bams.list -ref $CICHREF -r chr7:1-600000 -sites ~/ngs_intro/output/qc_sites.pos \
-remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -minQ 20 -minMapQ 20 -baq 1 -C 50 -skipTriallelic 1 \
-GL 1 -doMajorMinor 4 -doMaf 1 -out $DIR/output/calmas_region_af_littoral

# estimate benthic morph allele frequencies

$ANGSD -b $DATDIR/benthic_bams.list -ref $CICHREF -r chr7:1-600000 -sites ~/ngs_intro/output/qc_sites.pos \
-remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -minQ 20 -minMapQ 20 -baq 1 -C 50 -skipTriallelic 1 \
-GL 1 -doMajorMinor 4 -doMaf 1 -out $DIR/output/calmas_region_af_benthic
```
Take a look at these new maf files, you'll notice that there is an additional column specifying the reference allele since we included 
`-ref` in these runs. The major allele should match the ref allele since we used `doMajorMinor 4`.
<br>
Now we'll use these maf files as input to a program, dxyWindow, in order to calculate Dxy in XX kb windows. dxyWindow belongs to a 
suite of tools (https://github.com/tplinderoth/PopGenomicsTools) that can be used to perform various population genetic analyses.
You can enter `$DATDIR/prog/bin/dxyWindow` to see some help for the program:

	dxyWindow [options] <pop1 maf file> <pop2 maf file>
	
	Options:
	-winsize      INT     Window size in base pairs (0 for global calculation) [0]
	-stepsize     INT     Number of base pairs to progress window [0]
	-minind       INT     Minimum number of individuals in each population with data [1]
	-fixedsite    INT     (1) Use fixed number of sites from MAF input for each window (window sizes may vary) or (0) constant window size [0]
	-sizefile     FILE    Two-column TSV file with each row having (1) chromsome name (2) chromosome size in base pairs
	-skip_missing INT     Do not print windows with zero effective sites if INT=1 [0]
	
	Notes:
	* -winsize 1 -stepsize 1 calculates per site dxy
	* -sizefile is REQUIRED(!) with -fixedsite 0 (the default)
	* Both input MAF files need to have the same chromosomes in the same order
	* Assumes SNPs are biallelic across populations
	* For global Dxy calculations only columns 4, 5, and 6 below are printed
	* Input MAF files can contain all sites (including monomorphic sites) or just variable sites
	* -fixedsite 1 -winsize 500 would for example ensure that all windows contain 500 SNPs
	
	Output:
	(1) chromosome
	(2) Window start
	(3) Window end
	(4) dxy
	(5) number sites in MAF input that were analyzed
	(6) number of sites in MAF input that were skipped due to too few individuals

If you want to see what the `-sizefile` for the fAstCal1.2 reference genome looks like you can check it out: `cat $DATDIR/ref/fAstCal1.2_chr_lengths.txt`
Now, let's run dxyWindow.

```bash
$DATDIR/prog/bin/dxyWindow -winsize 10000 -stepsize 5000 -fixedsite 0 -sizefile $DATDIR/ref/fAstCal1.2_chr_lengths.txt -skip_missing 1 \
$DIR/output/calmas_region_af_littoral.mafs.gz $DIR/output/calmas_region_af_benthic.mafs.gz > $DIR/output/calmas_ecomorph_dxy.txt
``` 
The global Dxy printed to the screen states that there are on average ~758 nucleotide differences in 538,384 sites when comparing
littoral to benthic individuals. This means that there are approximately only 0.0014 pairwise differences on average per site, which 
is low and so indicates that these ecomorphs are likely quite genetically similar. This is only based on a small region and 
in practice you'd want to make this inference across the genome. Accordingly, looking at how Dxy is distributed across the genome, 
i.e. where the ecomorphs appear to be particularly divergent or similar could provide evidence for selection or introgression,
respectively.

You can see what Dxy looked like in 10 kb windows when sliding along our small example region in increments of 5 kb. Typically
it can be useful to visual these sliding window analyses as Manhattan plots.

	cat $DIR/output/calmas_ecomorph_dxy.txt

## call SNPs

At this point you're probably all amped up to actually call SNPs right? You're in luck because that's exactly what we're going to do.
You can statistically test for whether a site is variable in ANGSD using a likelihood ratio (LR) test, which compares the likelihood that the 
minor allele frequency (MAF) is zero (the null) to the likelihood of the estimated MAF (`-doMaf`). Under the null, -2log(LR statistic) 
is distributed according to a chi-square(1 d.f.), and so we can calculate a p-value for whether the estimated MAF is statistically 
different from zero, in which case the site is a SNP.
<br>
You can call SNPs based on a particular p-value cutoff using `SNP_pval`, so that's what we'll do for the first 1 MB on chromosome 7

```bash
$ANGSD -glf10_text $DIR/output/calmas_region.glf.gz -nInd 40 -fai $CICHREF.fai \
-doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -skipTriallelic 1 -out $DIR/output/calmas_region_snpcall
```
Look at the file `less $DIR/output/calmas_region_snpcall.mafs.gz`, which contains only sites called as variable.

The columns of this maf file are the same as before, except now there is an additional column 'pK-EM', which is the p-value corresponding
to the likelihood ratio test of whether a given site is variable.

Compare the number of called SNPs and distribution of allele frequencies to the case when you use a less stringent p-value cutoff
for whether a site is variable. This time we'll use a SNP p-value cutoff of 0.01.

$ANGSD -glf10_text $DIR/output/calmas_region.glf.gz -nInd 40 -fai $CICHREF.fai \
-doMajorMinor 1 -doMaf 1 -SNP_pval 0.01 -skipTriallelic 1 -out $DIR/output/calmas_region_snpcall_liberal

Count the number of SNPs in the two maf files of SNPs

``` bash
# p-value cutoff of 1e-6

echo "$(($(zcat $DIR/output/calmas_region_snpcall.mafs.gz | wc -l)-1))"

# p-value cutoff of 1e-2

echo "$(($(zcat $DIR/output/calmas_region_snpcall_liberal.mafs.gz | wc -l)-1))"
```
Plot the MAF distribution

```bash
$DATDIR/scripts/plotAFDist.R $DIR/output/calmas_region_snpcall.mafs.gz $DIR/output/calmas_region_snpcall_liberal.mafs.gz $DIR/output/snp_call_comparison

# Ignore the warnings about closing the unsed connections.

# View the plot

evince $DIR/output/snp_call_comparison.pdf
```
<details>

<summary> Click for plotAFDist.R code </summary>

```bash
#!/usr/bin/env Rscript

# plotAFDist.R <mafs.gz file 1> <mafz.gz file 2> <output prefix>

library(ggplot2)

# parse arguments
args <- commandArgs(trailingOnly=TRUE)

maf1 = read.table(gzfile(args[1],'rt'), head=TRUE)
maf2 = read.table(gzfile(args[2],'rt'), head=TRUE)
outprefix = args[3]

# combine the data so that allele frequency densities can be plotted together
maf1$SNP_pval <- "1e-6"
maf2$SNP_pval <- "1e-2"
maf.comb <- rbind(maf1, maf2)
maf.comb$SNP_pval <- factor(maf.comb$SNP_pval, levels=c("1e-6","1e-2"), order=TRUE)

# plot
pdf(file=paste0(outprefix,".pdf"))
ggplot(maf.comb, aes(knownEM, fill = SNP_pval)) + geom_density(alpha = 0.4, bw=0.015) + theme_classic(base_size=16) + theme(axis.line = element_line(size=0.5)) + xlab("MAF") + ylab("Density")
invisible(dev.off())
```
</details>

You can click below to view what you should have seen

<details>

<summary> snp calling comparison </summary>

`-SNP_pval 1e-6`: 1260 SNPs
`-SNP_pval 1e-2`: 1709 SNPs

![snp_calling_comparison](./outputs/snp_call_comparison.png)

</details>

Describe the difference between the two SNP p-value cutoffs.

## Genotype posteriors and calling	

## PCA

## SFS
-compare to expected SFS
-joint SFS?
