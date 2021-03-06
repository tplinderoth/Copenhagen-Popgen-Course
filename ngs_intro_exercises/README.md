INTRODUCTION TO NGS DATA EXERCISES
==================================

For these exercises we will use exome capture data from the Peruvian mimic poison dart frog, *Ranitomeya imitator*, which were 
sequenced to an average depth of 14x. The poison frog NGS data we mapped to a *de novo* exome assembly.
<br>

![imitator_morphs](../images/imitator.png)

<br>
We will also use whole genome sequencing data from the African cichlid fish *Astatotilapia callipeta*, which were sequenced to
a median depth of 5.7x. The cichlid NGS data were mapped to a chromosome-level *A. calliptera* assembly.
<br><br>

<img src="https://github.com/tplinderoth/Copenhagen-Popgen-Course/blob/main/images/calliptera.png" width="50%">

The goals for today are to:
* Become familiar with FASTQ format and know how to check the quality of raw sequencing reads
* Learn how to perform basic data quality control on sequencing reads
* Map sequencing reads to a reference genome
* Become familiar with SAM/BAM and pileup formats for mapped reads
* Learn how to use IGV to interactively explore mapping information
* Generate sequencing depth profiles
* Become familiar with VCF format and BCFtools
* Perform site-level quality control

Set some environmental variables
	
	# set up directories
	mkdir ~/ngs_intro
	mkdir ~/ngs_intro/output

	# set environment variables
	DIR=~/ngs_intro
	DATDIR=/ricco/data/tyler
	IGV=$DATDIR/prog/IGV_Linux_2.9.4/igv.sh
	BAMLIST=$DATDIR/cichlid_bams.list
	CICHREF=$DATDIR/ref/GCA_900246225.3_fAstCal1.2_genomic_chromnames_mt.fa
	BAMDIR=$DATDIR/bams
	SCRIPTS=$DATDIR/scripts
	BCFTOOLS=$DATDIR/prog/bin/bcftools
	ANGSD=/ricco/data/tyler/prog/bin/angsd

## fastq file format and quality

The first data we will take a look at is from exome capture data on poison dart frogs. These frogs 
have massive genomes and so WGS was not a practical approach for population genetics and gene mapping.

Have a look inside one of the frog fastq files

	less "$DATDIR/fastq/CH1401_R1.fastq"

visualize the data with [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). You can get some help information with `fastqc --help`.
FastQC will automatically detect file format type, but you can also specify it with `-f`.

	fastqc -outdir "$DIR/output/" $DATDIR/fastq/*.fastq
	ll "$DIR/output/"

You can scp the *.html files onto your local machine and view them in a browser. Here's an example of what the top of the web page would look like 
for the CH1401_R2.fastq file:

![fastq_preclean_example](./outputs/CH1401_R2_before.png)

Perhaps we have some adapter contamination to worry about so we can check this out without the browser interface

	unzip "$DIR/output/CH1401_R2_fastqc.zip" -d "$DIR/output/"
	
	# check out summary
	cat "$DIR/output/CH1401_R2_fastqc/summary.txt"

	# visualize the extent of adapter contamination (this can take a moment to show up, be patient...)
	eog "$DIR/output/CH1401_R2_fastqc/Images/adapter_content.png"

<details>

<summary> click here if you have trouble opening the fastqc image </summary>

![adapter_contamination](./outputs/CH1401_R2_preclean_adapter_content.png)

</details>

Can you see what FastQC is picking up on? How is the contamination distributed along the reads?

## Clean fastq files

We'll see if we can clean up the adapter contamination and low quality at the ends of the reads using [cutadapt](https://cutadapt.readthedocs.io/en/stable/). Normally for paired-end
data, you'd run in paired-end mode so that read pairs are maintained. Single-end mode can be run on the forward and reverse reads separately if it's desirable to retain unpaired reads 
(after one mate is removed due to quality issues). Cutadapt performs partial sums low-quality trimming. You can check out 
the [algorithm](https://cutadapt.readthedocs.io/en/stable/algorithms.html#quality-trimming-algorithm) to see how it works. We'll use a quality cutoff of `-q 15` and specify 
Illumina TruSeq [adapters](https://support.illumina.com/bulletins/2016/12/what-sequences-do-i-use-for-adapter-trimming.html), which is what were used for the frog sequencing.

	# RUN THIS ONE # paired-end mode #
	cutadapt -q 15 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 36 \
        -o "$DIR/output/CH1401_R1_clean.fastq.gz" -p "$DIR/output/CH1401_R2_clean.fastq.gz" \
	"$DATDIR/fastq/CH1401_R1.fastq" "$DATDIR/fastq/CH1401_R2.fastq"

	# single-end mode

	#cutadapt -q 15 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --minimum-length 36 \
        #-o "$DIR/output/CH1401_R1_clean_se.fastq.gz" "$DATDIR/fastq/CH1401_R1.fastq"

	#cutadapt3 -q 15 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 36 \
	#-o "$DIR/output/CH1401_R2_clean.fastq.gz" "$DATDIR/fastq/CH1401_R2.fastq"


use FastQC to check what the data quality looks like for CH1401_R2.fasq after trimming

	fastqc -outdir "$DIR/output/" "$DIR/output/CH1401_R2_clean.fastq.gz"

We can see that cutadapt did a good job of cleaning up the data:

![fastq_postclean_example](./outputs/CH1401_R2_after.png)
<br>
	
Note that the warning for the read length distribution is not of concern in this case. It's only alerting you that there are some short reads now, which we expect based on how we trimmed.

## Mapping

Use [bwa](https://github.com/lh3/bwa) to map the reads to the imitator exome assembly. 

	FROG_REF="$DATDIR/ref/imi_combined_targetedAndflanking_geneid.fasta"

	# normally you'd have to index the reference assembly (this has already been done for you, so skip the line below)
	# bwa index "$FROG_REF"

	# generated sorted bam
	bwa mem -R '@RG\tID:CH1401_capture1\tSM:CH1401' "$FROG_REF" "$DIR/output/CH1401_R1_clean.fastq.gz" "$DIR/output/CH1401_R2_clean.fastq.gz" \
	| samtools sort -O BAM > "$DIR/output/CH1401.bam"
	# consider option -I from bioanalyzer for bwa

	# index bam
	samtools index "$DIR/output/CH1401.bam"

## Working with mapped data

After you've mapped your cleaned-up fastq reads, you can have a look at the mapping information in SAM/BAM/CRAM files with samtools.
You can include the header with `-h`. When viewing the file type `/\@RG` and press `enter` to skip down to the bottom of the header to check that the
read group information that we intended to add is indeed there. Note also that you should see 'RG:Z:CH1401_capture1' associated with all of the reads
from this sample now as well.

	samtools view -h "$DIR/output/CH1401.bam" | less -S

The first 5 reads should look like:

	K00188:273:HG3VVBBXX:5:1102:21227:5077	99	contig1_combined	425	60	151M	=	522	248	GCCTCGGAGATGTGATGAATGTAACAGGTGCAGAGATTATTTCCAGGCCATGTGTGTGCGTCTGTGTGTAACAGCAAGAGGGGGAGAGAGAATGCAGAAAAGGAGAGAGCACATTGATGTCCCTGCTACGTCTCTGTAGCCTGTGAAACTT	AAFFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFFJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ	NM:i:0	MD:Z:151	MC:Z:151M	AS:i:151	XS:i:0	RG:Z:CH1401_capture1
	K00188:273:HG3VVBBXX:5:1102:21227:5077	147	contig1_combined	522	60	151M	=	425	-248	AAAAGGAGAGAGCACATTGATGTCCCTGCTACGTCTCTGTAGCCTGTGAAACTTACAGCAGATTGCCTTCGAGCGACCTAAATACTAAAGAAAAAAAGAGAAGAAAGATCGGAGAAGAGAAGAGGGTGGACACATAAAGGGTTATTTGCAT	JAAJJJFJJJJJFJJFAFJJJJFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJFJJJFJFJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJFJFFFAA	NM:i:0	MD:Z:151	MC:Z:151M	AS:i:151	XS:i:0	RG:Z:CH1401_capture1
	K00188:273:HG3VVBBXX:5:1101:22800:30556	163	contig1_combined	1850	60	151M	=	1879	180	TTTAATTTTTCAATCAGTCATGTGTAAATATGAATATAATTAAGATTACAATTTACATTTTCTTTCCACAGATCTAAACATTACTGCTGCTGTCCTGTCTCTTTTGAGTATAACTTTTATGGTAATGGGATCAATATGCATCACTATGGTT	AAFFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJAFJFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJAFJJJJJJJJJJJJJJJJJ	NM:i:5	MD:Z:0A4C4G20T8C110	MC:Z:151M	AS:i:130	XS:i:19	RG:Z:CH1401_capture1
	K00188:273:HG3VVBBXX:5:1101:22800:30556	83	contig1_combined	1879	60	151M	=	1850	-180	ATGAATATAATTAAGATTACAATTTACATTTTCTTTCCACAGATCTAAACATTACTGCTGCTGTCCTGTCTCTTTTGAGTATAACTTTTATGGTAATGGGATCAATATGCATCACTATGGTTCTCAGCAAAGGTGTGGAGTTCCTTCTGAA	JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFFFAA	NM:i:2	MD:Z:2T8C139	MC:Z:151M	AS:i:143	XS:i:0	RG:Z:CH1401_capture1
	K00188:273:HG3VVBBXX:5:1101:25621:48245	163	contig1_combined	2503	60	151M	=	2553	201	GGCTGGCTAATGGCAAGCAATCAGAACGTGAAGCTAGAATATGAGTATTCATGGTCCGTGGCCTGTGCGGCGGCAGCGGGAGGCGTCTTAATATTTGGTGGAATATGTTGCATTTTCTTGGTATTGCCGTCATTGCCCAAAAAGCCTTGGG	AAFFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJAJFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ<JJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJAJJJJJJJFJJJJJJJJJJJJJJJJJJAJJJ	NM:i:0	MD:Z:151	MC:Z:151M	AS:i:151	XS:i:0	RG:Z:CH1401_capture1

A useful resource for interpreting the bitwise flags is [here](https://broadinstitute.github.io/picard/explain-flags.html). What do they tell you about the mapping
of these first five reads?

Now, we are going to switch over to some WGS cichlid fish data and you'll see it's basically the same, but mapping to a more complete reference
will allow us to make some interesting inferences as you will see...

	# note the header won't show up because we omit -h
	samtools view "$BAMDIR/CMASS6169443.bam" | less -S

The first 5 lines, `samtools view $DATDIR/bams/CMASS6169443.bam | head -n 5` should look like

	HS30_18456:2:2102:20331:79677#7	99	chr7	73	57	125M	=	407	459	ACGAGCAAGAGGACAGTTTCAAAGCAGACATTATGGCACAATAAAAGCTGTGAGAGCAGAGATTCATTCATATCTTTACTATTTTTCTTTAGCTTCTACTCCCGCTGCTTACAATAAATCTGAGA	EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE	AS:i:125	XS:i:115	BC:Z:CAGATCTGTCTTTCCC	QT:Z:BBBBBFFF/BBBBFFF	MQ:i:40	MC:Z:125M	ms:i:4579	MD:Z:125	NM:i:0	RG:Z:18456_2#7
	HS30_18456:2:1213:20056:15084#7	163	chr7	75	0	125M	=	569	619	GAGCAAGAGGACAGTTTCAAAGCAGACATTATGGCACAATAAAAGCTGTAAGAGCAGAGATTCATTCATATCTTTACTATTTTTCTTTAGCTTCTACTCCCGCTGCTTACAATAAATCCGAGATA	DDDDDDDDDDDDDDD/AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA/AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA/DDDDDD/DDDD/FF/B5555F/B555FF	AS:i:115	XS:i:125	XA:Z:chr8,-23982502,125M,0;chr16,-2795516,125M,1;chr18,+31933308,125M,1;chr18,+30035828,125M,2;chr8,-23993408,125M,3;	MQ:i:9	MC:Z:125M	ms:i:4159	MD:Z:49G68T6	NM:i:2	RG:Z:18456_2#7
	HS30_18456:2:2211:1719:11546#7	163	chr7	103	51	125M	=	391	413	TTATGGCACAATAAAAGCTGTGAGAGCAGAGATTCATTCATATCTTTACTATTTTTCTTTAGCTTCTACTCCCGCTGCTTACAATAAATCTGAGATAAAGAACAAACAAGAACAACTCACATGGC	EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE	AS:i:125	XS:i:115	MQ:i:40	MC:Z:125M	ms:i:4597	MD:Z:125	NM:i:0	RG:Z:18456_2#7
	HS30_18456:2:2101:10435:101060#7	99	chr7	118	0	125M	=	499	506	AGCTGTAAGAGCAGAGATTCATTCATATCTTTACTATTTTTCTTTAGCTTCTACTCCCGCTGCTTACAATAAATCCGAGATAAAGAACAAACAAGAACAACTCACATGGCATTGATTGTTTAGTT	>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>/F	AS:i:115	XS:i:125	XA:Z:chr8,-23982459,125M,0;chr18,+30035871,125M,0;chr16,-2795473,125M,1;chr18,+31933351,125M,1;U_scaffold_23,-9073,125M,2;	BC:Z:CAGATCTGTCTTTCCC	QT:Z:BBBBBFFF/7<BB/B/	MQ:i:8	MC:Z:125M	ms:i:4556	MD:Z:6G68T49	NM:i:2	RG:Z:18456_2#7
	HS30_18456:2:1316:11207:16667#7	163	chr7	126	0	125M	=	459	458	GAGCAGAGATTCATTCATATCTTTACTATTTTTCTTTAGCTTCTACTCCCGCTGCTTACAATAAATCCGAGATAAAGAACAAACAAGAACAACTCACATGGCATTGATTGTTTAGTTCAGTGTCA	DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD/	AS:i:120	XS:i:125	MQ:i:6	MC:Z:125M	ms:i:4546	MD:Z:67T57	NM:i:1	RG:Z:18456_2#7

You can quickly view all reads that map to a specific region by specifying the region-of-interest after the bam.
For example, to see all reads that map to chr7:1000-1010 you'd use

```bash
samtools view $DATDIR/bams/CMASS6169443.bam chr7:1000-1010 | less -S
``` 

## SAMtools mpileup

Lets looks at the mapped data on chr7:10,000-10,015

	# check out the bam list
	cat "$BAMLIST"

	samtools mpileup -b $BAMLIST -f $CICHREF -r chr7:10000-10015 | less -S

Here's what the data looks like for the first 6 individuals in the bam file, `samtools mpileup -b $BAMLIST -f $CICHREF -r chr7:10000-10015 | cut -f1-21 | column -t`

	chr7  10000  T  2  .,      EB    7  .,,,,,^],   >IIBGGE   6  .....,       DEGEGG     9   ,.,,.,,.,     DABGIIIII   5  ,,.,.  >AB/A  6  .,,...  FGGBGB
	chr7  10001  C  3  .,^].   EBB   7  .,,,,,,     >IIBGGI   6  .....,       DEGEGG     9   ,.,,.,,.,     DABGIIIII   5  ,,.,.  >ABDA  6  .,,...  3GGBGB
	chr7  10002  G  3  .,.     EBB   7  .,,,,,,     >IIBGGI   6  .....,       DEGEGG     9   ,.,,.,,.,     DABGIIIII   5  ,,.,.  >ABDA  6  .,,...  3GGBGB
	chr7  10003  G  3  .,.     EBB   7  .,,,,,,     >IIBGGI   6  .....,       DEGEGG     10  ,.,,.,,.,^].  DABGIIIIIE  5  ,,.,.  >AB/A  6  .,,...  3GGBGB
	chr7  10004  A  3  .,.     EBB   8  .,,,,,,^].  >IIBGGIE  6  .....,       DEGEGG     10  ,.,,.,,.,.    DABGIIIIII  5  ,,.,.  >AB/A  6  .,,...  3GGBGB
	chr7  10005  G  3  .,.     EBB   8  .,,,,,,.    >IIBGGIG  6  .....,       DEGEGG     10  ,.,,.,,.,.    DABGIIIIII  5  ,,.,.  >ABDA  6  .,,...  BGGBGB
	chr7  10006  A  3  .,.     EBB   8  .,,,,,,.    >IIBGGIG  6  .....,       DEGEGG     10  ,.,,.,,.,.    DABGIIIIII  5  ,,.,.  >ABDA  6  C,,...  5G/BGB
	chr7  10007  G  3  .,.     EBB   8  .,,,,,,.    >IIBGGIG  7  .....,^].    DEGEGGE    10  ,.,,.,,.,.    DABGIIIIII  5  ,,.,.  >ABDA  6  .,,...  5GIBGB
	chr7  10008  C  3  .,.     EBB   8  .,,,,,,.    >IIBGGIG  7  .....,.      DEGEGGG    10  ,.,,.,,.,.    DABGIIIIII  5  ,,.,.  >ABDA  6  .,,...  5GIBGB
	chr7  10009  A  3  .,.     EBB   8  .,,,,,,.    >IIBGGIG  8  .....,.^],   DEGEGGGB   10  ,.,,.,,.,.    DABGIIIIII  5  ,,.,.  >ABDA  6  .,,...  5GIBGB
	chr7  10010  G  3  .,.     EBB   8  .,,,,,,.    >IIBGGIG  8  .....,.,     DEGEGGGB   10  ,.,,.,,.,.    DABGIIIIII  5  ,,.,.  >AB/A  6  .,,...  DGIBGB
	chr7  10011  C  3  .,.     EBB   8  .,,,,,,.    >IIBGGIG  8  .....,.,     DEGEGGGB   10  ,.,,.,,.,.    DABGIIIIII  5  ,,.,.  >ABAA  6  .,,...  DGIBGB
	chr7  10012  T  3  .,.     EBB   8  .,,,,,,.    >IIBGGIG  8  .....,.,     DEGEGGGB   10  ,.,,.,,.,.    DABGIIIIII  5  ,,.,.  >/BAA  6  .,,...  DGIBGB
	chr7  10013  T  3  .,.     EBB   8  .,,,,,,.    >IIBGGIG  9  .....,.,^],  DEGEGGGBE  10  ,.,,.,,.,.    DABGIIIIII  5  ,,.,.  >BBAA  6  .,,...  DGIBGB
	chr7  10014  A  4  .,.^].  EBBE  8  .,,,,,,.    >IIBGGIG  9  .G...gG,,    DEGEGGGBG  10  ,.,,.,,.,.    D3BGIIIIII  5  ,,.,C  >BBAA  6  .,,...  DGIBGB
	chr7  10015  G  4  .,..    EBBI  8  .,,,,,,.    >IIBGGIG  9  .....,.,,    DEGEGGGBG  10  ,.,,.,,.,.    D3BGIIIIII  5  ,,.,.  >BBAA  6  .,,...  /GIBGB

Looking at this data do you think you could confidently call genotypes?

Now figure out how to generate the same pileup but one that includes mapping quality and that only consider reads with a minimum quality/BAQ of 20.
Write the output to a file called 'mapq_example.pileup' in the 'output' ${DIR}/output/ <br>

You can see help for samtools mpileup with
	samtools mpileup --help

<details>

<summary>Click here for code</summary>

``` bash
samtools mpileup -b $BAMLIST -f $CICHREF -r chr7:10000-10015 -Q 20 --output-MQ > $DIR/output/mapq_example.pileup

# examine the file in an easy-to-read way (aligns the columns nicely)
cat $DIR/output/mapq_example.pileup | column -t | less -S
```

</details>

## IGV

A useful way to visualize mapping information is with the Integrative Genomics Viewer (IGV)

	# start up igv
	$IGV

	# SKIP THE FOLLOWING (it's already been done for you)
	# load reference genome
	# In top menu: 'Genomes' -> 'Load Genome from File...' and select /ricco/data/tyer/ref/GCA_900246225.3_fAstCal1.2_genomic_chromnames_mt.fa
	# load BAM
	# In top menu: 'File' -> 'Load from File...' select /ricco/data/tyler/bams/CMASS6607991.bam
	
	# START HERE
	# 'File -> Open session... -> '+ Other locations/Computer/ricco/data/tyler/igv_sessions/CMASS6607991.xml' 
	# In the box next to 'Go' type 'chr7:18,078,500-18,102,000' and press 'Go'
	# Right-click on the track containing the reads in the IGV window and select 'view as pairs' and 'Group alignment by' -> 'pair orientation'
	# click on one of the green reads and figure out where its mate maps. What orientation are these reads mapping in?
	# You can use the scroll bar along the top to explore more mapping along chr7.

<details>

<summary> Click here to see what you should be seeing </summary>

![bam_igv_region](./outputs/cichlid_igv_region.png)

</details>

<br> <br>

## coverage plot

You can also see that some samples do not have unproperly mapped reads, e.g. CMASS6169461.

![igv_region_comparison](./outputs/cichlid_comparison_igv_region.png)

Let's check the coverage spanning the region where reads were mapping in the RL orientation +/- ~20 kb. Do this for a subset of samples
showing different mapping patterns (4 strange, 2 normal controls).

``` bash
BAM1=$BAMDIR/CMASS6608026.bam # improper mapping
BAM2=$BAMDIR/CMASS6607991.bam # improper mapping
BAM3=$BAMDIR/CMASS6389722.bam # improper mapping
BAM4=$BAMDIR/CMASS6169445.bam # improper mapping
BAM5=$BAMDIR/CMASS6169461.bam # proper mapping
BAM6=$BAMDIR/CMASS6169443.bam # proper mapping

(( echo -e "CHR\tPOS\t$BAM1\t$BAM2\t$BAM3\t$BAM4\t$BAM5\t$BAM6" | sed "s;\($BAMDIR\/\|.bam\);;g" ); \
(samtools depth -a -r chr7:18059155-18120834 $BAM1 $BAM2 $BAM3 $BAM4 $BAM5 $BAM6)) > $DIR/output/cichlid_region_depth.txt

# calculate relative depth and plot a depth profile
# calmas_meta_sub.txt file contains cichlid metadata, including average genome-wide sequence depth

$SCRIPTS/plot_depth_region.R $DIR/output/cichlid_region_depth.txt $DATDIR/calmas_meta_sub.txt $DIR/output/region_depth_profile
```
<details>

<summary> Click here to see plot_depth_region.R code </summary>

```bash
#!/usr/bin/env Rscript

# plot_depth_region.R <depth file> <metadata file> <output prefix>

## parse input from the command line

args <- commandArgs(trailingOnly=TRUE)

depth <- read.table(args[1], head=TRUE)
meta <- read.table(args[2], head=TRUE)
outprefix <- args[3]

## normalize the depth of every sample by their genome-wide depth
for (i in 3:ncol(depth)) { depth[,i] = depth[,i]/meta$sequence.depth[which(meta$sample.ID == colnames(depth)[i])]}

# print relative depth of samples for the region spanning improperly mapped reads
start = 18079155 # start of region where reads map in RL orientation
end = 18100834 # end of region where reads map in RL orientation
cat("\nRelative sequencing depth\n\n")
for (i in 3:ncol(depth)) { cat(paste0(colnames(depth)[i],": ", round(mean(depth[,i][which(depth$POS > start & depth$POS < end)]),digits=3),"\n")) }

## take mean in 100 bp windows for smoothing (makes viewing the plots easier) ##
inc=100
win=matrix(nrow=ceiling(nrow(depth)/inc), ncol=ncol(depth))
a=1
b=inc
for (i in 1:nrow(win)) {
        win[i,1] = depth$POS[a]
        win[i,2] = depth$POS[b]
        for (j in 3:ncol(depth)) { win[i,j] = mean(depth[a:b,j])}
        a=b+1
        b=a+inc-1
        if (b > nrow(depth)) b = nrow(depth)
}

# plot relative coverage profile for 2 of the individuals

n=nrow(win)
mid = round((win[,1]+win[,2])/2)
pdf(file=paste0(outprefix,".pdf"),width=8,height=5)
plot(x=1:n, y=win[,3],type="l", ylim=c(0,4), col="royalblue4", xlab="Window midpoint", ylab="Relative sequence depth", cex.axis=1.2, cex.lab=1.2, lwd=1.5, xaxt='n')
axis(1, at=c(0,100,200,300,400,500,600),labels=mid[c(1,100,200,300,400,500,600)],cex.axis=1.2)
lines(x=1:n, y=win[,7], col="violet", lwd=1.5)
legend("topleft", c("CMASS6608026","CMASS6169461"),col=c("royalblue4","violet"), lwd=3, bty='n')
invisible(dev.off())
```

</details>

Lets look at the depth plot

```bash
evince "$DIR/output/region_depth_profile.pdf"
```
<details>

<summary> click here if you having trouble viewing the plot </summary>

![region_depth_profile](./outputs/region_depth_profile.png)

</details>

What can you infer about the genomes of these six cichlids samples from the mapping information
that you've seen in IGV and depth profiling?

<details>

<summary> Click here for answer </summary>
There is a tandem duplication on chr7:18,079,155-18,100,834 encompassing the *gsdf* gene. The tandem-duplicate allele masculinizes fish and so operates as a sex determiner. <br>
~1x relative depth = Homozygous for the non-duplicated *gsdf* allele (2 *gsdf* copies)<br>
~1.5x relative depth = Heterozygous for the duplication allele (3 *gsdf* copies)<br>
~2x = Homozygous for the duplication allele (4 *gsdf* copies)

</details>

<br>

## bcftools filtering

Generate a VCF with some info that we can filter on: the first 1 MB of chr7 (this will take around 1 minute).

```bash
$BCFTOOLS mpileup -f $CICHREF -b $BAMLIST \
-d 40000 -L 40000 -r chr7:1-1000000 -q 13 -Q 13 --ff UNMAP,SECONDARY,QCFAIL,DUP -a FORMAT/AD,FORMAT/DP,QS,FORMAT/SCR,INFO/AD,INFO/SCR -p O u \
| $BCFTOOLS call --ploidy 2 -a PV4,GQ,GP -m -P 0.001 -O u | $BCFTOOLS +fill-tags -O b -o $DIR/output/calmas_allsites.bcf.gz -- -t'AF,NS,ExcHet'

# index the bcf for rapid access to specific regions with 'bcftools view -r ...'
# not necessary now, but demonstrating how this is done

bcftools index $DIR/output/calmas_allsites.bcf.gz
```
The vcf is in compressed binary fomat to save space so we'll have to use bcftools to view it. Let's have a look at the annotations that we
can use to filter on.

	$BCFTOOLS view $DIR/output/calmas_allsites.bcf.gz | less -S

Lets see how we could extract some information for quantities we might want to examine the distribution of prior to filtering.

``` bash
((echo -e "CHROM\tPOS\tDP\tMQ\tSTRAND_BIAS\tBASEQ_BIAS\tMQ_BIAS\tPOS_BIAS\tEXCHET"); \
($BCFTOOLS query -f "%CHROM\t%POS\t%INFO/DP\t%INFO/MQ\t%INFO/PV4{0}\t%INFO/PV4{1}\t%INFO/PV4{2}\t%INFO/PV4{3}\t%ExcHet\n" $DIR/output/calmas_allsites.bcf.gz)) > $DIR/output/allsites_stats.txt

# Have a look ('.' indicates that the value is missing)
less $DIR/output/allsites_stats.txt

# let's plot these values and examine percentiles to get an idea of extremes
$SCRIPTS/plotStatDist.R $DIR/output/allsites_stats.txt $DIR/output/allsites_stats_plot

evince $DIR/output/allsites_stats_plot.pdf
```

What are these different quantities measuring? The info can be found in the VCF header.

<details>

<summary> Click here for the plotStatDist.R code </summary>

```bash
#!/usr/bin/env Rscript

# plotStatDist.R <file of vcf stats> <outfile prefix>

# parse arguments
args <- commandArgs(trailingOnly=TRUE)
df <- read.table(args[1],head=TRUE, na.strings=".")
suppressWarnings(df$EXCHET <- as.numeric(as.character(df$EXCHET))) # need this because R interprets this as a 'factor' due to many 'NA' at beginning
outprefix <- args[2]

# calculate percentiles
percdf=matrix(nrow=ncol(df)-2, ncol=13)
i=1
for (j in 3:ncol(df)) { percdf[i,] = quantile(df[,j], c(0.001, 0.01, seq(from=0.1, to=0.9, by=0.1), 0.99, 0.999), na.rm=TRUE); i=i+1 }
colnames(percdf) = names(quantile(df$DP, c(0.001, 0.01, seq(from=0.1, to=0.9, by=0.1), 0.99, 0.999), na.rm=TRUE))
rownames(percdf) = colnames(df[3:ncol(df)])

# print quantiles
cat("\n\nPercentiles\n\n")
print(percdf); cat("\n")

# plot stats
pdf(file=paste0(outprefix,".pdf"))
par(mfrow=c(2,2))
hist(df$DP, breaks=20, ylab="Number sites", xlab="Total site depth", main="DP", cex.lab=1.2)
hist(df$MQ, breaks=20, ylab="Number sites", xlab="Mapping quality", main="MQ", cex.lab=1.2)
hist(df$STRAND_BIAS, breaks=20, ylab="Number sites", xlab="Strand bias p-value", main="PV[0]", cex.lab=1.2)
hist(df$BASEQ_BIAS, breaks=20, ylab="Number sites", xlab="Base quality p-value", main="PV[1]", cex.lab=1.2)
hist(df$BASEQ_BIAS, breaks=10000, ylab="Number sites", xlab="Base quality p-value", main="PV[1] Zoom", cex.lab=1.2, xlim=c(0,0.001))
hist(df$MQ_BIAS, breaks=20, ylab="Number sites", xlab="Mapping quality bias p-value", main="PV[2]", cex.lab=1.2)
hist(df$POS_BIAS, breaks=20, ylab="Number sites", xlab="tail distance bias p-value", main="PV[3]", cex.lab=1.2)
hist(df$EXCHET, breaks=20, ylab="Number sites", xlab="Excess heterozygosity p-value", main="ExcHet", cex.lab=1.2)
invisible(dev.off())
```

</details>
<br>
<details>

<summary> click here to see quality statistic percentiles and plots </summary>

	Percentiles
	
	                    0.1%           1%          10%          20%          30%
	DP          1.000000e+00 1.200000e+01 1.060000e+02 1.570000e+02 1.780000e+02
	MQ          1.300000e+01 2.200000e+01 3.500000e+01 4.300000e+01 5.400000e+01
	STRAND_BIAS 4.761011e-07 1.219852e-02 3.061478e-01 4.393060e-01 4.701490e-01
	BASEQ_BIAS  5.224217e-31 6.404809e-21 1.103558e-12 3.795092e-10 1.480062e-08
	MQ_BIAS     0.000000e+00 2.802600e-45 1.185121e-01 1.000000e+00 1.000000e+00
	POS_BIAS    1.145862e-04 2.673683e-03 3.093082e-02 1.212964e-01 2.907322e-01
	EXCHET      1.022730e-11 5.025090e-06 2.091916e-01 4.348290e-01 5.858350e-01
	                    40%         50%          60%         70%          80%
	DP          1.90000e+02 2.00000e+02 2.090000e+02 2.19000e+02 2.300000e+02
	MQ          5.90000e+01 6.00000e+01 6.000000e+01 6.00000e+01 6.000000e+01
	STRAND_BIAS 4.91736e-01 7.68149e-01 1.000000e+00 1.00000e+00 1.000000e+00
	BASEQ_BIAS  2.56738e-07 2.99443e-06 3.206144e-05 3.75725e-04 9.199838e-03
	MQ_BIAS     1.00000e+00 1.00000e+00 1.000000e+00 1.00000e+00 1.000000e+00
	POS_BIAS    1.00000e+00 1.00000e+00 1.000000e+00 1.00000e+00 1.000000e+00
	EXCHET      7.00602e-01 7.97872e-01 8.622170e-01 9.16105e-01 9.589040e-01
	                   90% 99%   99.9%
	DP          246.000000 461 954.804
	MQ           60.000000  60  60.000
	STRAND_BIAS   1.000000   1   1.000
	BASEQ_BIAS    0.254081   1   1.000
	MQ_BIAS       1.000000   1   1.000
	POS_BIAS      1.000000   1   1.000
	EXCHET        0.987342   1   1.000

![quality_stat_distribution_1](./outputs/allsites_stat_plot_1.png)
<br>
![quality_stat_distribution_2](./outputs/allsites_stat_plot_2.png)

</details>

Filter the VCF. We'll avoid dumping another VCF with just sites that pass our quality controls by extracting just the sites.
We can use these sites with ANGSD tomorrow.

```bash
$BCFTOOLS norm -f $CICHREF -m +any $DIR/output/calmas_allsites.bcf.gz \
| $BCFTOOLS view -i 'N_PASS(FMT/DP[0-14] > 2) > 5 && N_PASS(FMT/DP[15-39] > 2) > 5' $DIR/output/calmas_allsites.bcf.gz \
| $BCFTOOLS view -e 'INDEL=1 || INFO/MQ < 25 || INFO/DP > 700 || INFO/PV4[0] < 1.2e-02 || INFO/PV4[1] < 6.4e-21 || INFO/PV4[2] < 2.8e-45 || INFO/PV4[3] < 2.7e-3 || INFO/ExcHet < 5e-06' -M 2 \
| $BCFTOOLS query -f "%CHROM\t%POS\n" > $DIR/output/qc_sites.pos

# Lets have a look at our list of quality-controlled sites
less $DIR/output/qc_sites.pos

# If we want to use these sites in ANGSD (tomorrow), we need to index them
$ANGSD sites index $DIR/output/qc_sites.pos
```

Can you describe what each of the filters was doing?
What pecentage of site was removed (discounting those with entirely missing data)?

<details>

<summary> Click here for help </summary>

```bash
$BCFTOOLS query -f "%POS\n" $DIR/output/calmas_allsites.bcf.gz | wc -l
# 995197 sites had some data initially

wc -l $DIR/output/qc_sites.pos
# 900037 sites passed quality-control
```

Therefore, we removed 95160 sites, which represents ~10% of the sites with data. This seems reasonable.

</details>

You could bypass dumping the initial VCF entirely using one large string of pipes. Give it a try.
