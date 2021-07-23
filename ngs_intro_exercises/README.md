INTRODUCTION TO NGS DATA EXERCISES
==================================

Set some environmental variables

	# from your home directory
	
	mkdir ngs_intro
	DIR=~/ngs_intro
	mkdir $DIR/output
	DATDIR=/ricco/data/tyler
	IGV=$DATDIR/prog/IGV_Linux_2.9.4/igv.sh
	BAMLIST=$DATDIR/cichlid_bams.list
	CICHREF=$DATDIR/ref/GCA_900246225.3_fAstCal1.2_genomic_chromnames_mt.fa
	BAMDIR=$DATDIR/bams
	SCRIPTS=$DATDIR/scripts

## fastq

Have a look inside one of the frog fastq files

	less "$DATDIR/fastq/CH1401_R1.fastq"

visualize data with [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). You can get some help info with `fastqc --help`
Will automatically detect file format type, but can specify with `-f`.

	fastqc -outdir "$DIR/output/" $DATDIR/fastq/*.fastq
	ll "$DIR/output/"

You can scp the html files onto your local machine and view them in the web. It would look something like:

![fastq_preclean_example](./outputs/CH1401_R2_before.png)

Perhaps we have some adapter contamination to worry about so we can check this out without the browser interface

	unzip "$DIR/output/CH1401_R2_fastqc.zip" -d "$DIR/output/"
	
	# check out summary
	cat "$DIR/output/CH1401_R2_fastqc/summary.txt"

	# visualize the extent of adapter contamination (this can take a moment to show up, be patient...)
	eog "$DIR/output/CH1401_R2_fastqc/Images/adapter_content.png"

Can you see what fastq is picking up on? How is the contamination distributed along the read?

## clean fastq

We'll see if we can clean up the adapter contamination and low quality at the ends of the read using cutadapt. Normally run in paired-end mode so that
read pairs are retained. Single-end mode can be run on the forward and reverse reads separately if it's desirable to retain unpaired reads. We'll run on
just the reverse read that we examined for sake of time. Partial sums low-quality trimming using -q 15 (algorithm: https://cutadapt.readthedocs.io/en/stable/algorithms.html#quality-trimming-algorithm)
Specifing Illumina TruSeq [adapters](https://support.illumina.com/bulletins/2016/12/what-sequences-do-i-use-for-adapter-trimming.html).

	# paired-end mode
	cutadapt -q 15 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 36 \
        -o "$DIR/output/CH1401_R1_clean.fastq.gz" -p "$DIR/output/CH1401_R2_clean.fastq.gz" \
	"$DATDIR/fastq/CH1401_R1.fastq" "$DATDIR/fastq/CH1401_R2.fastq"

	# single-end mode

	#cutadapt -q 15 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --minimum-length 36 \
        #-o "$DIR/output/CH1401_R1_clean_se.fastq.gz" "$DATDIR/fastq/CH1401_R1.fastq"

	#cutadapt3 -q 15 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 36 \
	#-o "$DIR/output/CH1401_R2_clean.fastq.gz" "$DATDIR/fastq/CH1401_R2.fastq"


use fastqc to check what the data quality looks like for CH1401_R2.fasq after trimming

	fastqc -outdir "$DIR/output/" "$DIR/output/CH1401_R2_clean.fastq.gz"

We can see that cutadapt did a good job of cleaning up the data

![fastq_postclean_example](./outputs/CH1401_R2_after.png)

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

Now, we are going to switch over to some WGS cihlid data and you'll see it's basically the same, but mapping to a more complete reference
will allow us to make some interesting inferences as you will see...

	# note the header won't show up because we omit -h
	samtools view "$BAMDIR/CMASS6169443.bam" | less -S

The first 5 lines, `samtools view "$DIR/data/bams/CMASS6169443.bam | head -n 5` should look like

	HS30_18456:2:2102:20331:79677#7	99	chr7	73	57	125M	=	407	459	ACGAGCAAGAGGACAGTTTCAAAGCAGACATTATGGCACAATAAAAGCTGTGAGAGCAGAGATTCATTCATATCTTTACTATTTTTCTTTAGCTTCTACTCCCGCTGCTTACAATAAATCTGAGA	EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE	AS:i:125	XS:i:115	BC:Z:CAGATCTGTCTTTCCC	QT:Z:BBBBBFFF/BBBBFFF	MQ:i:40	MC:Z:125M	ms:i:4579	MD:Z:125	NM:i:0	RG:Z:18456_2#7
	HS30_18456:2:1213:20056:15084#7	163	chr7	75	0	125M	=	569	619	GAGCAAGAGGACAGTTTCAAAGCAGACATTATGGCACAATAAAAGCTGTAAGAGCAGAGATTCATTCATATCTTTACTATTTTTCTTTAGCTTCTACTCCCGCTGCTTACAATAAATCCGAGATA	DDDDDDDDDDDDDDD/AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA/AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA/DDDDDD/DDDD/FF/B5555F/B555FF	AS:i:115	XS:i:125	XA:Z:chr8,-23982502,125M,0;chr16,-2795516,125M,1;chr18,+31933308,125M,1;chr18,+30035828,125M,2;chr8,-23993408,125M,3;	MQ:i:9	MC:Z:125M	ms:i:4159	MD:Z:49G68T6	NM:i:2	RG:Z:18456_2#7
	HS30_18456:2:2211:1719:11546#7	163	chr7	103	51	125M	=	391	413	TTATGGCACAATAAAAGCTGTGAGAGCAGAGATTCATTCATATCTTTACTATTTTTCTTTAGCTTCTACTCCCGCTGCTTACAATAAATCTGAGATAAAGAACAAACAAGAACAACTCACATGGC	EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE	AS:i:125	XS:i:115	MQ:i:40	MC:Z:125M	ms:i:4597	MD:Z:125	NM:i:0	RG:Z:18456_2#7
	HS30_18456:2:2101:10435:101060#7	99	chr7	118	0	125M	=	499	506	AGCTGTAAGAGCAGAGATTCATTCATATCTTTACTATTTTTCTTTAGCTTCTACTCCCGCTGCTTACAATAAATCCGAGATAAAGAACAAACAAGAACAACTCACATGGCATTGATTGTTTAGTT	>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>/F	AS:i:115	XS:i:125	XA:Z:chr8,-23982459,125M,0;chr18,+30035871,125M,0;chr16,-2795473,125M,1;chr18,+31933351,125M,1;U_scaffold_23,-9073,125M,2;	BC:Z:CAGATCTGTCTTTCCC	QT:Z:BBBBBFFF/7<BB/B/	MQ:i:8	MC:Z:125M	ms:i:4556	MD:Z:6G68T49	NM:i:2	RG:Z:18456_2#7
	HS30_18456:2:1316:11207:16667#7	163	chr7	126	0	125M	=	459	458	GAGCAGAGATTCATTCATATCTTTACTATTTTTCTTTAGCTTCTACTCCCGCTGCTTACAATAAATCCGAGATAAAGAACAAACAAGAACAACTCACATGGCATTGATTGTTTAGTTCAGTGTCA	DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD/	AS:i:120	XS:i:125	MQ:i:6	MC:Z:125M	ms:i:4546	MD:Z:67T57	NM:i:1	RG:Z:18456_2#7

## samtools mpileup

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

	# SKIP THE FOLLOWING
	# load reference genome
	# In top menu: 'Genomes' -> 'Load Genome from File...' and select XXX/GCA_900246225.3_fAstCal1.2_genomic_chromnames_mt.fa
	# load BAM
	# In top menu: 'File' -> 'Load from File...' select XXX/CMASS6607991.bam
	

	# 'File -> Open session... -> /ricco/data/tyler/igv_sessions/CMASS6607991.xml' 
	# In the box next to 'Go' type 'chr7:18,078,500-18,102,000' and press 'Go'
	# Right-click on the track containing the reads in the IGV window and select 'view as pairs' and 'Group alignment by' -> 'pair orientation'
	# click on one of the green reads and figure out where its mate maps

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
What can you infer about the genomes of these six cichlids samples from the mapping information
that you've seen in IGV and depth profiling?

<details>

<summary> Click here for answer </summary>
Tandem duplication of gsdf gene, which generates males. <br>
~1x relative depth = non-duplicated <br>
~1.5x = heterozygous for duplication <br>
~2x = homozygous for duplication

</details>

<br>

## bcftools filtering

First let's generate a distribution of the total site depth

Generate a VCF with some info that we can filter on

```bash
bcftools mpileup -f $CICHREF -b $BAMLIST \
-d 40000 -L 40000 -r chr7:1-50000 -q 13 -Q 13 --ff UNMAP,SECONDARY,QCFAIL,DUP -a FORMAT/AD,FORMAT/DP,QS,FORMAT/SCR,INFO/AD,INFO/SCR -p O u \
| bcftools call --ploidy 2 -a PV4,GQ,GP -m -P 0.001 -O u | bcftools +fill-tags -O b -o $DIR/output/calmas_allsites.bcf.gz -- -t'AF,NS,ExcHet'
```
The vcf is in compressed binary fomat to save space so we'll have to use bcftools to view it. Let's have a look at the annotations that we
can use to filter on.

	bcftools view $DIR/output/calmas_allsites.bcf.gz | less -S

Lets see how we could extract some information for quantities we might want to examine the distribution of prior to filtering.

``` bash
((echo -e "CHROM\tPOS\tDP\tMQ\tSTRAND_BIAS\tBASEQ_BIAS\tMQ_BIAS\tPOS_BIAS\tEXCHET"); \
(bcftools query -f "%CHROM\t%POS\t%INFO/DP\t%INFO/MQ\t%INFO/PV4{0}\t%INFO/PV4{1}\t%INFO/PV4{2}\t%INFO/PV4{3}\t%ExcHet\n" $DIR/output/calmas_allsites.bcf.gz)) > $DIR/output/allsites_stats.txt
```

Filter the VCF. We'll avoid dumping another VCF with just sites that pass our quality controls by extracting just the sites.
We can use these sites with ANGSD tomorrow.

| bcftools view -T ^/space/s2/diana/Rimitator/sites_290421/exclude_sites.bed -i 'N_PASS(FMT/DP[0-32] > 2) > 14 && N_PASS(FMT/DP[33-65] > 2) > 14 && N_PASS(FMT/DP[66-123] > 2) > 14' 
| bcftools view -e 'INFO/BandMisMap > 500 || INFO/StripeMisMap > 500 || INFO/AdmixMisMap > 500 || INDEL=1 || INFO/MQ < 30 || INFO/DP > 20000 || INFO/PV4[0] < 1e-100 || INFO/PV4[1] < 1e-100 || INFO/PV4[2] < 1e-100 || INFO/PV4[3] < 1e-100 || INFO/ExcHet < 1e-3' -M 2 | bcftools query -f "%CHROM\t%POS\n" > /space/s2/diana/Rimitator/sites_290421/imi_model_allsites_set1.pos

bcftools view -i 'N_PASS(FMT/DP[0-14] > 2) > 5 && N_PASS(FMT/DP[15-39] > 2) > 5' $DIR/output/calmas_allsites.bcf.gz | bcftools view -e 'INDEL=1 || INFO/MQ < 30 || INFO/DP > 4000 || INFO/PV4[0] < 1e-100 || INFO/PV4[1] < 1e-100 || INFO/PV4[2] < 1e-100 || INFO/PV4[2] < 1e-100 || INFO/PV4[3] < 1e-100 || INFO/ExcHet < 1e-3' -M 2 | bcftools query -f "%CHROM\t%POS\n" > $DIR/output/qc_sites.pos


You could bypass dumping the initial VCF entire using one large strong pipes. Give it a try.
