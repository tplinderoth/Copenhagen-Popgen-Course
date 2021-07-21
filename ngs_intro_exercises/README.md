INTRODUCTION TO NGS DATA EXERCISES
==================================

Set some environmental variables (assume bams are in ${DIR}/data/bams and fastq are in ${DIR}/data/fastq)

	DIR=/home/tyler/Desktop/workshop_tmp
	mkdir "$DIR/output"
	TRIMJAR=/home/tyler/install/Trimmomatic-0.39/trimmomatic-0.39.jar

## fastq

Have a look inside one of the frog fastq files

	less "$DIR/data/fastq/CH1401_R1.fastq"

visualize data with fastqc. You can get some help info with `fastqc --help`
Will automatically detect file format type, but can specify with `-f`.

	fastqc -outdir "$DIR/output/" $DIR/data/fastq/*.fastq
	ll "$DIR/output/"

normally you can scp the html files onto your local machine and view them in the web. It would look something like:

![fastq_preclean_example](./outputs/CH1401_R2_before.png)

Perhaps we have some adapter contamination to worry about so we can check this out without the browser interface

	unzip "$DIR/output/CH1401_R2_fastqc.zip" -d "$DIR/output/"
	
	# check out summary
	cat "$DIR/output/CH1401_R2_fastqc/summary.txt"

	# visualize the extent of adapter contamination (this can take a moment to show up, be patient...)
	eog "$DIR/output/CH1401_R2_fastqc/Images/adapter_content.png"

Can you see what fastq is picking up on? About what percentage of reads show contamination?

## clean fastq

We'll see if we can clean up the adapter contamination and low quality at the ends of the read using cutadapt. Normally run in paired-end mode so that
read pairs are retained. Single-end mode can be run on the forward and reverse reads separately if it's desirable to retain unpaired reads. We'll run on
just the reverse read that we examined for sake of time. Partial sums low-quality trimming using -q 15 (algorithm: https://cutadapt.readthedocs.io/en/stable/algorithms.html#quality-trimming-algorithm)
Specifing Illumina TruSeq adapters (https://support.illumina.com/bulletins/2016/12/what-sequences-do-i-use-for-adapter-trimming.html).

	cutadapt3 -q 15 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 36 \
	-o "$DIR/output/CH1401_R2_clean.fastq.gz" "$DIR/data/fastq/CH1401_R2.fastq"

use fastqc to check what the data quality looks like for CH1401_R2.fasq after trimming

	fastqc -outdir "$DIR/output/" "$DIR/output/CH1401_R2_clean.fastq.gz"

We can see that cutadapt did a good job of cleaning up the data

![fastq_postclean_example](./outputs/CH1401_R2_after.png)

## Working with mapped data

After you've mapped your cleaned-up fastq reads, you can have a look at the mapping information in SAM/BAM/CRAM files with samtools.
You can include the header with `-h`. Note that if they aren't already, you will want to sort your bam files with `samtools sort`, but ours are
already sorted so we'll have a look.

	samtools view "$DIR/data/bams/CMASS6169443.bam" | less -S

the first 5 lines, `samtools view "$DIR/data/bams/CMASS6169443.bam | head -n 5` should look like

	HS30_18456:2:2102:20331:79677#7	99	chr7	73	57	125M	=	407	459	ACGAGCAAGAGGACAGTTTCAAAGCAGACATTATGGCACAATAAAAGCTGTGAGAGCAGAGATTCATTCATATCTTTACTATTTTTCTTTAGCTTCTACTCCCGCTGCTTACAATAAATCTGAGA	EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE	AS:i:125	XS:i:115	BC:Z:CAGATCTGTCTTTCCC	QT:Z:BBBBBFFF/BBBBFFF	MQ:i:40	MC:Z:125M	ms:i:4579	MD:Z:125	NM:i:0	RG:Z:18456_2#7
	HS30_18456:2:1213:20056:15084#7	163	chr7	75	0	125M	=	569	619	GAGCAAGAGGACAGTTTCAAAGCAGACATTATGGCACAATAAAAGCTGTAAGAGCAGAGATTCATTCATATCTTTACTATTTTTCTTTAGCTTCTACTCCCGCTGCTTACAATAAATCCGAGATA	DDDDDDDDDDDDDDD/AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA/AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA/DDDDDD/DDDD/FF/B5555F/B555FF	AS:i:115	XS:i:125	XA:Z:chr8,-23982502,125M,0;chr16,-2795516,125M,1;chr18,+31933308,125M,1;chr18,+30035828,125M,2;chr8,-23993408,125M,3;	MQ:i:9	MC:Z:125M	ms:i:4159	MD:Z:49G68T6	NM:i:2	RG:Z:18456_2#7
	HS30_18456:2:2211:1719:11546#7	163	chr7	103	51	125M	=	391	413	TTATGGCACAATAAAAGCTGTGAGAGCAGAGATTCATTCATATCTTTACTATTTTTCTTTAGCTTCTACTCCCGCTGCTTACAATAAATCTGAGATAAAGAACAAACAAGAACAACTCACATGGC	EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE	AS:i:125	XS:i:115	MQ:i:40	MC:Z:125M	ms:i:4597	MD:Z:125	NM:i:0	RG:Z:18456_2#7
	HS30_18456:2:2101:10435:101060#7	99	chr7	118	0	125M	=	499	506	AGCTGTAAGAGCAGAGATTCATTCATATCTTTACTATTTTTCTTTAGCTTCTACTCCCGCTGCTTACAATAAATCCGAGATAAAGAACAAACAAGAACAACTCACATGGCATTGATTGTTTAGTT	>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>/F	AS:i:115	XS:i:125	XA:Z:chr8,-23982459,125M,0;chr18,+30035871,125M,0;chr16,-2795473,125M,1;chr18,+31933351,125M,1;U_scaffold_23,-9073,125M,2;	BC:Z:CAGATCTGTCTTTCCC	QT:Z:BBBBBFFF/7<BB/B/	MQ:i:8	MC:Z:125M	ms:i:4556	MD:Z:6G68T49	NM:i:2	RG:Z:18456_2#7
	HS30_18456:2:1316:11207:16667#7	163	chr7	126	0	125M	=	459	458	GAGCAGAGATTCATTCATATCTTTACTATTTTTCTTTAGCTTCTACTCCCGCTGCTTACAATAAATCCGAGATAAAGAACAAACAAGAACAACTCACATGGCATTGATTGTTTAGTTCAGTGTCA	DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD/	AS:i:120	XS:i:125	MQ:i:6	MC:Z:125M	ms:i:4546	MD:Z:67T57	NM:i:1	RG:Z:18456_2#7

A useful way to visualize mapping information is with the Integrative Genomics Viewer (IGV)

## software
fastqc: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
cutadapt: https://cutadapt.readthedocs.io/en/stable/index.html
