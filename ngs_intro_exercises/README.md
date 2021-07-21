INTRODUCTION TO NGS DATA EXERCISES
==================================

Set some environmental variables (assume bams are in ${DIR}/data/bams and fastq are in ${DIR}/data/fastq)

	DIR=/home/tyler/Desktop/workshop_tmp
	mkdir "$DIR/output"

Have a look inside one of the frog fastq files

	less "$DIR/data/fastq/CH1401_R1.fastq"

visualize data with fastqc. You can get some help info with `fastqc --help`
Will automatically detect file format type, but can specify with `-f`.

	fastqc -outdir "$DIR/output/" $DIR/data/fastq/*.fastq

normally you can scp the html files onto your local machine and view them in the web. It would look something like:

![fastq_preclean_example](CH1401_R2_before.png)

## software
fastqc: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
