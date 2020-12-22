#!/bin/bash

for FILE in $(ls fastq/*.fastq); do
	fastqc $FILE -o fastqc/ -t 16
done

multiqc fastqc/* -n fastqc/multiqc_raw;
