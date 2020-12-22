#!/bin/bash
for j in 00 01 03 30; do
	for m in 1 2; do
	
		NAME=DRG_WT_${j}_Rep${m};
		
		# CUTADAPT

		cutadapt -j 16 -m 1 -a AAGCAGTGGTATCAACGCAGAGTACTTTTTTTTTTTTTTTTTTTTTTTTT -a AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA -a GATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGACGGTGATCTCGTA -A AAGCAGTGGTATCAACGCAGAGTACTTTTTTTTTTTTTTTTTTTTTTTTT -A AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA -A GATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGACGGTGATCTCGTA -o fastq_trimmed/${NAME}_R1_trimmed.fastq -p fastq_trimmed/${NAME}_R2_trimmed.fastq fastq/${NAME}_R1.fastq fastq/${NAME}_R2.fastq;

		# FASTQC

		fastqc fastq_trimmed/${NAME}_R1_trimmed.fastq -o fastqc_trimmed/ -t 16;
		fastqc fastq_trimmed/${NAME}_R2_trimmed.fastq -o fastqc_trimmed/ -t 16;

		# STAR

		STAR --runThreadN 16 --genomeDir M_musculus_UCSC_mm10/STARIndex/ --readFilesIn fastq_trimmed/${NAME}_R1_trimmed.fastq fastq_trimmed/${NAME}_R2_trimmed.fastq --outFileNamePrefix bam_trimmed/${NAME} --sjdbGTFfile M_musculus_UCSC_mm10/Annotation/Genes/genes.gtf --outSAMstrandField intronMotif --outSAMtype BAM Unsorted --outReadsUnmapped Fastx;
		
	done
done
