#!/bin/bash

for FILE in fastq/*
	do
		NAME=${FILE##*/};
		echo "\n---" $NAME "---";
		echo ">> FASTQC"
		# FASTQC
		fastqc ${FILE} -o fastqc/;

		# STAR
		echo "\n>> STAR"
		STAR --runThreadN 16 --genomeDir S_cerevisiae_Ensembl_R64-1-1/STARIndex/ --readFilesIn ${FILE} --outFileNamePrefix bam/${NAME} --sjdbGTFfile S_cerevisiae_Ensembl_R64-1-1/Annotation/Genes/genes.gtf --outSAMstrandField intronMotif --outSAMtype BAM Unsorted --outReadsUnmapped Fastx;
	done
	
#FASTQ_FILES=$(ls fastq/ -xm | tr -d ' ')
#STAR --runThreadN 16 --genomeDir S_cerevisiae_Ensembl_R64-1-1/STARIndex/ --readFilesIn ${FASTQ_FILES} --readFilesPrefix fastq/ --outFileNamePrefix bam/ --sjdbGTFfile S_cerevisiae_Ensembl_R64-1-1/Annotation/Genes/genes.gtf --outSAMstrandField intronMotif --outSAMtype BAM Unsorted --outReadsUnmapped Fastx;
