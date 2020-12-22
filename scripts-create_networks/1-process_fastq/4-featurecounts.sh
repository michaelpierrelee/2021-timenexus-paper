#!/bin/bash

featureCounts -O -T 16 -F GTF -t exon -g gene_id -a S_cerevisiae_Ensembl_R64-1-1/Annotation/Genes/genes.gtf -o counts/yeast_WT_counts.txt bam/*Aligned.out.bam;
