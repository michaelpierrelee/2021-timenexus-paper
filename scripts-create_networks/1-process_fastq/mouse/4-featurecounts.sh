#!/bin/bash

featureCounts -O -p -T 16 -F GTF -t exon -g gene_id -a M_musculus_UCSC_mm10/Annotation/Genes/genes.gtf -o counts/DRG_WT_timecourse.txt bam_trimmed/*Aligned.out.bam
