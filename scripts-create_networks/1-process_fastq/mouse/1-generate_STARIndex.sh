#!/bin/bash
STAR --runThreadN 16 --runMode genomeGenerate --genomeDir M_musculus_UCSC_mm10/STARIndex/ --genomeFastaFiles M_musculus_UCSC_mm10/Sequence/WholeGenomeFasta/genome.fa --sjdbGTFfile M_musculus_UCSC_mm10/Annotation/Genes/genes.gtf;
