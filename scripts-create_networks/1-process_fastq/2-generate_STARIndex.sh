#!/bin/bash
STAR --runThreadN 16 --runMode genomeGenerate --genomeDir S_cerevisiae_Ensembl_R64-1-1/STARIndex/ --genomeFastaFiles S_cerevisiae_Ensembl_R64-1-1/Sequence/WholeGenomeFasta/genome.fa --sjdbGTFfile S_cerevisiae_Ensembl_R64-1-1/Annotation/Genes/genes.gtf --sjdbOverhang 49 --genomeSAindexNbases 10;
