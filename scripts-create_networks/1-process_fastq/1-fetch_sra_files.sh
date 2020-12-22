#!/bin/bash
filename='AccessionList_PRJNA319029_yeast_WT.txt'
while read p; do 
	echo "---" $p "---";
	#fasterq-dump $p -O fastq/ -e 8 -m 20000 -p;
done < $filename
