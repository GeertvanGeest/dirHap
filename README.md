# dirHap


## Counting haplotypes
The script haploCounter.py counts the haplotypes in a specific region.
All the reads in the bam file should overlap the entire region specified. 

Filter completely overlapping reads e.g. with bedtools:
```
bedtools intersect -F 1 -abam alignment.bam -b regions.bed > alignment.intsc.bam
```

Haplocounter works with a specific sample, on a specific region. 
A vcf file denotes the SNPs used for haplotyping. 

Example usage:
```
haploCounter.py \
-i example.bam \
-s s501 \
-vcf example.vcf.gz \
-r ST4.03ch01:29682248-29682338 \
-c 5 \
-wmode w \
-o haplotypes.txt
```
