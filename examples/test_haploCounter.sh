#!/bin/bash
../haploCounter.py \
-i example.bam \
-s s501 \
-vcf example.vcf.gz \
-bed region.bed \
-c 5 \
-o haplotypes.txt
