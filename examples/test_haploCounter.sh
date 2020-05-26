#!/bin/bash
../haploCounter.py \
-i example.bam \
-s s501 \
-vcf example.vcf.gz \
-r 'ST4.03ch01:29682248-29682338' \
-c 5 \
-wmode w \
-o haplotypes.txt
