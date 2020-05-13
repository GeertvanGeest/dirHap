# dirHap

## Installation
Make a conda environment with e.g. python3.7:
```
conda create -n pysam_env_py37 python=3.7.* 
conda activate pysam_env_py37
```

install pysam:
```
conda config --add channels bioconda
conda install pysam
``` 

and numpy using pip:
```
pip install numpy
```

or conda:
```
conda install -c anaconda numpy 
```

## Counting haplotypes
The script haploCounter.py counts the haplotypes in a specific region.
All the reads in the bam file should overlap the entire region specified. 

Filter completely overlapping reads e.g. with bedtools:
```
bedtools intersect -F 1 -abam alignment.bam -b regions.bed > alignment.intsc.bam
```

Haplocounter works with a specific sample, on a specific region. 
A vcf file denotes the SNPs used for haplotyping. 

Here's the usage:

```
usage: haploCounter.py [-h] -i INPUT.BAM -s SAMPLEID -vcf SNPS.VCF -r
                       CHROM:START-END [-c C] [-o O] [-wmode WMODE]

Processing of BAM files into haplotype counts.

optional arguments:
  -h, --help          show this help message and exit
  -i INPUT.BAM        Input BAM file. Reads should cover entire specified
                      region in -r
  -s SAMPLEID         Sample ID
  -vcf SNPS.VCF       Input (tabix-indexed) VCF/BCF file
  -r CHROM:START-END  Region chr:start-end
  -c C                Minimum number of reads supporting haplotype
  -o O                Output count file
  -wmode WMODE        Output write mode. E.g. 'a' for append or 'w' for write

Copyright GvG
```

And an example:
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

This gives the following output (haplotypes.txt)
```
s501    ST4.03ch01      29682248        29682338        CC      895
s501    ST4.03ch01      29682248        29682338        AC      408
s501    ST4.03ch01      29682248        29682338        CT      442
```
