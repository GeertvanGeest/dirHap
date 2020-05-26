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

and numpy and pandas using pip:
```
pip install numpy
pip install pandas
pip install scipy

```

or conda:
```
conda install -c anaconda numpy 
conda install -c anaconda scipy 
conda install -c anaconda pandas 

```

## Counting haplotypes
The script haploCounter.py counts the haplotypes in a specified region.
All the reads in the bam file should overlap the entire region specified. 

Filter completely overlapping reads e.g. with bedtools:
```
bedtools intersect -F 1 -abam alignment.bam -b regions.bed > alignment.intsc.bam
```

Haplocounter works with a specified sample, on a specified region. 
A vcf file denotes the SNPs used for haplotyping. 

Here's an example:
```
python haploCounter.py \
-i ./examples/example.bam \
-s s501 \
-vcf ./examples/example.vcf.gz \
-bed ./examples/region.bed \
-c 5 \
-o ./examples/haplotypes.txt
```

This gives the following output (haplotypes.txt), which is a tab-delimited txt file, where columns correspond to "sample  chrom start end seq count"

```
s501	ST4.03ch01	29682248	29682338	ACA	408
s501	ST4.03ch01	29682248	29682338	CCA	894
s501	ST4.03ch01	29682248	29682338	CTA	442
s501	ST4.03ch01	46274385	46274478	CCCC	538
s501	ST4.03ch01	46274385	46274478	CCCG	198
s501	ST4.03ch01	46274385	46274478	CACG	11
```

## Calculating dosages

The dosages are calculated based on the `haplotypes.txt` file that is extracted using the script `haploCounter.py`. 
To run the example with the default options execute the following command
```
python dosageEstimator.py \
--haplotype_data ./examples/haplotypes.txt \
--output_file ./examples/dosages.txt
```

For each sample and region we generate (haplotype) dosages. 
This gives the following output (haplotypes.txt)
```
sample	chrom	start	end	seq	count	prob_true	freq	used	BA_dosage	BA_LR_ratio	MA_dosage	MA_LR_ratio
s501	ST4.03ch01	29682248	29682338	CC	895	0.9999999999999999	0.5128939828080229	True	2.0	226.2838317191463	2	313.9956727936551
s501	ST4.03ch01	29682248	29682338	AC	408	0.9999999999999999	0.233810888252149	True	1.0	259.3027998721581	1	313.9956727936551
s501	ST4.03ch01	29682248	29682338	CT	442	0.9999999999999999	0.2532951289398281	True	1.0	221.94998205744247	1	313.9956727936551
```

where several columns are added to the original `haplotypes.txt` file:

1. **prob_true:** Probability that a haplotype is potentiallly a sequencing error
2. **freq:** read count of haplotype / total depth
3. **BA_dosage:** Dosage estimation under bi-allelic model
4. **BA_LR_ratio:** LR ratio for dosage estimation under bi-allelic model
5. **MA_dosage:** Dosage estimation under multi-allelic model
6. **MA_LR_ratio:**  LR ratio for dosage estimation under multi-allelic model

These results can be filtered in any way you like. 