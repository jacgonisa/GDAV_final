# TASK3

## 1.Sorting already done

## 2.Merging mappings

```
samtools merge 03variant_calling/all_merged.bam 02rna_seq/results/*.sorted.bam

samtools view -c 03variant_calling/all_merged.bam
```

2431012

Actually, if you sum the total reads from the individual .bam files, it makes 2431017, so 5 reads are missing

```
samtools sort -o 03variant_calling/all_merged.sorted.bam 03variant_calling/all_merged.bam

samtools index 03variant_calling/all_merged.sorted.bam
```

## 3. Variant calling
```
mkdir 03variant_calling/results
mkdir 03variant_calling/logs

bcftools mpileup --threads 1 -f genome.fasta 03variant_calling/all_merged.sorted.bam | bcftools call --threads 1 --ploidy 1 -mv -Ob -o  03variant_calling/results/calls.bcf 2> 03variant_calling/logs/bcftools_call.err


bcftools mpileup --threads 1 -f genome.fasta 03variant_calling/all_merged.sorted.bam | bcftools call --threads 1 -mv -Ob -o  03variant_calling/results/calls_diploid.bcf 2> 03variant_calling/logs/bcftools_call_diploid.log
```
Here I am producing the log, but nothing important, just the message related to missing --ploidy

```
bcftools mpileup --threads 1 -f genome.fasta 02rna_seq/results/*.sorted.bam | bcftools call --threads 1 --ploidy 1 -mv -Ob -o  03variant_calling/results/calls_4files.bcf 2> 03variant_calling/logs/bcftools_call_4files.log

bcftools mpileup --threads 1 -f genome.fasta 02rna_seq/results/*.sorted.bam | bcftools call --threads 1  -mv -Ob -o  03variant_calling/results/calls_4files_diploid.bcf 2> 03variant_calling/logs/bcftools_call_4files_diploid.err
```

## 3.1. Variants detected

First I make a file with IDs

```
find 03variant_calling/results/ -name '*.bcf' -exec basename {} .bcf \; > 03variant_calling/results/ids_vcalling.txt
```

I make now a summary table
```
(printf "%s\t%s\t%s\t%s\n" "file_name" "variants" "SNPs" "indels"; cat 03variant_calling/results/ids_vcalling.txt | parallel -j1 --colsep '\t' '
    sample={1};
    variants=$(bcftools view -H "03variant_calling/results/{}.bcf" | wc -l);
    snps=$(bcftools view -H -v snps "03variant_calling/results/{}.bcf" | wc -l);
    indels=$(bcftools view -H -v indels "03variant_calling/results/{}.bcf" | wc -l);
```

```
    printf "%s\t%s\t%s\t%s\n" "$sample" "$variants" "$snps" "$indels";
')> 03variant_calling/results/calling_summary.tsv
```
```
cat 03variant_calling/results/calling_summary.tsv | column -t
```

file_name             variants  SNPs  indels
calls                 2         2     0
calls_diploid         3         3     0
calls_4files          2         2     0
calls_4files_diploid  3         3     0


# 3.2 Quality greater than 100
```
cat 03variant_calling/results/ids_vcalling.txt | parallel -j 1 'echo {} ; bcftools filter -e "QUAL<100" 03variant_calling/results/{}.bcf ; echo " "' | grep -v "##"
```

I already see 2
However, we can do another thing to get the numerical value more easily

```
cat 03variant_calling/results/ids_vcalling.txt | parallel -j 1 'bcftools filter -e "QUAL<100" 03variant_calling/results/{}.bcf > 03variant_calling/results/{}_filter_Q100.bcf'

cat 03variant_calling/results/ids_vcalling.txt | parallel -j 1 bcftools view -H 03variant_calling/results/filter_Q100_{}.bcf | wc -l
```

The output is 2

## 3.3 Depth of coverage greater than 100
```
cat 03variant_calling/results/ids_vcalling.txt | parallel -j 1 'bcftools filter -e "DP<100" 03variant_calling/results/{}.bcf > 03variant_calling/results/filter_DP100_{}.bcf'
 
cat 03variant_calling/results/ids_vcalling.txt | parallel -j 1 bcftools view -H 03variant_calling/results/filter_DP100_{}.bcf | wc -l
```

Again the output is 2, and the variants are the same actually


## 3.4 Differences en ploidy

```
cat 03variant_calling/results/ids_vcalling.txt | parallel -j 1 'echo {}; bcftools view -H 03variant_calling/results/filter_Q100_{}.bcf | wc -l; echo " " '
```

calls
0
 
calls_diploid
1
 
calls_4files
0
 
calls_4files_diploid
1
 

## 3.5 The best quality variant

```
parallel -j 1 'echo {}; bcftools view -H -m2 -M2 03variant_calling/results/filter_Q100_{}.bcf ; echo " " ' ::: calls_diploid calls_4files_diploid
```

calls_diploid
Aquifex_genome	1265060	.	C	T	188	PASS	DP=249;VDB=3.5208e-10;SGB=-0.693147;RPB=0.944775;MQB=1;BQB=0.984769;MQ0F=0;ICB=1;HOB=0.5;AC=1;AN=2;DP4=136,0,112,0;MQ=60	GT:PL	0/1:221,0,249
 
calls_4files_diploid
Aquifex_genome	1265060	.	C	T	950	PASS	DP=674;VDB=1.01192e-06;SGB=-101.196;RPB=0.969374;MQB=1;MQSB=1;BQB=0.985692;MQ0F=0;ICB=0.166667;HOB=0.5;AC=4;AN=8;DP4=294,55,272,45;MQ=60	GT:PL	0/1:223,0,249	0/1:255,0,255	0/1:255,0,255	0/1:255,0,255

calls_4files_diploid has clearly the best quality variant (950 quality and DP=674)
It seems heterozygous 

I am converting into vcf to visualize in IGV

```
parallel -j 1 'bcftools convert -O v -o 03variant_calling/results/{}.vcf 03variant_calling/results/{}.bcf' ::: calls_diploid calls_4files_diploid
```

Yes, it could be affecting to a gene whose product is "Nif-specific regulatory protein". I inspected this visually but annotation can be done:
```
mkdir 03variant_calling/results/intersection_gff
cat 03variant_calling/results/ids_vcalling.txt | parallel -j 1 "bedtools intersect -loj -a <(grep 'CDS' genome.gff) -b <(bcftools view 03variant_calling/results/filter_Q100_{}.bcf) | awk -F $'\t' '\$10!=\".\" {print \$0}' > 03variant_calling/results/intersection_gff/filter_Q100_{}_intersection_with_bcf.gff"

cat 03variant_calling/results/intersection_gff/*
```
Aquifex_genome	Prodigal:002006	CDS	1264237	1265730	.	-	0	ID=AQUIFEX_01423;Name=nifA;gene=nifA;inference=ab initio prediction:Prodigal:002006,similar to AA sequence:UniProtKB:P03027;locus_tag=AQUIFEX_01423;product=Nif-specific regulatory protein	Aquifex_genome	1265060	.	C	T	950	PASS	DP=674;VDB=1.01192e-06;SGB=-101.196;RPB=0.969374;MQB=1;MQSB=1;BQB=0.985692;MQ0F=0;ICB=0.166667;HOB=0.5;AC=4;AN=8;DP4=294,55,272,45;MQ=60	GT:PL	0/1:223,0,249	0/1:255,0,255	0/1:255,0,255	0/1:255,0,255
Aquifex_genome	Prodigal:002006	CDS	1264237	1265730	.	-	0	ID=AQUIFEX_01423;Name=nifA;gene=nifA;inference=ab initio prediction:Prodigal:002006,similar to AA sequence:UniProtKB:P03027;locus_tag=AQUIFEX_01423;product=Nif-specific regulatory protein	Aquifex_genome	1265060	.	C	T	188	PASS	DP=249;VDB=3.5208e-10;SGB=-0.693147;RPB=0.944775;MQB=1;BQB=0.984769;MQ0F=0;ICB=1;HOB=0.5;AC=1;AN=2;DP4=136,0,112,0;MQ=60	GT:PL	0/1:221,0,249




