# RNAseq
# TASK1
```
mkdir 02rna_seq
```

## Question 2.1 
The directory 'RNAseq' is quite simple, no need of "fancy" stuff as one can easily see there are 8 files, corresponding to 4 samples (2 conditions and 2 reps per condition). Anyway:
```
ls RNAseq | wc -l ##to count number of files
```
## Question 2.2
I could do it using awk as we have seen in class: 
```
zcat RNAseq/hightemp01.r1.fq.gz | awk 'END{print NR/4}'
```

or using another similar strategy:
```
echo $(zcat RNAseq/hightemp01.r1.fq.gz | wc -l) /4|bc      
```
from https://www.biostars.org/p/139006/

However, I will use a high-efficient specialized package called "seqkit"

```
conda create -n seqkit
conda activate seqkit
conda install -c bioconda seqkit
```
```
seqkit stats RNAseq/* > 02rna_seq/summary_samples.txt

cat 02rna_seq/summary_samples.txt 
```

file                        format  type  num_seqs     sum_len  min_len  avg_len  max_len
RNAseq/hightemp01.r1.fq.gz  FASTQ   DNA    318,719  31,871,900      100      100      100
RNAseq/hightemp01.r2.fq.gz  FASTQ   DNA    318,719  31,871,900      100      100      100
RNAseq/hightemp02.r1.fq.gz  FASTQ   DNA    318,719  31,871,900      100      100      100
RNAseq/hightemp02.r2.fq.gz  FASTQ   DNA    318,719  31,871,900      100      100      100
RNAseq/normal01.r1.fq.gz    FASTQ   DNA    288,742  28,874,200      100      100      100
RNAseq/normal01.r2.fq.gz    FASTQ   DNA    288,742  28,874,200      100      100      100
RNAseq/normal02.r1.fq.gz    FASTQ   DNA    288,742  28,874,200      100      100      100
RNAseq/normal02.r2.fq.gz    FASTQ   DNA    288,742  28,874,200      100      100      100

## Question 2.3. Paired-end reads

```
conda create env -n rnseqc python=3.10
mamba install -c bioconda rnseqc
mamba install -c bioconda bedops
gff2bed < genome.gff > genome.bed
infer_experiment.py -r genome.bed -i 02rna_seq/results/hightemp01.bam > 02rna_seq/log/infer_experiment.out
```


## Question 2.4. 
```
cat 02rna_seq/summary_samples.txt
```

## Question 2.5
Answer found in final report (.pdf file).

## Question 2.6
```
mamba install -c bioconda fastqc

mkdir 02rna_seq/fastqc

find RNAseq/ | parallel -j 8 "fastqc {} --outdir 02rna_seq/fastqc/."
```

in another environment with multiqc
```
multiqc 02rna_seq/fastqc/. --outdir 02rna_seq/fastqc/
```
with rsync I transfer to my laptop the multiqc report. Relevant image is attached in my report. 
Anyways, it is not absolutely necessary to do this, because with one quick look:

```
zcat RNAseq/hightemp01.r1.fq.gz | head -n 10
```
one can notice that all the PHRED scores are the same (9, which in ASCII is "I")

#######

# TASK2

index the genome
```
bwa index genome.fasta
```
# TASK 3 and 4 - I am doing them in one line
Map sample to the reference genome and create BAM files
I create first an ids.txt to use parallel
```
ls RNAseq/* | cut -d '.' -f1 | sed 's#.*\/##'| uniq  > 02rna_seq/ids_rnaseq.txt   #(it would have been easier to echo the four names...)

cat 02rna_seq/ids_rnaseq.txt | parallel --compress '(bwa mem genome.fasta RNAseq/{}.r1.fq.gz RNAseq/{}.r2.fq.gz | samtools sort -o 02rna_seq/results/{}.bam) 2> 02rna_seq/log/{}_run_bwamem.err'
```

## Question 2.7
sort and index BAM files
```
cat 02rna_seq/ids_rnaseq.txt | parallel -j 1 'samtools sort 02rna_seq/results/{}.bam > 02rna_seq/results/{}.sorted.bam; samtools index 02rna_seq/results/{}.sorted.bam' 
```

Using samtools
```
cat 02rna_seq/ids_rnaseq.txt | parallel -j 1 '(samtools flagstat 02rna_seq/results/{}.sorted.bam) 1> 02rna_seq/results/{}.samtools.flagstat

cat 02rna_seq/results/hightemp01.samtools.flagstat
```

638004 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
566 + 0 supplementary
0 + 0 duplicates
638004 + 0 mapped (100.00% : N/A)
637438 + 0 paired in sequencing
318719 + 0 read1
318719 + 0 read2
634444 + 0 properly paired (99.53% : N/A)
637438 + 0 with itself and mate mapped
0 + 0 singletons (0.00% : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)

Using bamtools

```
cat 02rna_seq/ids_rnaseq.txt | parallel -j 1 '(bamtools stats -in 02rna_seq/results/{}.sorted.bam) 1> 02rna_seq/results/{}.bamtools.stats'

cat 02rna_seq/results/hightemp01.bamtools.stats
```
Total reads:       638004
Mapped reads:      638004	(100%)
Forward strand:    318998	(49.9994%)
Reverse strand:    319006	(50.0006%)
Failed QC:         0	(0%)
Duplicates:        0	(0%)
Paired-end reads:  638004	(100%)
'Proper-pairs':    634651	(99.4745%)
Both pairs mapped: 638004	(100%)
Read 1:            319024
Read 2:            318980
Singletons:        0	(0%)

Another way of just getting the mapped reads (I include the supplementary mappings)
```
cat 02rna_seq/ids_rnaseq.txt | parallel -j 1 'printf "{}\t" &&  samtools view -c -F 4  02rna_seq/results/{}.sorted.bam'
```
normal01	577508
normal02	577503
hightemp01	638004
hightemp02	637997

```
cat 02rna_seq/ids_rnaseq.txt | parallel 'printf "{}\t" &&  samtools view -c -F 4 -F 2048 02rna_seq/results/{}.sorted.bam'
```
normal01	577484
normal02	577484
hightemp01	637438
hightemp02	637438



## Question 2.8
```
samtools coverage 02rna_seq/results/hightemp01.sorted.bam | column -t
```
we can check many the coverage and depth stats and they are OK but...
Copy number variations with RNA-seq? It seems difficult
