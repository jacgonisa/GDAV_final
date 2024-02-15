##TASK4

ln 02rna_seq/results/*.sorted.bam 04de_analysis

ln 02rna_seq/ids_rnaseq.txt 04de_analysis

mkdir 04de_analysis/results
mkdir 04de_analysis/log

## This command create the .count and .log files successfully
cat 04de_analysis/ids_rnaseq.txt | parallel '(htseq-count -i locus_tag -t CDS 04de_analysis/{}.sorted.bam genome.gff > 04de_analysis/results/{}.count) 2> 04de_analysis/log/{}.count.err'

## Now let's start R scripting
mkdir 04de_analysis/scripts

vi 04de_analysis/scripts/DESeq2.R #and I customize it

Rscript 04de_analysis/scripts/DESeq2.R 1> 04de_analysis/log/Rscript_DESeq2.out

cat 04de_analysis/results/deseq2_results_padj_0.05.csv | column -t

embl_gene_id","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","sig"
"1","AQUIFEX_01423",10127.8125528592,5.05681872544562,0.055800558144404,90.6230850300686,0,0,"yes"
"2","AQUIFEX_01759",3004.72718330008,4.97020342086057,0.0887472820291642,56.0040071900711,0,0,"yes"
"3","AQUIFEX_01761",2793.70139675557,4.95252945473872,0.09156565942088,54.0871925792015,0,0,"yes"
"4","AQUIFEX_01723",153.830592373295,-10.7143610248187,1.45250300474677,-7.37648114310558,1.62527893113547e-13,6.12730157038071e-11,"yes"
"5","AQUIFEX_01749",71.9091835372515,-9.61731838535804,1.46325507821799,-6.57255083445218,4.94604928247163e-11,1.49172846359344e-08,"yes"

#4.1. The description of the fields is in the document

#4.2. Visually, it was clear that only 5 genes whose p-adj < 0.05
##The total of genes?  

wc -l 04de_analysis/results/deseq2_results_total.csv
1509 -1(header)= 1508


##I generate a file with gene IDs
cat 04de_analysis/results/deseq2_results_padj_0.05.csv | awk -F',' '{print $2"\t"$4}' | tail -n +2 | sed 's#\"##g' | sed 's#\(.*\)#ID=\1#g' > 04de_analysis/results/degs_ids_logfc.txt

##The intersection
grep -Ff <(cut -f1 04de_analysis/results/degs_ids_logfc.txt) genome.gff > 04de_analysis/results/degs_togff.gff

##Check the functions(the protein products)
cut -f9 04de_analysis/results/degs_togff.gff | grep -o -P 'product=\K\.*'
Nif-specific regulatory protein
hypothetical protein
hypothetical protein
FeMo cofactor biosynthesis protein NifB
Nitrogenase iron protein 1

#4.2
##first I get the clean IDs
cat 04de_analysis/results/deseq2_results_padj_0.05.csv | awk -F',' '{print $2}' | sed 's#\"##g' > 04de_analysis/results/degs_ids.txt
