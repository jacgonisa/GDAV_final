#!/usr/bin/env Rscript

library(DESeq2);
library(tibble);
cat(" ")


cat("###SET UP###\n")

# Create an object with the directory containing HTseq counts:
directory <- "/home/2023/gdav/jacob.gonzalez.isa.gdav2023/final_project/04de_analysis/results"

cat("Looking for .count files in " ,list.files(directory))

sampleFiles <- list.files(directory, pattern = "\\.count$")
sampleFiles

# create a vector of sample names. Ensure these are in the same order as the "sampleFiles" object!
sampleNames <- c( "hightemp01", "hightemp02","normal01","normal02")

# create a vector of conditions. again, mind that they are ordered correctly!
replicate <- c("Rep1","Rep2","Rep1","Rep2")

# create a vector of conditions. again, mind that they are ordered correctly!
sampleCondition <- c("hightemp","hightemp","normal", "normal")

# now create a data frame from these three vectors. 
sampleTable <- data.frame(
		sampleName = sampleNames,
		fileName = sampleFiles,
		condition = sampleCondition,
		replicate = replicate)

sampleTable

## Make DESeq2 object from counts and metadata

ddsHTSeq <- DESeqDataSetFromHTSeqCount(
		sampleTable = sampleTable, 
		directory = directory, 
		design = ~condition)

# specify the reference level:

ddsHTSeq$condition <- relevel(ddsHTSeq$condition, ref = "normal")


cat("###MINIMAL PREFILTERING###\n") 

# sum counts for each gene across samples

sumcounts <- rowSums(counts(ddsHTSeq))

# get genes with summed counts greater than 10; remove lowly expressed genes

keep <- sumcounts > 10

# keep only the genes for which the vector "keep" is TRUE

ddsHTSeq_filter <- ddsHTSeq[keep,]

cat("###RUNNING DESEQ2###\n")

dds <- DESeq(ddsHTSeq_filter)

# get results table

res <- results(dds, pAdjustMethod="BH") 

summary(res)


# check out the first few lines

head(res)

mcols(res, use.names = T)

resultsNames(dds)

cat("###CREATING OUTPUT###\n")

##Create normalized read counts

normalized_counts <- counts(dds, normalized=TRUE)

normalized_counts_mad <- apply(normalized_counts, 1, mad)
normalized_counts <- normalized_counts[order(normalized_counts_mad, decreasing=T), ]

#DESeq get results table

Res_A_X_total <- results(dds, name="condition_hightemp_vs_normal", pAdjustMethod="BH")
Res_A_X_total <- Res_A_X_total[order(Res_A_X_total$padj),]
Res_A_X_total <- data.frame(Res_A_X_total)
Res_A_X_total  <- rownames_to_column(Res_A_X_total, var = "ensembl_gene_id")

Res_A_X_total$sig <- ifelse(Res_A_X_total$padj <= 0.05, "yes", "no")

Res_A_X_total_0.05 <- subset(Res_A_X_total, padj <= 0.05)

# Export output files

write.csv(normalized_counts,"/home/2023/gdav/jacob.gonzalez.isa.gdav2023/final_project/04de_analysis/results/deseq2_normcounts.csv")
write.csv(Res_A_X_total, "/home/2023/gdav/jacob.gonzalez.isa.gdav2023/final_project/04de_analysis/results/deseq2_results_total.csv")
write.csv(Res_A_X_total_0.05, "/home/2023/gdav/jacob.gonzalez.isa.gdav2023/final_project/04de_analysis/results/deseq2_results_padj_0.05.csv")

cat("###WE'RE DONE!###\n")
