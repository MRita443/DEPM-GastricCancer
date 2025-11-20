library(BiocGenerics) 
library(DESeq2)
library(psych) 
library(NetworkToolbox)
library(ggplot2) 
library(GGally)
library(sna)
library(network)
library(TCGAbiolinks)
library(GenomicRanges)
library(SummarizedExperiment)
library(DT)
library(readr)
library(dplyr)

# PART 1 

#1: Downloading data from the TGCA -------

# create directory with the data
proj <- "TCGA-STAD" #stomach adenocarcinoma

#RNAseq data primary tumor
rna.query.C <- TCGAbiolinks::GDCquery(project = proj, data.category = "Transcriptome Profiling", 
                                      data.type = "Gene Expression Quantification",
                                      workflow.type = "STAR - Counts",
                                      sample.type = "Primary Tumor")

#download data and prepare it to use it
GDCdownload(query = rna.query.C, directory = "TCGA-STAD", method = "api", files.per.chunk = 10)
rna.data.C <- GDCprepare(rna.query.C, directory = "TCGA-STAD")
rna.expr.data.C <- assay(rna.data.C)

#inspect it
View(BiocGenerics::as.data.frame(rowRanges(rna.data.C)))
genes.info <- BiocGenerics::as.data.frame(rowRanges(rna.data.C))

#RNAseq control patients
rna.query.N <- GDCquery(project = proj, data.category = "Transcriptome Profiling", 
                        data.type = "Gene Expression Quantification", 
                        workflow.type = "STAR - Counts", 
                        sample.type = "Solid Tissue Normal")

#same 
GDCdownload(query = rna.query.N, directory = "TCGA-STAD", method = "api")
rna.data.N <- GDCprepare(rna.query.N, directory = "TCGA-STAD"  )
rna.expr.data.N <- assay(rna.data.N)
genes.info2 <- BiocGenerics::as.data.frame(rowRanges(rna.data.N))
all(na.omit(genes.info2) == na.omit(genes.info))

View(rna.expr.data.N)

dim(rna.expr.data.C)
dim(rna.expr.data.N)

#2: Data cleaning -----

#find duplicates 
ncol(rna.expr.data.N)
length(unique(substr(colnames(rna.expr.data.N), 1,12))) #no duplicates 
ncol(rna.expr.data.C)
length(unique(substr(colnames(rna.expr.data.C), 1,12))) #no duplicates

expr.C <- as.data.frame(rna.expr.data.C)
expr.N <- as.data.frame(rna.expr.data.N)

#let's rename patients in a shorter way
colnames(expr.C) <- substr(colnames(expr.C), 1,12)
colnames(expr.N) <- substr(colnames(expr.N), 1,12)

#align cancer and control
matched <- intersect(colnames(expr.N), colnames(expr.C))
length(matched) #33
setdiff(colnames(expr.N), colnames(expr.C))

# drop samples with no pair
expr.C <- expr.C[, matched, drop = FALSE]
expr.N <- expr.N[, matched, drop = FALSE]

length(intersect(colnames(expr.N), colnames(expr.C)))
setdiff(colnames(expr.N), colnames(expr.C))

dim(expr.C)
dim(expr.N)

#let's check the actual counts
typeof(expr.C[1,1]) #ok
any(is.na(expr.C)) #ok
any(is.nan(as.matrix(expr.C))) #ok

typeof(expr.N[1,1]) #ok
any(is.na(expr.N)) #ok
any(is.nan(as.matrix(expr.N))) #ok

#3: Normalizing data with Deseq2 ----- 

all(rownames(expr.C) == rownames(expr.N))
full.data <- cbind(expr.N, expr.C)

#full dataset cancer + control
dim(full.data)
full.data <- data.frame(full.data)

# prepare dataset. assign to each column normal vs cancerous condition

condition <- factor(
  c(
    rep("Normal", ncol(expr.N)),   # repeats "Normal" for each normal sample
    rep("Tumor", ncol(expr.C))  # repeats "Tumor" for each tumor sample
  )
)

metad <- data.frame(condition)
rownames(metad) <- colnames(full.data)
colnames(metad)[1] <- "condition"
metad[,1] <- as.factor(metad[,1])
full.data <- cbind(rownames(full.data), full.data)

dds <- DESeqDataSetFromMatrix(countData = full.data,
                              colData = metad,
                              design = ~ condition,
                              tidy=TRUE)

# remove low counts genes
keep <- rowSums(counts(dds) >= 10) >= (0.9 * ncol(expr.C)) # over 10 counts on 90% of patients
dds <- dds[keep,]
dim(counts(dds))

#normalize
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
sum(rowSums(normalized_counts == 0) == ncol(expr.C)) #no gene is 0 on all samples

filtr.expr.N <- as.data.frame(normalized_counts[, 1:ncol(expr.C)])
filtr.expr.C <- as.data.frame(normalized_counts[, (ncol(expr.C) + 1):ncol(normalized_counts)])

#cancerous sample names were added a ".1" in full.data because  
#they had the same names as the normal samples
colnames(filtr.expr.C) <- substr(colnames(filtr.expr.C), 1,12)

#set genes to names instead of IDs
# ASK: COMMENTED DUE TO DUPLICATE GENE NAMES
# genes.C <- intersect(rownames(filtr.expr.C), 
#                      genes.info[ , "gene_id"]   ) 
# 
# genes.N <- intersect(rownames(filtr.expr.N),  
#                      genes.info2[ , "gene_id"]   )  
# 
# setdiff(genes.C, genes.N)
# 
# length(genes.C)
# length(genes.N)
# 
# rownames(filtr.expr.N) <- genes.info2[genes.N, "gene_name"]
# rownames(filtr.expr.C) <- genes.info[genes.C, "gene_name"]
