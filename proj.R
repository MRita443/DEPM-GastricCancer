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

###############################################
# 1. Download data from TCGA
###############################################

proj <- "TCGA-STAD"   # Stomach adenocarcinoma
data_dir <- "TCGA-STAD"

## --- Primary Tumor RNA-seq data --------------------------------
rna.query.C <- GDCquery(
  project      = proj,
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type   = "Primary Tumor"
)

GDCdownload(rna.query.C, directory = data_dir, method = "api", files.per.chunk = 10)
rna.data.C <- GDCprepare(rna.query.C, directory = data_dir)
rna.expr.data.C <- assay(rna.data.C)
genes.info <- as.data.frame(rowRanges(rna.data.C))

## --- Normal Tissue RNA-seq data --------------------------------
rna.query.N <- GDCquery(
  project      = proj,
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type   = "Solid Tissue Normal"
)

GDCdownload(rna.query.N, directory = data_dir, method = "api")
rna.data.N <- GDCprepare(rna.query.N, directory = data_dir)
rna.expr.data.N <- assay(rna.data.N)
genes.info2 <- as.data.frame(rowRanges(rna.data.N))

dim(rna.expr.data.C)
dim(rna.expr.data.N)

## Check gene annotation consistency
all(na.omit(genes.info2) == na.omit(genes.info))

###############################################
# 2. Data cleaning and sample alignment
###############################################

## Check for duplicates
ncol(rna.expr.data.N)
length(unique(substr(colnames(rna.expr.data.N), 1,12))) # No duplicates 
ncol(rna.expr.data.C)
length(unique(substr(colnames(rna.expr.data.C), 1,12))) # No duplicates

## Convert to data frames
expr.C <- as.data.frame(rna.expr.data.C)
expr.N <- as.data.frame(rna.expr.data.N)

## Shorten patient identifiers
short_names_C <- substr(colnames(expr.C), 1, 12)
short_names_N <- substr(colnames(expr.N), 1, 12)
colnames(expr.C) <- short_names_C
colnames(expr.N) <- short_names_N

## Identify matched samples (patients present in both groups)
matched <- intersect(short_names_N, short_names_C)
length(matched) # 33
setdiff(short_names_N, short_names_C) # 3 normal samples without cancer pair

## Subset to paired samples only
expr.C <- expr.C[, matched, drop = FALSE]
expr.N <- expr.N[, matched, drop = FALSE]

length(intersect(colnames(expr.N), colnames(expr.C)))
setdiff(colnames(expr.N), colnames(expr.C))

dim(expr.C)
dim(expr.N)

## Basic integrity checks
stopifnot(all(rownames(expr.C) == rownames(expr.N)))
stopifnot(!any(is.na(expr.C)), !any(is.na(expr.N)))
stopifnot(!any(is.nan(as.matrix(expr.C))), !any(is.nan(as.matrix(expr.N))))

###############################################
# 3. Normalization using DESeq2
###############################################

## Construct full count matrix
full.data <- cbind(expr.N, expr.C) # Normal first, Cancer after
full.data <- data.frame(full.data)
gene_ids <- rownames(full.data)

## Build condition metadata
condition <- factor(c(
  rep("Normal", ncol(expr.N)), # Repeats "Normal" for each normal sample
  rep("Tumor",  ncol(expr.C)) # Repeats "Tumor" for each tumor sample
))

metad <- data.frame(condition = condition)
rownames(metad) <- colnames(full.data)

## DESeq2 object
full.data.tidy <- cbind(gene = gene_ids, full.data)

dds <- DESeqDataSetFromMatrix(
  countData = full.data.tidy,
  colData   = metad,
  design    = ~ condition,
  tidy      = TRUE
)

## Filter low-count genes
keep <- rowSums(counts(dds) >= 10) >= (0.9 * ncol(expr.C)) # Over 10 counts on 90% of patients
dds <- dds[keep,]

## Normalize counts
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized = TRUE)

any(sum(rowSums(normalized_counts == 0) == ncol(expr.C))) # Verify no gene is 0 on all samples ASK: Should this check across all conjoint samples (66)?

## Split normalized data back into Normal and Tumor
n_norm <- ncol(expr.N)

filtr.expr.N <- as.data.frame(normalized_counts[, 1:n_norm])
filtr.expr.C <- as.data.frame(normalized_counts[, (n_norm + 1):ncol(normalized_counts)])

## Cancerous sample names were added a ".1" because they had the same names as the normal samples
## Remove unwanted ".1" suffixes
colnames(filtr.expr.C) <- substr(colnames(filtr.expr.C), 1, 12)