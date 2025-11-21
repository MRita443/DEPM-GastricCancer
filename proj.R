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

setwd('/home/yunbao/Bureau/Sapienza/DEPML/DEPM-GastricCancer')

# PART 1

# QUESTION 1
#1: Downloading data from the TCGA -------

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

#####################
# On récupère les données cliniques des associées aux patients
clinical.query<- GDCquery_clinic(project = proj, type = "clinical", save.csv = FALSE)
#write.csv(clinical.query, file = file.path(proj,paste(proj, "_clinical_data.txt",sep="")), row.names = FALSE, quote = FALSE)

View(clinical.query)
table(clinical.query$ajcc_pathologic_stage)

boxplot(age_at_index ~ ajcc_pathologic_stage, data = clinical.query,
        col = "gold", main = "Title", xlab = "", ylab= "age", las=2 )
########################################################################

View(rna.expr.data.N)
View(rna.expr.data.C)

dim(rna.expr.data.C)
dim(rna.expr.data.N)

#2: Data cleaning -----

#find duplicates 
ncol(rna.expr.data.N)
length(unique(substr(colnames(rna.expr.data.N), 1,12))) #no duplicates 
ncol(rna.expr.data.C)
length(unique(substr(colnames(rna.expr.data.C), 1,12))) #no duplicates

patients.C <- substr(colnames(rna.expr.data.C), 1,12)
sort(table(patients.C)) 
#keep patients with ony a sample
unique.patients.C <- names(which(table(patients.C) == 1)) 
#let's get their index in the list of patients
idx.unique.pats <- match(patients.C, substr(colnames(rna.expr.data.C), 1,12) )

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

filtr.expr.N <- as.data.frame(normalized_counts[, 1:ncol(expr.N)])
filtr.expr.C <- as.data.frame(normalized_counts[, (ncol(expr.N) + 1):ncol(normalized_counts)])

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


# QUESTION 2 :

# TODO : ask the professor if we need a set of genes of interest

#5: Differentially expressed genes (DEGs)

fc <-  log2(rowMeans(filtr.expr.C) / rowMeans(filtr.expr.N) ) 
names(fc) <- rownames(filtr.expr.C)
head(fc)

pval.fc <- sapply(1:nrow(filtr.expr.C), function(i) (t.test(as.numeric(filtr.expr.C[i,]), as.numeric(filtr.expr.N[i,]), paired = T ))$p.value)
pval.fc.fdr <- p.adjust(pval.fc, method="fdr")

expr.table <- data.frame(cbind(fc, pval.fc.fdr))
expr.table[,1] <- round(expr.table[,1],2)

# we apply the thresholds
deg.genes <- rownames(expr.table[abs(expr.table$fc) >= 1.2 & expr.table$pval.fc.fdr <=0.05,]) 
head(expr.table[deg.genes,], 10)
write.table(expr.table[deg.genes,], file = "DEG.csv", sep = ";")

#volcano plot
expr.table$diffexpressed <- "/";
expr.table$diffexpressed[expr.table$fc >= 1.2 & expr.table$pval.fc.fdr <= 0.05] <- "UP"
expr.table$diffexpressed[expr.table$fc <= -1.2 & expr.table$pval.fc.fdr <= 0.05] <- "DOWN"
head(expr.table, 5)

expr.table$diffexpressed <- as.factor(expr.table$diffexpressed)
summary(expr.table$diffexpressed)

ggplot(data=expr.table, aes(x=fc, y=-log10(pval.fc.fdr), col=diffexpressed))+  
  geom_point() +
  xlab("fold change (log2)") + 
  ylab("-log10 adjusted p-value") +
  geom_hline(yintercept=-log10(0.01), col="red")+
  geom_vline(xintercept=1.5, col="red")+
  geom_vline(xintercept=-1.5, col="red")

#cat(deg.genes , sep = "\n")


#6: Adjacency (correlations) matrices of co-expression networks

# CANCER NETWORK
# correlation between genes
cor.mat.c <- corr.test(t(filtr.expr.C), use="pairwise", method="spearman", adjust="fdr", ci=FALSE)
rho.c <- cor.mat.c$r # put the correlations into a matrix
diag(rho.c) <- 0
qval.c <- cor.mat.c$p #padj matrix
qval.c[lower.tri(qval.c)] <- t(qval.c)[lower.tri(qval.c)]
#retain only links with q <0.01
adj.mat.c <- rho.c * (qval.c <= 0.01)

# NORMAL NETWORK
cor.mat.n <- corr.test(t(filtr.expr.N), use="pairwise", method="spearman", adjust="fdr", ci=FALSE)
rho.n <- cor.mat.n$r
diag(rho.n) <- 0
qval.n <- cor.mat.n$p
qval.n[lower.tri(qval.n)] <- t(qval.n)[lower.tri(qval.n)]
adj.mat.n <- rho.n * (qval.n <= 0.01)

#7 : Co-expression networks

# CANCER NETWORK :
# TODO

# NORMAL NETWORK :
# TODO

