###############################################################################
# DEPM Project – TCGA-STAD Gene Expression, DEG Identification and Network Analysis
###############################################################################

## --- Load Required Libraries --------------------------------------------------
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

## ============================================================================== 
## PART 1 – DATA ACQUISITION
## ============================================================================== 

# ------------------------------------------------------------------------------
# 1. Download and prepare TCGA-STAD RNA-seq data
# ------------------------------------------------------------------------------

proj <- "TCGA-STAD"

# --- Tumor samples -------------------------------------------------------------
rna.query.C <- GDCquery(
  project = proj,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor"
)

GDCdownload(rna.query.C, directory = proj, method = "api", files.per.chunk = 10)
rna.data.C <- GDCprepare(rna.query.C, directory = proj)
rna.expr.C  <- assay(rna.data.C)
genes.info  <- as.data.frame(rowRanges(rna.data.C))

# --- Normal samples ------------------------------------------------------------
rna.query.N <- GDCquery(
  project = proj,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = "Solid Tissue Normal"
)

GDCdownload(rna.query.N, directory = proj, method = "api")
rna.data.N <- GDCprepare(rna.query.N, directory = proj)
rna.expr.N  <- assay(rna.data.N)
genes.info2 <- as.data.frame(rowRanges(rna.data.N))

# gene order check
stopifnot(all(na.omit(genes.info2) == na.omit(genes.info)))

# --- Clinical data -------------------------------------------------------------
clinical <- GDCquery_clinic(project = proj, type = "clinical")
table(clinical$ajcc_pathologic_stage)

# Simple QC visualization
boxplot(age_at_index ~ ajcc_pathologic_stage, data = clinical,
        col = "gold", main = "Age vs Stage",
        xlab = "", ylab = "Age", las = 2)

## ============================================================================== 
## PART 2 – CLEANING AND ALIGNMENT
## ============================================================================== 

# ------------------------------------------------------------------------------
# 2. Basic QC and filtering
# ------------------------------------------------------------------------------

expr.C <- as.data.frame(rna.expr.C)
expr.N <- as.data.frame(rna.expr.N)

# unify column names to patient IDs
colnames(expr.C) <- substr(colnames(expr.C), 1, 12)
colnames(expr.N) <- substr(colnames(expr.N), 1, 12)

# keep matched tumor/normal patients
matched <- intersect(colnames(expr.N), colnames(expr.C))
expr.C <- expr.C[, matched, drop = FALSE]
expr.N <- expr.N[, matched, drop = FALSE]

stopifnot(all(colnames(expr.C) == colnames(expr.N)))

# QC NA checks
stopifnot(!any(is.na(expr.C)), !any(is.na(expr.N)))

## ============================================================================== 
## PART 3 – NORMALIZATION (DESeq2)
## ============================================================================== 

# ------------------------------------------------------------------------------
# 3. Normalize using DESeq2
# ------------------------------------------------------------------------------

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

###############################################################################
# QUESTION 2 : Differentially Expressed Genes (DEGs)
###############################################################################

# 5: Compute fold change (FC) for each gene
# FC = log2(mean expression in cancer / mean expression in normal)
fc <- log2(rowMeans(filtr.expr.C) / rowMeans(filtr.expr.N))
names(fc) <- rownames(filtr.expr.C)
head(fc)

# Compute paired t-test p-values gene-by-gene
# (same patients in tumor vs normal -> paired test)
pval.fc <- sapply(
  1:nrow(filtr.expr.C),
  function(i) {
    t.test(
      as.numeric(filtr.expr.C[i, ]),
      as.numeric(filtr.expr.N[i, ]),
      paired = TRUE
    )$p.value
  }
)

# Adjust p-values using FDR correction
pval.fc.fdr <- p.adjust(pval.fc, method = "fdr")

# Build table with FC + adjusted p-values
expr.table <- data.frame(cbind(fc, pval.fc.fdr))
expr.table[, 1] <- round(expr.table[, 1], 2)

# Select DEGs using FC and FDR thresholds
deg.genes <- rownames(
  expr.table[
    abs(expr.table$fc) >= 1.2 & expr.table$pval.fc.fdr <= 0.01,
  ]
)

# Inspect first 10 DEGs
head(expr.table[deg.genes, ], 10)

# Save DEG table
write.table(expr.table[deg.genes, ], file = "DEG.csv", sep = ";")

# ------------------------------------------------------------------------------
# Volcano plot preparation
# ------------------------------------------------------------------------------

expr.table$diffexpressed <- "/"
expr.table$diffexpressed[
  expr.table$fc >= 1.2 & expr.table$pval.fc.fdr <= 0.01
] <- "UP"
expr.table$diffexpressed[
  expr.table$fc <= -1.2 & expr.table$pval.fc.fdr <= 0.01
] <- "DOWN"

head(expr.table, 5)

expr.table$diffexpressed <- as.factor(expr.table$diffexpressed)
summary(expr.table$diffexpressed)

# Volcano plot
ggplot(
  data = expr.table,
  aes(x = fc, y = -log10(pval.fc.fdr), col = diffexpressed)
) +
  geom_point() +
  xlab("fold change (log2)") +
  ylab("-log10 adjusted p-value") +
  geom_hline(yintercept = -log10(0.01), col = "red") +
  geom_vline(xintercept = 1.2, col = "red") +
  geom_vline(xintercept = -1.2, col = "red")

# Extract cancer/normal matrices only for DEGs
filtr.expr.C.DEGs <- filtr.expr.C[rownames(filtr.expr.C) %in% deg.genes, ]
filtr.expr.N.DEGs <- filtr.expr.N[rownames(filtr.expr.N) %in% deg.genes, ]

# cat(deg.genes, sep="\n")


###############################################################################
# QUESTION 3 : Gene Co-expression Networks
###############################################################################

# 6: Build adjacency matrices from pairwise correlations
# ------------------------------------------------------

# ===================== CANCER NETWORK ===================== #

# Spearman correlation matrix (with FDR correction)
cor.mat.c <- corr.test(
  t(filtr.expr.C.DEGs),
  use = "pairwise",
  method = "spearman",
  adjust = "fdr",
  ci = FALSE
)

# Correlation matrix
rho.c <- cor.mat.c$r
diag(rho.c) <- 0  # remove self-correlation

# Corresponding adjusted p-values
qval.c <- cor.mat.c$p
qval.c[lower.tri(qval.c)] <- t(qval.c)[lower.tri(qval.c)]  # symmetrize matrix

# Start adjacency matrix from rho
adj.mat.c <- rho.c

# Apply correlation threshold |rho| ≥ 0.7
adj.mat.c <- adj.mat.c * (abs(rho.c) >= 0.7)

# Apply significance threshold (FDR)
adj.mat.c <- adj.mat.c * (qval.c <= 1e-3)

# Binary adjacency matrix
adj.bin.c <- adj.mat.c
adj.bin.c[abs(adj.mat.c) > 0] <- 1


# ===================== NORMAL NETWORK ===================== #

cor.mat.n <- corr.test(
  t(filtr.expr.N.DEGs),
  use = "pairwise",
  method = "spearman",
  adjust = "fdr",
  ci = FALSE
)

rho.n <- cor.mat.n$r
diag(rho.n) <- 0

qval.n <- cor.mat.n$p
qval.n[lower.tri(qval.n)] <- t(qval.n)[lower.tri(qval.n)]

adj.mat.n <- rho.n

# Slightly stricter correlation threshold for normal samples
adj.mat.n <- adj.mat.n * (abs(rho.n) >= 0.85)

# Apply p-value threshold
adj.mat.n <- adj.mat.n * (qval.n <= 1e-5)

adj.bin.n <- adj.mat.n
adj.bin.n[abs(adj.mat.n) > 0] <- 1


###############################################################################
# 7 : Network Construction and Analysis (Scale-free check, hubs, etc.)
###############################################################################

# ===================== CANCER NETWORK ANALYSIS ===================== #

# Build network object
net.c <- network(
  adj.mat.c,
  matrix.type = "adjacency",
  ignore.eval = FALSE,
  names.eval = "weights",
  directed = FALSE
)

# --- Basic network statistics ---
network.density(net.c)
network.size(net.c)
network.edgecount(net.c)
clustcoeff(adj.mat.c, weighted = FALSE)$CC

# Count edges
sum(adj.mat.c != 0) / 2
sum(adj.mat.c > 0) / 2   # positive correlations
sum(adj.mat.c < 0) / 2   # negative correlations

# Degrees
degree.c <- rowSums(adj.mat.c != 0)
names(degree.c) <- rownames(adj.mat.c)
degree.c <- sort(degree.c, decreasing = TRUE)

head(degree.c, 10)
sum(degree.c == 0)  # number of isolated nodes

# Scale-free check
hist(degree.c)

# --- Identify hubs (top 5% degree) ---
x <- quantile(degree.c[degree.c > 0], 0.95)
x  # threshold value
abline(v = x, col = "red")

hubs.c <- degree.c[degree.c >= x]
names(hubs.c)

write.table(hubs.c, file = "hubs_c.csv", sep = ";")


# --- Annotate and plot network ---
net.c %v% "type"  <- ifelse(network.vertex.names(net.c) %in% names(hubs.c), "hub", "non-hub")
net.c %v% "color" <- ifelse(net.c %v% "type" == "hub", "tomato", "deepskyblue3")

set.edge.attribute(
  net.c, "edgecolor",
  ifelse(net.c %e% "weights" > 0, "red", "blue")
)

ggnet2(
  net.c, color = "color", alpha = 0.7, size = 2,
  edge.color = "edgecolor", edge.alpha = 1, edge.size = 0.15
) + guides(size = "none")


# ===================== NORMAL NETWORK ANALYSIS ===================== #

net.n <- network(
  adj.mat.n,
  matrix.type = "adjacency",
  ignore.eval = FALSE,
  names.eval = "weights"
)

network.density(net.n)
network.size(net.n)
network.edgecount(net.n)
clustcoeff(adj.mat.n, weighted = FALSE)$CC

sum(adj.mat.n != 0) / 2
sum(adj.mat.n > 0) / 2
sum(adj.mat.n < 0) / 2

degree.n <- rowSums(adj.mat.n != 0)
names(degree.n) <- rownames(adj.mat.n)
degree.n <- sort(degree.n, decreasing = TRUE)

head(degree.n, 10)
sum(degree.n == 0)

hist(degree.n)

# Identify hubs
y <- quantile(degree.n[degree.n > 0], 0.95)
y
hist(degree.n)
abline(v = y, col = "red")

hubs.n <- degree.n[degree.n >= y]
names(hubs.n)

write.table(hubs.n, file = "hubs_n.csv", sep = ";")


# Annotate and plot
net.n %v% "type"  <- ifelse(network.vertex.names(net.n) %in% names(hubs.n), "hub", "non-hub")
net.n %v% "color" <- ifelse(net.n %v% "type" == "hub", "tomato", "deepskyblue3")

set.edge.attribute(
  net.n, "edgecolor",
  ifelse(net.n %e% "weights" > 0, "red", "blue")
)

ggnet2(
  net.n, color = "color", alpha = 0.7, size = 2,
  edge.color = "edgecolor", edge.alpha = 1, edge.size = 0.15
) + guides(size = "none")


# --- Overlap between cancer and normal hubs ---
length(names(hubs.c))          # 59
length(names(hubs.n))          # 56
length(intersect(names(hubs.c), names(hubs.n)))  # 18

hubs_intersect_coexpr <- intersect(names(hubs.c), names(hubs.n))
hubs_intersect_coexpr

gained_hubs_coexpr <- setdiff(names(hubs.c), names(hubs.n))
gained_hubs_coexpr

lost_hubs_coexpr <- setdiff(names(hubs.n), names(hubs.c))
lost_hubs_coexpr

# TODO in report: KEGG enrichment for gained/lost hubs


###############################################################################
# OPTIONAL VISUALIZATION – HUB SUBNETWORK
###############################################################################

hubs.c
hubs.c.ids <- vector("integer", length(hubs.c))

for (i in 1:length(hubs.c)) {
  hubs.c.ids[i] <- match(names(hubs.c)[i], rownames(adj.mat.c))
}

# Identify neighbors of each cancer hub
hubs.c.neigh <- c()
for (f in hubs.c.ids) {
  hubs.c.neigh <- append(hubs.c.neigh, get.neighborhood(net.c, f))
}

hubs.c.neigh <- unique(hubs.c.neigh)
hubs.c.neigh

hubs.c.neigh.names <- rownames(adj.mat.c[hubs.c.neigh, ])
subnet <- unique(c(names(hubs.c), hubs.c.neigh.names))

# Build subnetwork adjacency matrix
hub.c.adj <- adj.mat.c[subnet, subnet]

head(rownames(hub.c.adj))
head(colnames(hub.c.adj))

net.hub <- network(
  hub.c.adj, matrix.type = "adjacency",
  ignore.eval = FALSE, names.eval = "weights"
)

network.density(net.hub)

sum(hub.c.adj > 0)
sum(hub.c.adj < 0)

net.hub %v% "type"  <- ifelse(network.vertex.names(net.hub) %in% names(hubs.c), "hub", "non-hub")
net.hub %v% "color" <- ifelse(net.hub %v% "type" == "non-hub", "deepskyblue3", "tomato")

set.edge.attribute(
  net.hub, "ecolor",
  ifelse(net.hub %e% "weights" > 0, "red", "blue")
)

ggnet2(
  net.hub, color = "color", alpha = 0.9, size = 2,
  edge.color = "ecolor", edge.alpha = 0.9, edge.size = 0.15,
  node.label = names(hubs.c), label.color = "black", label.size = 3
) + guides(size = "none")


###############################################################################
# QUESTION 4 : Differential Co-expressed Network
###############################################################################

# Fisher z-transform for cancer and normal adjacency matrices
z.mat.C <- 0.5 * log((1 + adj.mat.c) / (1 - adj.mat.c))
z.mat.N <- 0.5 * log((1 + adj.mat.n) / (1 - adj.mat.n))

# Compute Z-scores for differential correlations
Z.mat <- (z.mat.C - z.mat.N) /
  sqrt((1 / (length(filtr.expr.C.DEGs) - 3)) +
         (1 / (length(filtr.expr.N.DEGs) - 3)))

# Threshold |Z| < 3 -> keep only significant changes
Z.mat[abs(Z.mat) < 3] <- 0

# Binary differential adjacency matrix
Z.bin <- (abs(Z.mat) != 0) * 1


# Build differential network
net.z <- network(
  Z.bin,
  matrix.type = "adjacency",
  ignore.eval = FALSE,
  names.eval = "weights",
  directed = FALSE
)

# --- Network analysis ---
network.density(net.z)
network.size(net.z)
network.edgecount(net.z)
clustcoeff(Z.bin, weighted = FALSE)$CC

degree.z <- rowSums(Z.bin != 0)
names(degree.z) <- rownames(Z.bin)
degree.z <- sort(degree.z, decreasing = TRUE)

head(degree.z, 10)
sum(degree.z == 0)

hist(degree.z)

# --- Identify differential hubs ---
x.z <- quantile(degree.z[degree.z > 0], 0.95)
x.z
hist(degree.z)
abline(v = x.z, col = "red")

hubs.z <- degree.z[degree.z >= x.z]
names(hubs.z)

write.table(hubs.c, file = "hubs_z.csv", sep = ";")


# Annotate and plot differential network
net.z %v% "type"  <- ifelse(network.vertex.names(net.z) %in% names(hubs.z), "hub", "non-hub")
net.z %v% "color" <- ifelse(net.z %v% "type" == "hub", "tomato", "deepskyblue3")

network::set.edge.attribute(
  net.z, "edgecolor",
  ifelse(net.z %e% "weights" > 0, "red", "blue")
)

ggnet2(
  net.z, color = "color", alpha = 0.7, size = 2,
  edge.color = "edgecolor", edge.alpha = 1, edge.size = 0.15
) + guides(size = "none")


# --- Hub comparisons ---
length(names(hubs.c))  # 59
length(names(hubs.n))  # 56
length(names(hubs.z))  # 70

length(intersect(names(hubs.c), names(hubs.n)))  # 18
length(intersect(names(hubs.z), names(hubs.n)))  # 49
length(intersect(names(hubs.c), names(hubs.z)))  # 21
length(intersect(names(hubs.z), hubs_intersect_coexpr))  # 14

hubs_intersect_diffcoexpr_and_coexpr <- intersect(names(hubs.z), hubs_intersect_coexpr)
hubs_intersect_diffcoexpr_and_coexpr

gained_hubs_diffcoexpr <- setdiff(names(hubs.z), names(hubs.c))
gained_hubs_diffcoexpr <- setdiff(gained_hubs_diffcoexpr, names(hubs.n))
gained_hubs_diffcoexpr


###############################################################################
# QUESTION 5 : Patient Similarity Network (PSN)
###############################################################################

# (Code continues…)
=======
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
>>>>>>> ad4960cb3627acdc136cd298730f78408fd64716
