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
library(DGCA) # package for differential co-expression network
library(igraph)
library(maftools)

set.seed(123)
# NOTE: Set your working directory here if needed
# setwd('/your/path/here')

###############################################
# 1. Download data from TCGA
###############################################

proj <- "TCGA-STAD"   # Stomach adenocarcinoma
data_dir <- "TCGA-STAD"

## --- Primary Tumor RNA-seq data --------------------------------
rna.query.C <- GDCquery(
  project       = proj,
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
  project       = proj,
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

any(sum(rowSums(normalized_counts == 0) == ncol(expr.C))) 

## Split normalized data back into Normal and Tumor
n_norm <- ncol(expr.N)

filtr.expr.N <- as.data.frame(normalized_counts[, 1:n_norm])
filtr.expr.C <- as.data.frame(normalized_counts[, (n_norm + 1):ncol(normalized_counts)])

## Cancerous sample names were added a ".1" because they had the same names as the normal samples
## Remove unwanted ".1" suffixes
colnames(filtr.expr.C) <- substr(colnames(filtr.expr.C), 1, 12)

###############################################################################
# QUESTION 2 : Differentially Expressed Genes (DEGs)
###############################################################################

# 5: Compute fold change (FC) for each gene
# FC = log2(mean expression in cancer / mean expression in normal)
fc <- log2(rowMeans(filtr.expr.C) / rowMeans(filtr.expr.N))
names(fc) <- rownames(filtr.expr.C)
head(fc)

# Compute paired t-test p-values gene-by-gene
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

network::set.edge.attribute(
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

network::set.edge.attribute(
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


###############################################################################
# OPTIONAL VISUALIZATION – HUB SUBNETWORK
###############################################################################

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
hubs.c.neigh.names <- rownames(adj.mat.c[hubs.c.neigh, ])
subnet <- unique(c(names(hubs.c), hubs.c.neigh.names))

# Build subnetwork adjacency matrix
hub.c.adj <- adj.mat.c[subnet, subnet]

net.hub <- network(
  hub.c.adj, matrix.type = "adjacency",
  ignore.eval = FALSE, names.eval = "weights"
)

network.density(net.hub)

net.hub %v% "type"  <- ifelse(network.vertex.names(net.hub) %in% names(hubs.c), "hub", "non-hub")
net.hub %v% "color" <- ifelse(net.hub %v% "type" == "non-hub", "deepskyblue3", "tomato")

network::set.edge.attribute(
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

write.table(hubs.z, file = "hubs_z.csv", sep = ";")


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
# QUESTION 5 : Patient Similarity Network (PSN) & Fusion
###############################################################################

# Helper function to run the Python Bridge
run_python_community <- function(matrix_data) {
  # 1. Save matrix to CSV format expected by python script (sep=";", dec=",")
  # NOTE: Update path if not running locally with this specific structure
  write.table(matrix_data, file = "input-matrix.csv", 
              sep = ";", dec = ",", row.names = TRUE, col.names = NA)
  
  # 2. Run the Python script
  # NOTE: Ensure 'python3' and 'bctpy' are installed in your terminal
  # Update path to your python executable if necessary
  system("python3 btc-community.py input-matrix.csv")
  
  # 3. Read the output
  if(file.exists("output.txt")){
    comm <- read.table("output.txt", header = FALSE)
    return(as.factor(comm$V1))
  } else {
    stop("Python script did not generate output.txt")
  }
}

# ---------------------------------------------------------
# TASK 5a & 5b: Expression PSN + Community Detection
# ---------------------------------------------------------

# 1. Transpose: Patients as Rows
expr.C.transposed <- t(filtr.expr.C.DEGs)

# 2. Compute Similarity (Spearman)
# We use absolute correlation for similarity strength
W_expr <- abs(cor(expr.C.transposed, method = "spearman"))

# 3. Thresholding (Optional but recommended for Louvain stability)
# Keep only strong connections (e.g., > 0.6) to reduce noise
W_expr_clean <- W_expr
W_expr_clean[W_expr_clean < 0.6] <- 0 

# 4. Run Python Community Detection
print("Running Python for Expression Communities...")
comm_expr <- run_python_community(W_expr_clean)

# 5. Visualization (Using ggnet2 style from template)
# Create network object
net.psn.expr <- network(W_expr_clean, matrix.type = "adjacency", directed = FALSE, ignore.eval = FALSE, names.eval = "weights")

# Assign communities to nodes
net.psn.expr %v% "community" = as.character(comm_expr)

# Plot
ggnet2(net.psn.expr, color = "community", size = 3, 
       edge.alpha = 0.2, edge.size = 0.1,
       main = "PSN: Gene Expression (Cancer)") +
  guides(color = "none")


# ---------------------------------------------------------
# TASK 5c: Similarity Network Fusion (Expr + Mutation)
# ---------------------------------------------------------

# --- Step 1 & 2: Get Mutation Data & Matrix ---
# Ensuring proper naming
mut.query <- GDCquery(
  project = proj, 
  data.category = "Simple Nucleotide Variation", 
  access = "open", 
  data.type = "Masked Somatic Mutation", 
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
  # REMOVED: directory = data_dir (This argument is not valid for GDCquery)
)

# 2. Prepare the data (Loads it from your hard drive)
GDCdownload(mut.query, directory = data_dir) 
maf <- GDCprepare(mut.query, directory = data_dir)

# 3. Process into a Binary Matrix
maf.obj <- read.maf(maf)
mut.matrix.count <- mutCountMatrix(maf.obj)

# Create the 'mut.matrix.binary' object
mut.matrix.binary <- (mut.matrix.count > 0) * 1 
mut.matrix.binary <- t(mut.matrix.binary) # Transpose so Patients are Rows

rownames(mut.matrix.binary) <- substr(rownames(mut.matrix.binary), 1, 12)

# =========================================================
# CORRECTED STEP 3: Align Patients (Fixing Separators)
# =========================================================

# 1. Fix separators in Expression Data (Replace '.' with '-')
rownames(expr.C.transposed) <- gsub("\\.", "-", rownames(expr.C.transposed))

# 2. Ensure Mutation Data is also formatted (Trim to 12 chars if needed)
# (Your diagnostics show they are already 12 chars, but this is safe to keep)
rownames(mut.matrix.binary) <- substr(rownames(mut.matrix.binary), 1, 12)

# 3. Now Find Intersection
common.patients <- intersect(rownames(expr.C.transposed), rownames(mut.matrix.binary))
print(paste("Number of common patients:", length(common.patients)))

# 4. Stop if still empty (Safety check)
if(length(common.patients) == 0) {
  stop("No common patients! Please check if the patient IDs (e.g. TCGA-XX-XXXX) actually overlap between the two datasets.")
}

# 5. Subset both matrices to shared patients
dat_expr_sub <- expr.C.transposed[common.patients, ]
dat_mut_sub  <- mut.matrix.binary[common.patients, ]

print("Step 3 Complete: Patients aligned successfully.")

# --- Step 4: Compute Similarities ---

# 1. Expression Similarity (Patient vs Patient)
W_expr_sub <- abs(cor(t(dat_expr_sub), method = "spearman")) 

# Normalize to 0-1
W_expr_sub <- (W_expr_sub - min(W_expr_sub)) / (max(W_expr_sub) - min(W_expr_sub))

# 2. Mutation Similarity (Jaccard)
# Jaccard Distance = 1 - Similarity
dist_mut <- dist(dat_mut_sub, method = "binary") # binary = Jaccard distance
W_mut <- 1 - as.matrix(dist_mut)

# HANDLE NA: If a patient has 0 mutations, Jaccard might return NaN (0/0). 
# We assume 0 similarity in that case.
W_mut[is.na(W_mut)] <- 0

# Normalize to 0-1 (Check if max > min to avoid division by zero)
if(max(W_mut) > min(W_mut)){
  W_mut <- (W_mut - min(W_mut)) / (max(W_mut) - min(W_mut))
}

# --- Step 5: FUSION (Mean Fusion - Old Way) ---
# Now both matrices should be Patients x Patients (e.g., 30x30)
if(all(dim(W_expr_sub) == dim(W_mut))) {
  W_fused <- (W_expr_sub + W_mut) / 2
  
  diag(W_fused) <- 0
  
  W_fused[W_fused < 0.6] <- 0
  
  print("Fusion successful!")
} else {
  stop("Dimensions still do not match!")
}

# --- Step 6: Community Detection on Fused Network ---
# Use the function you defined earlier
comm_fused <- run_python_community(W_fused)

# --- Step 7: Visualizing Fused Network ---
net.psn.fused <- network(W_fused, matrix.type = "adjacency", 
                         directed = FALSE, ignore.eval = FALSE, names.eval = "weights")

net.psn.fused %v% "community" = as.character(comm_fused)

ggnet2(net.psn.fused, color = "community", size = 4, 
       edge.alpha = 0.2, edge.size = 0.1,
       node.label = TRUE, label.size = 2.5,
       main = "PSN: Fused (Expression + Mutation)") +
  guides(color = "none")

print("Expression Clusters:"); table(comm_expr)
print("Fused Clusters:"); table(comm_fused)