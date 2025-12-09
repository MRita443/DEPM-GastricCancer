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
library(igraph)
library(TCGAbiolinks)
library(maftools)
library(network)
library(GGally)

set.seed(123)
setwd('/home/rafal/Documents/studia magisterskie/semestr_2/DEPM-GastricCancer')

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

#5: Differentially expressed genes (DEGs)
# compute fc
fc <-  log2(rowMeans(filtr.expr.C) / rowMeans(filtr.expr.N) ) 
names(fc) <- rownames(filtr.expr.C)
head(fc)

# compute p-value for threshold
pval.fc <- sapply(1:nrow(filtr.expr.C), function(i) (t.test(as.numeric(filtr.expr.C[i,]), as.numeric(filtr.expr.N[i,]), paired = T ))$p.value)
# we adjust p-value with fdr method to control false positive
pval.fc.fdr <- p.adjust(pval.fc, method="fdr")
# table to resume
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
  geom_hline(yintercept=-log10(0.05), col="red")+
  geom_vline(xintercept=1.2, col="red")+
  geom_vline(xintercept=-1.2, col="red")

filtr.expr.C.DEGs <- filtr.expr.C[rownames(filtr.expr.C) %in% deg.genes, ]
filtr.expr.N.DEGs <- filtr.expr.N[rownames(filtr.expr.N) %in% deg.genes, ]

#cat(deg.genes , sep = "\n")

# Question 3

#6: Adjacency (correlations) matrices of co-expression networks

# CANCER NETWORK
# correlation between genes
cor.mat.c <- corr.test(t(filtr.expr.C.DEGs), use="pairwise", method="spearman", adjust="fdr", ci=FALSE)
# rho.c : matrix containing the correlations
rho.c <- cor.mat.c$r
diag(rho.c) <- 0 # put 0 on the diagonal of the correlations matrix
# qval.c : matrix containing p-values of the correlation matrix
qval.c <- cor.mat.c$p
qval.c[lower.tri(qval.c)] <- t(qval.c)[lower.tri(qval.c)] #(qvals are reported on the upper triangle only to have the matrix symetric)
#retain only links with q <0.01 : not asked in the project
adj.mat.c <- rho.c * (qval.c <= 0.01)
# keep only |correlation| > 0.7 : asked in the project
adj.mat.c  <- adj.mat.c * (abs(rho.c) >= 0.7)
adj.bin.c  <- adj.mat.c * 1 # get a binary version of the same matrix
# check if this works

# NORMAL NETWORK
cor.mat.n <- corr.test(t(filtr.expr.N.DEGs), use="pairwise", method="spearman", adjust="fdr", ci=FALSE)
rho.n <- cor.mat.n$r
diag(rho.n) <- 0
qval.n <- cor.mat.n$p
qval.n[lower.tri(qval.n)] <- t(qval.n)[lower.tri(qval.n)]
adj.mat.n <- rho.n * (qval.n <= 0.01)
adj.mat.n  <- adj.mat.n * (abs(rho.n) >= 0.7)
adj.bin.n  <- adj.mat.n * 1 # get a binary version of the same matrix

#7 : Co-expression networks
# ANALYSIS PART
# Aims to check if network is scale-free or not

# CANCER NETWORK :
#build the network
net.c <- network(adj.mat.c, matrix.type="adjacency",ignore.eval = FALSE, names.eval = "weights", directed = F)
### Check
network.density(net.c) #how much is the network connected
network.size(net.c) #n nodes
network.edgecount(net.c)  #n of edges
###
clustcoeff(adj.mat.c, weighted = FALSE)$CC # clustcoeff measures how connected are my neighbors between them
### Check
sum(adj.mat.c != 0) / 2
#how many positive/negative correlations? 
sum(adj.mat.c > 0) / 2
sum(adj.mat.c < 0) / 2
###
degree.c <- rowSums(adj.mat.c != 0)
names(degree.c) <- rownames(adj.mat.c)
degree.c <- sort(degree.c, decreasing = T)
head(degree.c,10)
sum(degree.c == 0) # check how many unconnected nodes
# Plot the histogram to see if we have a scale free network or not
hist(degree.c)
# Now we want find the hubs (5% of nodes with highest degree values)
x <- quantile(degree.c[degree.c>0],0.95) #top 5% of nodes
x
hist(degree.c)
abline(v=x, col="red")
hubs.c <- degree.c[degree.c>=x]
names(hubs.c) 
#let's enrich them
write.table(hubs.c, file = "hubs.csv", sep = ";")
# Anotate the hubs
net.c %v% "type" = ifelse(network.vertex.names(net.c) %in% names(hubs.c),"hub", "non-hub")
net.c %v% "color" = ifelse(net.c %v% "type" == "hub", "tomato", "deepskyblue3")
network::set.edge.attribute(net.c, "edgecolor", ifelse(net.c %e% "weights" > 0, "red", "blue"))
# Visualizing
ggnet2(net.c, color = "color", alpha = 0.7, size = 2,  #mode= c("x","y"),
       edge.color = "edgecolor", edge.alpha = 1, edge.size = 0.15)+
  guides(size = "none")
##### CHECK if needed define adj.mat.c with good p-value threshold
#this is extremely dense... what if we lower the pval threshold?
# TODO : check good threshold for p-value
adj.mat.c1 <- rho.c * (qval.c <= 1e-3)
adj.mat.c2 <- rho.c * (qval.c <= 1e-4) #too much?
#might be useful to look at negative and postive edges separately:
net.c1 <- network(adj.mat.c1* (adj.mat.c1 > 0), matrix.type="adjacency",ignore.eval = FALSE, names.eval = "weights")
ggnet2(net.c1, color = "deepskyblue3", alpha = 0.7, size = 2, 
       edge.color = "red", edge.alpha = 1, edge.size = 0.15)+
  guides(size = "none") 
net.c2 <- network(adj.mat.c2* (adj.mat.c2 < 0),
                  matrix.type="adjacency",ignore.eval = FALSE, names.eval = "weights")
ggnet2(net.c2, color = "deepskyblue3", alpha = 0.7, size = 2, 
       edge.color = "blue", edge.alpha = 1, edge.size = 0.15)+
  guides(size = "none") 
#####

# NORMAL NETWORK :
net.n <- network(adj.mat.n, matrix.type="adjacency",ignore.eval = FALSE, names.eval = "weights")
### Check
network.density(net.n)
network.size(net.n)
network.edgecount(net.n)
###
clustcoeff(adj.mat.n, weighted = FALSE)$CC
### Check
sum(adj.mat.n != 0) /2
sum(adj.mat.n > 0) /2
sum(adj.mat.n < 0) /2
###
degree.n <- rowSums(adj.mat.n != 0)
names(degree.n) <- rownames(adj.mat.n)
degree.n <- sort(degree.n, decreasing = T)
head(degree.n,10)
sum(degree.n == 0) #unconnected nodes 
hist(degree.n)
y <- quantile(degree.n[degree.n>0],0.95) 
y
hist(degree.n)
abline(v=y, col="red")
hubs.n <- degree.n[degree.n>=y]
names(hubs.n)
net.n %v% "type" = ifelse(network.vertex.names(net.n) %in% names(hubs.n),"hub", "non-hub")
net.n %v% "color" = ifelse(net.n %v% "type" == "hub", "tomato", "deepskyblue3")

network::set.edge.attribute(net.n, "edgecolor", ifelse(net.n %e% "weights" > 0, "red", "blue"))
# Visualizing
ggnet2(net.n, color = "color", alpha = 0.7, size = 2,
       edge.color = "edgecolor", edge.alpha = 1, edge.size = 0.15)+
  guides(size = "none") 

# See which genes are hubs both in cancer and normal tissue
intersect(names(hubs.c), names(hubs.n))
# TODO : identify the hubs selectively characterizing each network

##### CHECK if needed define adj.mat.n with good p-value threshold
# TODO : check good threshold for p-value
adj.mat.n1 <- rho.n * (qval.n <= 1e-3)
adj.mat.n2 <- rho.n * (qval.n <= 1e-4) #too much?
##### 


# NOT ASKED IN THE PROJECT, BUT USEFUL FOR VISUALIZATION
####
#8: Plotting the hub subnetwork -----

hubs.c
hubs.c.ids <- vector("integer",length(hubs.c))
for (i in 1:length(hubs.c)){hubs.c.ids[i] <- match(names(hubs.c)[i],rownames(adj.mat.c))}
hubs.c.ids

#identifying the neighborhood
hubs.c.neigh <- c()
for (f in hubs.c.ids){
  hubs.c.neigh <- append(hubs.c.neigh, get.neighborhood(net.c, f))
}

hubs.c.neigh <- unique(hubs.c.neigh)
hubs.c.neigh
hubs.c.neigh.names <- rownames(adj.mat.c[hubs.c.neigh,])
subnet <- unique(c(names(hubs.c), hubs.c.neigh.names))

#creating the subnetwork
hub.c.adj <- adj.mat.c[subnet, subnet]

head(rownames(hub.c.adj))
head(colnames(hub.c.adj))

net.hub <- network(hub.c.adj, matrix.type="adjacency",ignore.eval = FALSE, names.eval = "weights")
network.density(net.hub)

sum(hub.c.adj > 0 )
sum(hub.c.adj < 0)

net.hub %v% "type" = ifelse(network.vertex.names(net.hub) %in% names(hubs.c),"hub", "non-hub")
net.hub %v% "color" = ifelse(net.hub %v% "type" == "non-hub", "deepskyblue3", "tomato")
network::set.edge.attribute(net.hub, "ecolor", ifelse(net.hub %e% "weights" > 0, "red", "blue"))

ggnet2(net.hub,  color = "color",alpha = 0.9, size = 2, 
       edge.color = "ecolor", edge.alpha = 0.9,  edge.size = 0.15, 
       node.label = names(hubs.c), label.color = "black", label.size = 3)+
  guides(size = "none") 
####
####

# Question 4 : Differential Co-expressed Network
# We chose to use adj.mat.c which is the correlation matrix after filtering with the p-value the non-significative correlations
# Get z for each conditions (Cancer and Normal)
z.mat.C <- 0.5 * log((1+adj.mat.c)/(1-adj.mat.c))
z.mat.N <- 0.5 * log((1+adj.mat.n)/(1-adj.mat.n))
# Get Z-scores
Z.mat <- (z.mat.C - z.mat.N)/sqrt((1/(length(filtr.expr.C.DEGs)-3)) + (1/(length(filtr.expr.N.DEGs)-3)))
Z.mat[abs(Z.mat) < 3] <- 0 # apply threshold
# We make it binary
Z.bin <- (abs(Z.mat) != 0) * 1 # any value != 0 takes 1 as value

# ANALYSIS PART (exactly based on analysis question 3)
#build the network
net.z <- network(Z.bin, matrix.type="adjacency",ignore.eval = FALSE, names.eval = "weights", directed = F)
### Check
network.density(net.z) #how much is the network connected
network.size(net.z) #n nodes
network.edgecount(net.z)  #n of edges
###
clustcoeff(Z.bin, weighted = FALSE)$CC # clustcoeff measures how connected are my neighbors between them

degree.z <- rowSums(Z.bin != 0)
names(degree.z) <- rownames(Z.bin)
degree.z <- sort(degree.z, decreasing = T)
head(degree.z,10)
sum(degree.z == 0) # check how many unconnected nodes
# Plot the histogram to see if we have a scale free network or not
hist(degree.z)
# Now we want find the hubs (5% of nodes with highest degree values)
x.z <- quantile(degree.z[degree.z>0],0.95) #top 5% of nodes
x.z
hist(degree.z)
abline(v=x.z, col="red")
hubs.z <- degree.z[degree.z>=x.z]
names(hubs.z) 
#let's enrich them
#write.table(hubs.c, file = "hubs.csv", sep = ";")
# Anotate the hubs
net.z %v% "type" = ifelse(network.vertex.names(net.z) %in% names(hubs.z),"hub", "non-hub")
net.z %v% "color" = ifelse(net.z %v% "type" == "hub", "tomato", "deepskyblue3")
network::set.edge.attribute(net.z, "edgecolor", ifelse(net.z %e% "weights" > 0, "red", "blue"))
# Visualizing
ggnet2(net.z, color = "color", alpha = 0.7, size = 2,  #mode= c("x","y"),
       edge.color = "edgecolor", edge.alpha = 1, edge.size = 0.15)+
  guides(size = "none")

# TODO : compare the identified hubs set with those obtained in task 3.

# Question 5 : Patient Similarity Network (PSN)
# TODO








# =========================================================
# QUESTION 5: Patient Similarity Network (PSN) & Fusion
# =========================================================


# Helper function to run the Python Bridge
run_python_community <- function(matrix_data) {
  # 1. Save matrix to CSV format expected by python script (sep=";", dec=",")
  write.table(matrix_data, file = "/home/rafal/Documents/studia magisterskie/semestr_2/DEPM-GastricCancer/input-matrix.csv", 
              sep = ";", dec = ",", row.names = TRUE, col.names = NA)
  
  # 2. Run the Python script
  # Ensure 'python3' and 'bctpy' are installed in your terminal
  system("cd '/home/rafal/Documents/studia magisterskie/semestr_2/DEPM-GastricCancer/' && /home/rafal/anaconda3/bin/python3 btc-community.py input-matrix.csv")
  
  # 3. Read the output
  if(file.exists("/home/rafal/Documents/studia magisterskie/semestr_2/DEPM-GastricCancer/output.txt")){
    comm <- read.table("/home/rafal/Documents/studia magisterskie/semestr_2/DEPM-GastricCancer/output.txt", header = FALSE)
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

# --- Step 1 & 2: Get Mutation Data & Matrix (Already done by you) ---
# Ensuring proper naming
mut.query <- GDCquery(
  project = "TCGA-STAD", 
  data.category = "Simple Nucleotide Variation", 
  access = "open", 
  data.type = "Masked Somatic Mutation", 
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)

# 2. Prepare the data (Loads it from your hard drive)
GDCdownload(mut.query) 
maf <- GDCprepare(mut.query)

# 3. Process into a Binary Matrix
maf.obj <- read.maf(maf)
mut.matrix.count <- mutCountMatrix(maf.obj)

# Create the 'mut.matrix.binary' object
mut.matrix.binary <- (mut.matrix.count > 0) * 1 
mut.matrix.binary <- t(mut.matrix.binary) # Transpose so Patients are Rows



rownames(mut.matrix.binary) <- substr(rownames(mut.matrix.binary), 1, 12)

# --- Step 3: Align Patients ---
common.patients <- intersect(rownames(expr.C.transposed), rownames(mut.matrix.binary))
print(paste("Number of common patients:", length(common.patients)))

if(length(common.patients) == 0) {
  stop("No common patients found between RNA-seq and Mutation data! Check your barcode naming.")
}

# Subset both matrices to shared patients
dat_expr_sub <- expr.C.transposed[common.patients, ]
dat_mut_sub  <- mut.matrix.binary[common.patients, ]

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

# --- Step 5: FUSION (Mean Fusion) ---
# Now both matrices should be Patients x Patients (e.g., 30x30)
if(all(dim(W_expr_sub) == dim(W_mut))) {
  W_fused <- (W_expr_sub + W_mut) / 2
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

# ------------------------
# LOAD REQUIRED LIBRARIES
# ------------------------
library(clusterProfiler)
library(org.Hs.eg.db)
library(pathview)
library(dplyr)

# Extract ENSEMBL IDs for up/down regulated genes
res_up_ensembl   <- rownames(expr.table[expr.table$diffexpressed == "UP", ])
res_down_ensembl <- rownames(expr.table[expr.table$diffexpressed == "DOWN", ])

# ------------------------
# GO ENRICHMENT (Biological Process)
# ------------------------

# Remove version numbers from ENSEMBL IDs
res_up_clean   <- sub("\\..*$", "", res_up_ensembl)
res_down_clean <- sub("\\..*$", "", res_down_ensembl)

# Convert ENSEMBL → SYMBOL
up_symbols   <- bitr(res_up_clean, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Hs.eg.db)$SYMBOL
down_symbols <- bitr(res_down_clean, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Hs.eg.db)$SYMBOL

# Build comparison list for clusterProfiler
comparison <- list(
  Upregulated   = na.omit(up_symbols),
  Downregulated = na.omit(down_symbols)
)

# GO enrichment (Biological Process)
CompareGO_BP <- compareCluster(
  comparison,
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  readable = TRUE
)

# Visualize GO enrichment
dotplot(CompareGO_BP, title = "GO Biological Process Enrichment")

# ------------------------
# KEGG ENRICHMENT
# ------------------------

# Convert ENSEMBL → ENTREZ IDs
up_entrez   <- bitr(res_up_clean, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
down_entrez <- bitr(res_down_clean, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID

# Combine all DEGs and remove NAs
all_entrez <- na.omit(c(up_entrez, down_entrez))

# KEGG enrichment
result_KEGG <- enrichKEGG(gene = all_entrez, organism = "hsa")

# Visualize KEGG enrichment
dotplot(result_KEGG, title = "KEGG Pathway Enrichment")

# Optional: view table
kegg_table <- as.data.frame(result_KEGG@result)
View(kegg_table)

# ------------------------
# PREPARE DATA FOR PATHVIEW
# ------------------------

# Combine up- and down-regulated DEGs with fold change
matrice <- rbind(
  expr.table[expr.table$diffexpressed == "UP", , drop = FALSE],
  expr.table[expr.table$diffexpressed == "DOWN", , drop = FALSE]
)

# Clean Ensembl IDs
matrice$ENSEMBL <- sub("\\..*$", "", rownames(matrice))

# Map Ensembl → ENTREZ
mapping <- bitr(matrice$ENSEMBL, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Merge mapping and keep fold change
matrice <- merge(matrice, mapping, by = "ENSEMBL", all.x = TRUE)
matrice <- matrice[, c("fc", "ENTREZID")]
matrice <- na.omit(matrice)

# Aggregate duplicates: take max absolute fold change
matrice <- aggregate(fc ~ ENTREZID, data = matrice, FUN = function(x) x[which.max(abs(x))])

# Set ENTREZ IDs as rownames
rownames(matrice) <- matrice$ENTREZID
matrice <- matrice[, "fc", drop = FALSE]

# Preview
head(matrice)

# ------------------------
# RUN PATHVIEW
# ------------------------
pathway_id <- "hsa04020"  # Calcium signaling pathway

pathview(
  gene.data  = matrice,
  species    = "hsa",
  pathway    = pathway_id,
  lowcol     = "blue",   # down-regulated
  highcol    = "red",    # up-regulated
  limit      = c(min(matrice[,1]), max(matrice[,1])),
  out.suffix = "Calcium_DEGs"
)

# ------------------------
# HUB ENRICHMENT
# ------------------------

# Function to clean ENSEMBL IDs and convert to SYMBOLs
convert_to_symbols <- function(genes) {
  genes_clean <- sub("\\..*$", "", names(genes))
  na.omit(bitr(genes_clean, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Hs.eg.db)$SYMBOL)
}

# Convert hub genes
hub_list <- list(
  Cancer_Hubs = convert_to_symbols(hubs.c),
  Normal_Hubs = convert_to_symbols(hubs.n),
  DiffNet_Hubs = convert_to_symbols(hubs.z)
)

# GO enrichment for hubs
hub_GO_BP <- compareCluster(
  hub_list,
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  readable = TRUE
)

# Visualize
dotplot(hub_GO_BP, title = "GO Biological Process Enrichment for Hubs")

# Convert SYMBOL → ENTREZ IDs for KEGG enrichment
hub_entrez <- lapply(hub_list, function(x) bitr(x, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID)

# KEGG enrichment for hubs
hub_KEGG <- compareCluster(
  hub_entrez,
  fun = "enrichKEGG",
  organism = "hsa"
)

# Visualize
dotplot(hub_KEGG, title = "KEGG Pathway Enrichment for Hubs")
