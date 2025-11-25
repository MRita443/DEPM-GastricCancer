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
BiocManager::install("DGCA")
library(DGCA) # package for differential co-expression network

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
  geom_hline(yintercept=-log10(0.01), col="red")+
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
#what if it's even lower?
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
set.edge.attribute(net.n, "edgecolor", ifelse(net.n %e% "weights" > 0, "red", "blue"))
# Visualizing
ggnet2(net.n, color = "color", alpha = 0.7, size = 2,
       edge.color = "edgecolor", edge.alpha = 1, edge.size = 0.15)+
  guides(size = "none") 

# See which genes are hubs both in cancer and normal tissue
intersect(names(hubs.c), names(hubs.n))
##### CHECK if needed define adj.mat.n with good p-value threshold
adj.mat.n1 <- rho.n * (qval.n <= 1e-3)
adj.mat.n2 <- rho.n * (qval.n <= 1e-4) #too much?
##### 


# NOT ASKED IN THE PROJECT, BUT USEFUL FOR VISUALIZATION
####
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
set.edge.attribute(net.hub, "ecolor", ifelse(net.hub %e% "weights" > 0, "red", "blue"))

ggnet2(net.hub,  color = "color",alpha = 0.9, size = 2, 
       edge.color = "ecolor", edge.alpha = 0.9,  edge.size = 0.15, 
       node.label = names(hubs.c), label.color = "black", label.size = 3)+
  guides(size = "none") 
####
####

# Question 4 
