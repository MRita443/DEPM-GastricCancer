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

#setwd('/home/yunbao/Bureau/Sapienza/DEPML/DEPM-GastricCancer')

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
deg.genes <- rownames(expr.table[abs(expr.table$fc) >= 1.2 & expr.table$pval.fc.fdr <=0.01,]) 
head(expr.table[deg.genes,], 10)
write.table(expr.table[deg.genes,], file = "DEG.csv", sep = ";")

#volcano plot
expr.table$diffexpressed <- "/";
expr.table$diffexpressed[expr.table$fc >= 1.2 & expr.table$pval.fc.fdr <= 0.01] <- "UP"
expr.table$diffexpressed[expr.table$fc <= -1.2 & expr.table$pval.fc.fdr <= 0.01] <- "DOWN"
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
adj.mat.c <- rho.c
# keep only |correlation| > 0.7
adj.mat.c  <- adj.mat.c * (abs(rho.c) >= 0.7)
adj.mat.c <- adj.mat.c * (qval.c <= 1e-3)
adj.bin.c  <- adj.mat.c
adj.bin.c[abs(adj.mat.c)>0] <- 1 # get a binary version of the same matrix

# NORMAL NETWORK
cor.mat.n <- corr.test(t(filtr.expr.N.DEGs), use="pairwise", method="spearman", adjust="fdr", ci=FALSE)
rho.n <- cor.mat.n$r
diag(rho.n) <- 0
qval.n <- cor.mat.n$p
qval.n[lower.tri(qval.n)] <- t(qval.n)[lower.tri(qval.n)]
adj.mat.n <- rho.n
adj.mat.n  <- adj.mat.n * (abs(rho.n) >= 0.85) # we use a more selective threshold to have comparable networks (between cancer and normal)
adj.mat.n <- adj.mat.n * (qval.n <= 1e-5)
adj.bin.n  <- adj.mat.n
adj.bin.n[abs(adj.mat.n)>0] <- 1 # get a binary version of the same matrix

#7 : Co-expression networks
# ANALYSIS PART
# Aims to check if network is scale-free or not

# CANCER NETWORK :
#build the network
net.c <- network(adj.mat.c, matrix.type="adjacency",ignore.eval = FALSE, names.eval = "weights", directed = F)
### ANALYSIS Q3.1
network.density(net.c) #how much is the network connected
network.size(net.c) #n nodes
network.edgecount(net.c)  #n of edges
clustcoeff(adj.mat.c, weighted = FALSE)$CC # clustcoeff measures how connected are my neighbors between them
sum(adj.mat.c != 0) / 2
#how many positive/negative correlations? 
sum(adj.mat.c > 0) / 2
sum(adj.mat.c < 0) / 2
degree.c <- rowSums(adj.mat.c != 0)
names(degree.c) <- rownames(adj.mat.c)
degree.c <- sort(degree.c, decreasing = T)
head(degree.c,10)
sum(degree.c == 0) # check how many unconnected nodes
# check if scale-free
hist(degree.c) # This is obviously a scale-free network !
###
### ANALYSIS Q3.2
# Now we want find the hubs (5% of nodes with highest degree values)
x <- quantile(degree.c[degree.c>0],0.95) #top 5% of nodes
x # x=50 means the 95% highest degree node is >= 50 links
abline(v=x, col="red")
hubs.c <- degree.c[degree.c>=x]
names(hubs.c) 
#let's enrich them
write.table(hubs.c, file = "hubs_c.csv", sep = ";")
###

# Anotate the hubs
net.c %v% "type" = ifelse(network.vertex.names(net.c) %in% names(hubs.c),"hub", "non-hub")
net.c %v% "color" = ifelse(net.c %v% "type" == "hub", "tomato", "deepskyblue3")
network::set.edge.attribute(net.c, "edgecolor", ifelse(net.c %e% "weights" > 0, "red", "blue"))
# Visualizing
ggnet2(net.c, color = "color", alpha = 0.7, size = 2,  #mode= c("x","y"),
       edge.color = "edgecolor",
       edge.alpha = 1, edge.size = 0.15)+
      guides(size = "none")

# NORMAL NETWORK :
net.n <- network(adj.mat.n, matrix.type="adjacency",ignore.eval = FALSE, names.eval = "weights")
### ANALYSIS Q3.1
network.density(net.n)
network.size(net.n)
network.edgecount(net.n)
clustcoeff(adj.mat.n, weighted = FALSE)$CC
sum(adj.mat.n != 0) /2
sum(adj.mat.n > 0) /2
sum(adj.mat.n < 0) /2
degree.n <- rowSums(adj.mat.n != 0)
names(degree.n) <- rownames(adj.mat.n)
degree.n <- sort(degree.n, decreasing = T)
head(degree.n,10)
sum(degree.n == 0) #unconnected nodes 
hist(degree.n) # This is obviously a scale-free network !
###
### ANALYSIS Q3.2
# Now we want find the hubs (5% of nodes with highest degree values)
y <- quantile(degree.n[degree.n>0],0.95) 
y
hist(degree.n)
abline(v=y, col="red")
hubs.n <- degree.n[degree.n>=y]
names(hubs.n)
#let's enrich them
write.table(hubs.n, file = "hubs_n.csv", sep = ";")
###

# Annotate this
net.n %v% "type" = ifelse(network.vertex.names(net.n) %in% names(hubs.n),"hub", "non-hub")
net.n %v% "color" = ifelse(net.n %v% "type" == "hub", "tomato", "deepskyblue3")
set.edge.attribute(net.n, "edgecolor", ifelse(net.n %e% "weights" > 0, "red", "blue"))
# Visualizing
ggnet2(net.n, color = "color", alpha = 0.7, size = 2,
       edge.color = "edgecolor", edge.alpha = 1, edge.size = 0.15)+
  guides(size = "none") 

### ANALYSIS Q3.2
# See which genes are hubs both in cancer and normal tissue
length(names(hubs.c)) # 59
length(names(hubs.n)) # 56
length(intersect(names(hubs.c), names(hubs.n))) # 18 -> they share ~1/3 of their hubs

hubs_intersect_coexpr <- intersect(names(hubs.c), names(hubs.n))
hubs_intersect_coexpr
gained_hubs_coexpr <- setdiff(names(hubs.c),names(hubs.n))
gained_hubs_coexpr
lost_hubs_coexpr <- setdiff(names(hubs.n),names(hubs.c))
lost_hubs_coexpr
# TODO Report : with KEGG look to fonctions which desapeared (hubs in normal not hub in cancer anymore). Also to the fonctions which appeared (hubs in cancer but not in normal). 
###

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
set.edge.attribute(net.hub, "ecolor", ifelse(net.hub %e% "weights" > 0, "red", "blue"))

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
### ANALYSIS Q4.1
network.density(net.z) #how much is the network connected
network.size(net.z) #n nodes
network.edgecount(net.z)  #n of edges
clustcoeff(Z.bin, weighted = FALSE)$CC # clustcoeff measures how connected are my neighbors between them
degree.z <- rowSums(Z.bin != 0)
names(degree.z) <- rownames(Z.bin)
degree.z <- sort(degree.z, decreasing = T)
head(degree.z,10)
sum(degree.z == 0) # check how many unconnected nodes
# Plot the histogram to see if we have a scale free network or not
hist(degree.z) # This is obviously a scale-free network !
###
### ANAYSIS Q4.2
# Now we want find the hubs (5% of nodes with highest degree values)
x.z <- quantile(degree.z[degree.z>0],0.95) #top 5% of nodes
x.z
hist(degree.z)
abline(v=x.z, col="red")
hubs.z <- degree.z[degree.z>=x.z]
names(hubs.z) 
#let's enrich them
write.table(hubs.c, file = "hubs_z.csv", sep = ";")
###

# Anotate the hubs
net.z %v% "type" = ifelse(network.vertex.names(net.z) %in% names(hubs.z),"hub", "non-hub")
net.z %v% "color" = ifelse(net.z %v% "type" == "hub", "tomato", "deepskyblue3")
network::set.edge.attribute(net.z, "edgecolor", ifelse(net.z %e% "weights" > 0, "red", "blue"))
# Visualizing
ggnet2(net.z, color = "color", alpha = 0.7, size = 2,  #mode= c("x","y"),
       edge.color = "edgecolor", edge.alpha = 1, edge.size = 0.15)+
  guides(size = "none")

### ANALYSIS Q4.2
length(names(hubs.c)) # 59
length(names(hubs.n)) # 56
length(names(hubs.z)) # 70
# TODO Report : with KEGG look to functions in hubs of names(hubs.z) because it shows which biological process are really perturbated

length(intersect(names(hubs.c), names(hubs.n))) # 18 -> ~half of co_expr cancer genes hubs are also hubs in co_expr normal tissue !
length(intersect(names(hubs.z), names(hubs.n))) # 49 co_expr normal hubs are also in diff_coexpr hubs !
length(intersect(names(hubs.c), names(hubs.z))) # 21 co_expr cancer hubs are also in diff_coexpr hubs !
length(intersect(names(hubs.z), hubs_intersect_coexpr)) # 14 genes are always hubs in coexpr tissues and diffcoexpr newtworks

hubs_intersect_diffcoexpr_and_coexpr <- intersect(names(hubs.z), hubs_intersect_coexpr)
hubs_intersect_diffcoexpr_and_coexpr
# FNBP1 TMTC1 PNCK, CSRP1 and others were all hubs in cancer and normal, but still after diffcoexpr

gained_hubs_diffcoexpr <- setdiff(names(hubs.z), names(hubs.c))
gained_hubs_diffcoexpr <- setdiff(gained_hubs_diffcoexpr, names(hubs.n))
gained_hubs_diffcoexpr
# TODO Report : with KEGG look to fonctions which are in diff_coexpr but non hubs neither in cancer or normal (cf gained_hubs_diffcoexpr) those are the relations which changed the more !
###


# Question 5 : Patient Similarity Network (PSN)