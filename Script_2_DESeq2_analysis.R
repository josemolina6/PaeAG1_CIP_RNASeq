# Transcriptomic analysis of response of Pseudomonas aeruginosa AG1 after exposure to Ciprofloxacin

#SCRIPT FOR DESEQ2 ANALYSIS
#----------------------------------------------
#Implemented by Jose Arturo Molina Mora
#University of Costa Rica
#----------------------------------------------

#Install R and DESeq2. Upon installing R, install DESeq2 on R:
# source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
#Import DESeq2 library in R


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install("DESeq2")
BiocManager::install("vsn")
BiocManager::install("hexbin")

library("vsn")
library("DESeq2")
library("gplots") 
library("readr")
library("reshape2")
library("ggplot2")
library("magrittr")
library("pheatmap")
library("genefilter")
library("grDevices")
library("RColorBrewer")
library("data.table")
library("hexbin")

#Load gene(/transcript) count matrix and labels
# Set the prefix for each output file name
outputPrefix <- "EDGEpro"
precountData <- as.matrix(read.csv("deseqFileAG1.csv", row.names="ID"))
precolData <- as.matrix(read.csv("Experimental_design.csv", row.names=1))

countData <- precountData[,1:9]
colData <- precolData[1:9,]

#Check all sample IDs in colData are also in CountData and match their orders
all(rownames(colData) %in% colnames(countData))
countData <- countData[, rownames(colData)]
all(rownames(colData) == colnames(countData))

dim(countData)
print(colData)

# Basic description of the data: number of reads per sample

barplot(colSums(countData)/1000000, main="Total reads per sample (million)")

# percentage of null counts per sample
prop.null <- apply(countData, 2, function(x) 100*mean(x==0))
print(head(prop.null))
barplot(prop.null, main="Percentage of null counts per sample", las=1)


#Create a DESeqDataSet from count matrix and labels. With "strain" se cargan los datos de todos las condiciones, aunque se puede seleccionar cuales cargar.
predds <- DESeqDataSetFromMatrix(countData = countData, 
                                colData = colData, design = ~ Sample)

# Filter data 
nrow(predds)
dds <- predds[ rowSums(counts(predds)) > 1, ]
nrow(dds)


#Run the default analysis for DESeq2
dds <- DESeq(dds)
head(dds)

print(sizeFactors(dds))

mean.counts <- rowMeans(counts(dds))
variance.counts <- apply(counts(dds), 1, var)
plot(x=mean.counts, y=variance.counts, pch=16, cex=0.3, main="Mean-variance relationship", log="xy")
abline(a=0, b=1)

# transform raw counts into normalized values, and visualize

par(mfrow=c(1,2))

samplecolors=c(rep("palevioletred1",3),rep("skyblue1",3),rep("yellowgreen",3),rep("orange",3))
boxplot(log2(counts(dds, normalized=FALSE)+1), main="Raw counts", col=samplecolors)
boxplot(log2(counts(dds, normalized=TRUE)+1), main="Normalized counts", col=samplecolors)

ntd <- normTransform(dds) # Este es el mismo que el anterior normalizado pero log2(n + 1)
rld <- rlogTransformation(dds, blind=T)
vsd <- varianceStabilizingTransformation(dds, blind=T)

par(mfrow=c(1,3))
boxplot(assay(ntd), main="ntd normalized counts (regular)", col=samplecolors)
boxplot(assay(rld), main="rld normalized counts", col=samplecolors)
boxplot(assay(vsd), main="vsd normalized counts", col=samplecolors)

par(mfrow=c(1,2))
boxplot(log2(counts(dds, normalized=FALSE)+1), main="Log(Raw counts)", col=samplecolors)
boxplot(assay(rld), main="Normalized counts (rlog)", col=samplecolors)

par(mfrow=c(3,1))

group <- c(1,1,1,2,2,2,3,3,3)
plot(0,xlim=c(-10,10),ylim=c(0,1),type="n",
     ylab="Density",xlab="log2(Counts+1)", main="Density ntd")
for(i in 1:ncol(assay(ntd))){
  lines(density(assay(ntd)[,i]),col=group[i])
}

plot(0,xlim=c(-10,10),ylim=c(0,1),type="n",
     ylab="Density",xlab="rlogt(Counts)", main="Density rld")
for(i in 1:ncol(assay(rld))){
  lines(density(assay(rld)[,i]),col=group[i])
}

plot(0,xlim=c(-10,10),ylim=c(0,1),type="n",
     ylab="Density",xlab="vst(Counts)", main="Density vsd")
for(i in 1:ncol(assay(vsd))){
  lines(density(assay(vsd)[,i]),col=group[i])
}

meanSdPlot(assay(ntd))

meanSdPlot(assay(rld))

meanSdPlot(assay(vsd))

par(mfrow=c(1,1))
plotDispEsts(dds)

# scatter plot of rlog transformations between Sample conditions
# nice way to compare control and experimental samples
# uncomment plots depending on size of array
# head(assay(rld))
plot(log2(1+counts(dds,normalized=T)[,c(1,4)]),col='black',pch=20,cex=0.3, main='Log2 transformed')
plot(assay(rld)[,c(2,5)],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")
plot(assay(rld)[,c(3,9)],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")

#Scatterplot
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.4/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
pairs(assay(rld), lower.panel = panel.smooth, upper.panel = panel.cor)

#Principal components plot shows clustering of samples

plotPCA(rld, intgroup=c("Experiment"))
plotPCA(rld, intgroup=c("Sample","Time","Antibiotic","Experiment"))
pcaData <- plotPCA(rld, intgroup=c("Sample","Time","Antibiotic","Experiment"), returnData=TRUE)

# Other option for PCA
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Experiment, shape=Time)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

# heatmap of data

# 1000 top expressed genes with heatmap.2
select <- order(rowMeans(counts(dds,normalized=T)),decreasing=T)[1:1000]
my_palette <- colorRampPalette(c("navy", "white", "firebrick3"))(n=1000)
heatmap.2(assay(rld)[select,], col=my_palette,
          scale="row", key=T, keysize=1, symkey=T,
          density.info="none", trace="none",
          cexCol=0.6, labRow=F,
          main=paste0(outputPrefix, "Clustering"))


# heatmap of data (type 2)

df <- data.frame(Strain = SummarizedExperiment::colData(dds)[,c("Sample","Antibiotic","Time")], row.names = rownames(SummarizedExperiment::colData(dds)))
select2 <- order(rowMeans(counts(dds,normalized=T)),decreasing=T)[1:100]
#Sin rownames:
pheatmap::pheatmap(assay(rld)[select2,], cluster_rows=TRUE, show_rownames=FALSE, cluster_cols=TRUE, annotation_col = df)
#Con rownames:
#pheatmap::pheatmap(assay(rld)[select2,], cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE, annotation_col = df)

df <- data.frame(SummarizedExperiment::colData(dds)[,c("Sample","Antibiotic","Time")], row.names = rownames(SummarizedExperiment::colData(dds)))
select2 <- order(rowMeans(counts(dds,normalized=T)),decreasing=T)[1:30]
colors<-colorRampPalette(rev(brewer.pal(n=7,name="RdYlGn")))(255)
pheatmap::pheatmap(assay(rld)[select2,], color = colors, cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE, annotation_col = df)

# Distance heatmap
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix( sampleDists ) 
rownames(sampleDistMatrix) <- paste(rld$Sample, vsd$Experiment, sep="-") 
colnames(sampleDistMatrix) <- paste(rld$Sample, vsd$Experiment, sep="-") 
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
hc <- hclust(sampleDists)
heatmap.2( sampleDistMatrix, Rowv=as.dendrogram(hc),symm=TRUE, trace="none", col=colors, margins=c(2,10), labCol=FALSE)


# Count outlier detection: Cook distance
par(mfrow=c(1,1))
res <- results(dds) # preliminar
W <- res$stat
maxCooks <- apply(assays(dds)[["cooks"]],1,max)
idx <- !is.na(W)
plot(rank(W[idx]), maxCooks[idx], xlab="rank of Wald statistic",
     ylab="maximum Cook's distance per gene",
     ylim=c(0,5), cex=.4, col=rgb(0,0,0,.3)) 
m <- ncol(dds)
p <- 3
abline(h=qf(.99, p, m - p))


# Clusters
df <- assay(rld)
df <- df - rowMeans(df)
df <- as.data.frame(df)
colnames(df)
res_km <- kmeans(df, 6, nstart = 10)
data_plot <- data.table(melt(data.table(class = as.factor(res_km$cluster), df)))
data_plot[, Time := rep(1:ncol(df), each = nrow(df))]
data_plot[, ID := rep(1:nrow(df), ncol(df))]
head(data_plot)

# prepare centroids
centers <- data.table(melt(res_km$centers))
setnames(centers, c("Var1", "Var2"), c("class", "Time"))
centers[, ID := class]
centers[, gr := as.numeric(as.factor(Time))]
head(centers)
head(data_plot)

# plot the results
ggplot(data_plot, aes(variable, value, group = ID)) +
  facet_wrap(~class, ncol = 2, scales = "free_y") +
  geom_line(color = "grey10", alpha = 0.65) +
  geom_line(data = centers, aes(gr, value),
            color = "firebrick1", alpha = 0.80, size = 1.2) +
  labs(x = "Samples", y = "Load (normalised)") +
  theme_bw()

# MDS

mds <- as.data.frame(colData(rld))  %>%
  cbind(cmdscale(sampleDistMatrix))
ggplot(mds, aes(x = `1`, y = `2`, color = Sample, shape = Time)) +
  geom_point(size = 3) + coord_fixed()

# SAVE GENERAL DATA 

# save normalized values (raw, ntd and rlog) 
write.csv(as.data.frame(counts(dds, normalized=FALSE)),file = paste0(outputPrefix, "-original-counts.txt"))
write.csv(as.data.frame(counts(dds, normalized=TRUE)),file = paste0(outputPrefix, "-ntd-normalized-counts.txt"))
#Allgenes (no min-filter)
predds <- DESeq(predds)
write.csv(as.data.frame(counts(predds, normalized=TRUE)),file = paste0(outputPrefix, "-ntd-normalized-counts-ALL-genes.txt"))

write.csv(as.data.frame(assay(rld)),file = paste0(outputPrefix, "-rlog-transformed-counts.txt"))
write.csv(as.data.frame(assay(vsd)),file = paste0(outputPrefix, "-vsd-transformed-counts.txt"))


# Differential expression analysis

resultsNames(dds)

Comparacion<-matrix(c("Comparison","BvrsA","CvrsB","CvrsA","Number of DEG","-","-","-"),  nrow = 2,ncol = 4, byrow = TRUE)
Comparacion

# -----------------------------------

#CASE 1: 2.5h (B) vrs 0 h (A) 
textconstrast<- "B-vrs-A"
contrast <- c("Sample", "B","A")

res <- results(dds, contrast)  
head(res)

#Sort by adjusted p-value and display
resorder <- res[order(res$padj),]
resdata <- merge(as.data.frame(resorder), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata)[1] <- 'gene'
write.csv(resdata, file = paste0(outputPrefix,"_",textconstrast,"_all_data.csv"))

hist(res$pvalue, col="lightblue", main=paste0(textconstrast, ": Histogram of raw P-values"), breaks=20, xlab="P-value")


# MA plot of RNAseq data for entire dataset
# http://en.wikipedia.org/wiki/MA_plot
# genes with padj < 0.1 are colored Red

plotMA(res, ylim=c(-8,8),main = paste0(textconstrast,": MAplot"))

## Draw a volcano plot
plot(x=res$log2FoldChange, y=-log10(res$padj), 
     xlab="log2(Fold-Change)", ylab="-log10(adjusted P-value)",
     col=ifelse(res$padj<=0.05, "red", "black"), main=paste0(textconstrast,": Volcano plot"))
grid() ## Add a grid to the plot
abline(v=0) ## plot the X axis

#DEG number
Comparacion[2,2]<- sum(res$padj <= 0.05, na.rm=TRUE)
Comparacion

# filter results by p value
# order results by padj value (most significant to least)

ressub= subset(res, padj<0.05)
ressub <- ressub[order(ressub$padj),]

# should see DataFrame of baseMean, log2Foldchange, stat, pval, padj

# save data results and normalized reads to csv
resdatasub <- merge(as.data.frame(ressub), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdatasub)[1] <- 'gene'
write.csv(resdatasub, file = paste0(outputPrefix,"_",textconstrast,"_only_DEG_data.csv"))

#Compare
par(mfrow=c(1,3))
plotCounts(dds, gene=rownames(ressub)[1], intgroup="Sample", normalized=TRUE)
plotCounts(dds, gene=rownames(ressub)[2], intgroup="Sample", normalized=TRUE)
plotCounts(dds, gene=rownames(ressub)[3], intgroup="Sample", normalized=TRUE)

# Save gene list
genes_BA<-resdatasub$gene

# -----------------------------------

#CASE 2: 5 h (C) vrs 0 h (A) 
textconstrast<- "C-vrs-A"
contrast <- c("Sample", "C","A")

res <- results(dds, contrast)
head(res)

#Sort by adjusted p-value and display
resorder <- res[order(res$padj),]
resdata <- merge(as.data.frame(resorder), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata)[1] <- 'gene'
write.csv(resdata, file = paste0(outputPrefix,"_",textconstrast,"_all_data.csv"))

par(mfrow=c(1,1))
hist(res$pvalue, col="lightblue", main=paste0(textconstrast, ": Histogram of raw P-values"), breaks=20, xlab="P-value")

# MA plot of RNAseq data for entire dataset
# http://en.wikipedia.org/wiki/MA_plot
# genes with padj < 0.1 are colored Red

plotMA(res, ylim=c(-8,8),main = paste0(textconstrast,": MAplot"))

## Draw a volcano plot
plot(x=res$log2FoldChange, y=-log10(res$padj), 
     xlab="log2(Fold-Change)", ylab="-log10(adjusted P-value)",
     col=ifelse(res$padj<=0.05, "red", "black"), main=paste0(textconstrast,": Volcano plot"))
grid() ## Add a grid to the plot
abline(v=0) ## plot the X axis

# DEG number
Comparacion[2,4]<- sum(res$padj <= 0.05, na.rm=TRUE)
Comparacion

# filter results by p value
# order results by padj value (most significant to least)

ressub= subset(res, padj<0.05)
ressub <- ressub[order(ressub$padj),]

# should see DataFrame of baseMean, log2Foldchange, stat, pval, padj

# save data results and normalized reads to csv
resdatasub <- merge(as.data.frame(ressub), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdatasub)[1] <- 'gene'
write.csv(resdatasub, file = paste0(outputPrefix,"_",textconstrast,"_only_DEG_data.csv"))

#Compare
par(mfrow=c(1,3))
plotCounts(dds, gene=rownames(ressub)[1], intgroup="Sample", normalized=TRUE)
plotCounts(dds, gene=rownames(ressub)[2], intgroup="Sample", normalized=TRUE)
plotCounts(dds, gene=rownames(ressub)[3], intgroup="Sample", normalized=TRUE)

# Save gene list
genes_CA<-resdatasub$gene



# -----------------------------------


