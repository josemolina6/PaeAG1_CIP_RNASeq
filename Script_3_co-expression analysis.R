# Transcriptomic analysis of response of Pseudomonas aeruginosa AG1 after exposure to Ciprofloxacin

#SCRIPT FOR CO-EXPRESSION ANALYSIS
#----------------------------------------------
#Implemented by Jose Arturo Molina Mora
#University of Costa Rica
#----------------------------------------------

#----------------------------------------------
#
#  PART 1: MODULES AND CO-EXPRESSION NETWORKS
#
#----------------------------------------------

# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA. This helps speed up certain calculations.
# At present this call is necessary for the code to work.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments. 
# See note above.
enableWGCNAThreads()
# Load the data saved in the first part

datExpr <- read.csv("deseqFileAG1_518_DEGs.csv", row.names = 1)
datTraits <- read.csv("PaeAG1_traits.csv", row.names = 1)

#----------------------------------------------


# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#----------------------------------------------

net = blockwiseModules(log2(datExpr+1), power =9,
                       deepSplit = 2, 
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.15,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "AG1TOM", 
                       verbose = 3)
dev.off()

table(net$colors)

#----------------------------------------------

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

me = moduleEigengenes(log2(datExpr+1),colors=net$colors)$eigengenes
dev.off()


par(mfrow = c(2,3));
for (i in 1:5)
        plot(me[,i])


par(mfrow = c(1,1));
D <- 1-(1+cor(me,use="p"))/2
hcm = hclust(as.dist(D), method="average")
plot(hcm)

#----------------------------------------------

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree, 
     file = "networkConstruction-auto.RData")

TOM = TOMsimilarityFromExpr(datExpr, power = 9);

inmod = which(net$colors==6)

modgenes = rownames(data)

colnames(TOM) <- modgenes; rownames(TOM) <- modgenes;

exportNetworkToCytoscape(TOM[inmod,inmod],
                         edgeFile='wgcna_edges.txt',
                         nodeFile='wgcna_nodes.txt',
                         weighted=FALSE, threshold = 0.985,
                         nodeNames=modgenes)

#----------------------------------------------
#
#  PART 2: TRAITS VRS MODULES
#
#----------------------------------------------

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

#----------------------------------------------

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

moduleTraitCor

#----------------------------------------------
#
#  PART 3: ANNOTATIONS
#
#----------------------------------------------

#
# Define variable weight containing the weight column of datTrait
Time = as.data.frame(datTraits$Time);
names(Time) = "Time"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, Time, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(Time), sep="");
names(GSPvalue) = paste("p.GS.", names(Time), sep="");

annot = read.csv(file = "PaeAG1_final_annot.csv");
dim(annot)
names(annot)
probes = names(datExpr)
probes2annot = match(probes, annot$LocusAG1)
# The following is the number or probes without annotation:
sum(is.na(probes2annot))
# Should return 0.

#----------------------------------------------


# Create the starting data frame
geneInfo0 = data.frame(LocusAG1 = probes,
                       geneSymbol = annot$gene_symbol[probes2annot],
                       LocusLinkID = annot$LocusLinkID[probes2annot],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)

# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, Time, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
        oldNames = names(geneInfo0)
        geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                               MMPvalue[, modOrder[mod]]);
        names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                             paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Time));
geneInfo = geneInfo0[geneOrder, ]


#----------------------------------------------

Write.csv(geneInfo, file = "Results_COexprss_AllgeneInfo.csv")

