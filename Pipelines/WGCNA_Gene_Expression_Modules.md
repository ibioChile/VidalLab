# Gene Expression Modules using WGCNA

This pipeline explains how to create gene expression modules (clusters) using the WGCNA R package. The pipeline's input correspond to a list of differentially expressed (DE) genes (it could also be all genes) along with thir gene expression (log scale). In this case, we wil use a set of DE genes from *Botrytis cinerea* when in contact with *Trichoderma*. Input data can be found [here](https://github.com/ibioChile/VidalLab/tree/master/Data/Botrytis_Trichoderma_RNAseq)

1. Install the following libraries

```
install.packages("WGCNA")
install.packages("readxl")
library(WGCNA)
allowWGCNAThreads()
library(readxl)
library(pheatmap)
options(stringsAsFactors = FALSE)
```

2. Import gene expression data

```
expr_data_bot <- read.csv("NormData_botrytis.log.txt", sep="\t")
```

3. Import list of DE genes (T: Trichoderma treatment, G: genotype (*Botrytis* wild type or mutant)). 

```
T_DE <- read_xlsx("botrytis.results.01.a2v.coef.xlsx",sheet = "T",col_names = "DE")
G_DE <- read_xlsx("botrytis.results.01.a2v.coef.xlsx",sheet = "G",col_names = "DE")
TG_DE <- read_xlsx("botrytis.results.01.a2v.coef.xlsx",sheet = "TG",col_names = "DE")
DE_list <- rbind(T_DE,G_DE,TG_DE) 
```

```
#remove duplicated DE genes
DE_data <- expr_data_bot[expr_data_bot$X %in% unique(DE_list$DE),]

#row.names = gene Ids
datExpr0 = as.data.frame(t(DE_data[, -1]))
names(datExpr0) = DE_data$X
```

4. Import trait data

```
traitData = read.csv("Botrytis_treatment.csv")
Samples = rownames(datExpr0);
traitRows = match(Samples, traitData$Sample)
datTraits = traitData[traitRows, -1]
rownames(datTraits) = traitData[traitRows, 1]
```

5. Check if samples cluster according to experiment conditions (not needed for modules estimation, but recommended).

```
sampleTree = hclust(dist(datExpr0), method = "average")
# Convert traits to a color representation.
traitColors = labels2colors(datTraits,colorSeq=c("brown1","darkblue","darkgoldenrod1","forestgreen"))
# Plot the sample dendrogram and the colors underneath.
par(cex = 0.6);
par(mar = c(0,4,2,0))
plotDendroAndColors(sampleTree, traitColors,
                    groupLabels = names(datTraits), 
                    rowText = datTraits,
                    rowTextAlignment = c("left"),
                    main = "Sample dendrogram and trait heatmap")
dev.off()
```

This is how the tree should look after clusterization:

![SampleClustering_DE](https://user-images.githubusercontent.com/53570955/65054502-07994700-d944-11e9-8dad-8e1e4802f6db.jpg)

6. Select parameters for module clusterization.

```
# Choose a set of soft-thresholding powers. Iterate this number until finding a plato in the curve.
powers = c(seq(from = 1, to=70, by=4))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5, networkType ="signed")
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
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
```

![soft-thresholding powers_parameters](https://user-images.githubusercontent.com/53570955/65054812-91e1ab00-d944-11e9-8b6a-df0e036f6ca8.jpg)

In this case, we picked a soft thresholding value of 65 because it allows an r2 close to 0.8 (it is a local peak for the r2) and the mean connectivity is still above 0.


7. Modules generation. 

```
# Calculate adjacencies
softPower = 65
adjacency = adjacency(datExpr0, power = softPower,type = "signed")
# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");

# Set the minimum module size:
minModuleSize = 20;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

# Calculate eigengenes
MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
```

8. Module merging

```
#We choose a height cut of 0.20 to merge modules (modules with a correlation > 0.8 are merged).
MEDissThres = 0.20
# Plot the cut line into the dendrogram
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
abline(h=MEDissThres, col = "red")
dev.off()
```
The histogram generated shows how different modules are clustered and draws a line for merging.

![clustering_modules-DE](https://user-images.githubusercontent.com/53570955/65058585-dff9ad00-d94a-11e9-878d-9021a726b42a.jpg)

```
# Call an automatic merging function
merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
```

This gene dendogram shows how modules are clustered based on distance:

![geneDendro-DE](https://user-images.githubusercontent.com/53570955/65058552-d4a68180-d94a-11e9-91ae-d7c622bae074.jpg)

```
# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
```

9. Plot modules gene expression, represented by the eigengenes, at the different conditions. 

```
# Heatmap of eigengenes per module
merged_MEs <- merge(MEs,datTraits, by=0)
merged_MEs_sorted <- merged_MEs[order(merged_MEs$G,merged_MEs$T),] 
merged_MEs_sorted <- subset(merged_MEs_sorted, select=-c(MEgrey))
  
pheatmap_colour = list(
  G = c("bc∆/∆"="darkred","bcWT"="darkblue"),
  T = c("bc∆/∆"="darkred","ta∆/∆"="darkgoldenrod1","taWT" ="forestgreen","bcWT"="darkblue"),
  module=c("MEorange" = "orange","MEbrown" = "brown","MEpurple" = "purple"))

modules_colour = data.frame(colnames(merged_MEs_sorted[2:4]),row.names = colnames(merged_MEs_sorted[2:4]))
colnames(modules_colour) <- "module"

#Change the columns to plot according to results.
pheatmap(t(merged_MEs_sorted[,2:4]), cluster_rows = FALSE, main = "Eigengene expression of each module",
         cluster_cols= FALSE,show_colnames = F, 
         annotation_col = merged_MEs_sorted[,5:6], 
         annotation_colors = pheatmap_colour,
         #annotation_row = modules_colour,
         cellheight = 30,cellwidth = 30)
```

![modules_expression_DE](https://user-images.githubusercontent.com/53570955/65058536-c9535600-d94a-11e9-8864-5b40956b7812.jpg)


