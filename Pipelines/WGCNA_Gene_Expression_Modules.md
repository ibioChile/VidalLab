# Gene Expression Modules using WGCNA

This pipeline explains how to create gene expression modules (clusters) using the WGCNA R package. The pipeline's input correspond to a list of differentially expressed (DE) genes (it could also be all genes) along with thir gene expression (log scale). In this case, we wil use a set of DE genes from *Botrytis cinerea* when in contact with *Trichoderma*. 

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

#remove duplicated DE genes

```
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




