# Import expresion data (normalized/log) -------------
expr_data_bot <- read.csv("/Users/pamelacamejo/Documents/IBIO/Elena_Vidal/Projects/Botritis-Trichoderma/Data/redatosrnaseqconsu/NormData_botrytis.log.txt", sep="\t")
expr_data_bot_t <- data.frame(t(expr_data_bot[,-1]))

# Remove columns with zero variance
expr_data_bot_filt <- expr_data_bot_t[,apply(expr_data_bot_t, 2, var, na.rm=TRUE) != 0]

#Trait data -------------
traitData = read.csv("/Users/pamelacamejo/Documents/IBIO/Elena_Vidal/Projects/Botritis-Trichoderma/Data/Botrytis_treatment.csv")
Samples = rownames(expr_data_bot_filt);
traitRows = match(Samples, traitData$Sample)
datTraits = traitData[traitRows, -1]
rownames(datTraits) = traitData[traitRows, 1]

# PCA plot, option 1: ggbiplot -------------

install.packages("devtools")
library(devtools)
install_github("vqv/ggbiplot")
library(ggbiplot)

expr_data_bot_pca <- prcomp(expr_data_bot_filt,center = TRUE,scale. = TRUE)
ggbiplot(expr_data_bot_pca, var.axes = FALSE,labels=rownames(expr_data_bot_filt), groups=datTraits$G, ellipse=TRUE)

# PCA plot, option2: ggplot2 -------------

library(ggplot2)

dat.pca <- data.frame(expr_data_bot_pca$x) 
percentVar <- data.frame(summary(expr_data_bot_pca)$importance)

ggplot(dat.pca, aes(PC1,PC2)) + geom_point(aes(shape=datTraits$T, colour=datTraits$G),size=4) + 
  stat_ellipse(aes(fill=factor(datTraits$G)), type = "norm", geom="polygon", level=0.9, ,alpha=0.2, show.legend=F) +
  scale_shape_manual(values=c(15,16,17,18)) +
  scale_colour_manual(values = c("#00AFBB", "#E7B800")) +
  scale_fill_manual(values = c("#00AFBB", "#E7B800"))  + 
  theme(
    axis.line = element_line(colour = "black", size = 0.5, linetype = "solid"),
    panel.background = element_rect(fill = "white", colour = "black", size = 1, linetype = "solid"),
    panel.grid.major = element_line(colour ="grey", size = , linetype = "dashed"),
    panel.grid.minor = element_line(colour ="black", size = , linetype = "dashed")
  ) +
  xlab(paste0("PC1: (",round(percentVar[2,1] *100,2),"%)")) +
  ylab(paste0("PC2: (",round(percentVar[2,2] * 100,2),"%)")) + 
  coord_cartesian() +
  labs(shape = "treatment", colour = "genotype")

# PCA plot, option3: FactoMineR & factoextra -------------

install.packages(c("FactoMineR", "factoextra"))
library("FactoMineR")
library("factoextra")

res.pca <- PCA(expr_data_bot_filt, scale.unit = TRUE, ncp = 5, graph = TRUE)

# Useful plot to check how much variance is explained by each PC
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))

# PCA plot
fviz_pca_ind(res.pca, col.ind = datTraits$G, pointsize = 2,
             palette = c("#00AFBB", "#E7B800"), addEllipses = TRUE, legend.title = "Groups") +
  theme_minimal() +
  scale_shape_manual(values=c(19,20,21))


