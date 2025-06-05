```r

# Load required libraries
library(dplyr)
library(GENIE3)

# Set working directory
setwd(".")

# Load expression matrix
expr_df <- read.table("data/tissue1_counts.txt", header = TRUE, sep = "\t")
rownames(expr_df) <- expr_df$Gene
expr_df$Gene <- NULL

# Convert to matrix
expr_matrix <- data.matrix(expr_df)

# Load TF list
tf_list <- read.table("data/CuratedlistofTFs.txt", header = FALSE)
regulators <- tf_list$V1

# Run GENIE3
set.seed(123)
weight_matrix <- GENIE3(expr_matrix, regulators = regulators, nTrees = 2000, nCores = 40, verbose = TRUE)

# Get regulatory links
link_list <- getLinkList(weight_matrix, reportMax = NULL)

# Create results directory if not present
if (!dir.exists("results")) {
  dir.create("results")
}

# Save raw pGRN
write.table(link_list, "results/raw_pGRN_tissue1.txt", sep = " ", quote = FALSE, row.names = FALSE)

```
