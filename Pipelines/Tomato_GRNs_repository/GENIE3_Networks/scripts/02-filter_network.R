```r

# Load required libraries
library(dplyr)
library(data.table)

# Set working directory
setwd(".")

# Read raw pGRN output
net <- fread("results/raw_pGRN_tissue1.txt", header = TRUE, sep = " ")

# Rename columns if needed
colnames(net) <- c("row_id", "TF", "TARGET", "weight")

# Sort by TARGET and weight descending
sorted_net <- net %>%
  arrange(desc(TARGET), desc(weight)) %>%
  group_by(TARGET) %>%
  mutate(percent = round(row_number() * 100 / n(), 2)) %>%
  filter(percent <= 2)

# Save filtered network
write.table(sorted_net, "results/filtered-pGRN_tissue1.txt", sep = "\t", row.names = FALSE, quote = FALSE)

```
