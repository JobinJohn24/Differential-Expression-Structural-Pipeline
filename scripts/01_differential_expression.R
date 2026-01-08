suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
})

# 2. Load Data 
# This allows the script to run on any machine without hard coded paths
print("Select your COUNTS file (TCGA-LUAD.star_counts.tsv)...")
counts_file <- file.choose()

print("Loading data... (This may take a moment)")
raw_data <- read.table(counts_file, header = TRUE, row.names = 1, sep = "\t")

# 3. Pre-Process: Un-log and Round
# TCGA data is log2(x+1), DESeq2 requires integers
print("Processing raw counts...")
count_matrix <- round(2^raw_data - 1)

# 4. Extract Metadata from Column Names
# TCGA format: TCGA-38-7271-01A (01-09 = Tumor, 10+ = Normal)
sample_ids <- colnames(count_matrix)
sample_codes <- substr(sample_ids, 14, 15)
conditions <- ifelse(as.numeric(sample_codes) < 10, "Tumor", "Normal")

metadata <- data.frame(
  row.names = sample_ids,
  condition = factor(conditions, levels = c("Normal", "Tumor"))
)

print(paste("Detected:", sum(metadata$condition == "Tumor"), "Tumors and", 
            sum(metadata$condition == "Normal"), "Normals"))

# 5. Run DESeq2 Analysis
print("Running DESeq2 statistical model... (Grab a coffee)")
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = metadata,
                              design = ~ condition)

# Set "Normal" as the baseline reference
dds$condition <- relevel(dds$condition, ref = "Normal")
dds <- DESeq(dds)

# 6. Filter results
print("Filtering for Significant Genes...")
res <- results(dds, contrast = c("condition", "Tumor", "Normal"))
res_df <- as.data.frame(res)

# Criteria: Significant (padj < 0.05) AND High Change (Log2FC > 1.5)
# Sorted by P-value (Most reliable first)
top_candidates <- res_df %>%
  filter(padj < 0.05 & log2FoldChange > 1.5) %>%
  arrange(padj) %>%
  head(2000) # Increased limit to ensure structural targets

# Add Gene IDs for the next pipeline stage
top_candidates$Ensembl_ID <- row.names(top_candidates)

output_filename <- "results/phase1_targets.csv"
write.csv(top_candidates, output_filename, row.names = FALSE)

print(paste("✅ SUCCESS! Top candidates saved to:", output_filename))
