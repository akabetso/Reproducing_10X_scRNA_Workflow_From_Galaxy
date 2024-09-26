# Load the necessary libraries
library(DropletUtils)
library(Matrix)
library(SingleCellExperiment)
library(ggplot2)


# Directories
results_dir <- "../results/DropletUtils"
matrix_path <- "STARsolo_STARsolo/results/STARsoloGene/raw/matrix.mtx"
barcodes_path <- "STARsolo_STARsolo/results/STARsoloGene/raw/barcodes.tsv"
genes_path <- "STARsolo_STARsolo/results/STARsoloGene/raw/features.tsv"
counts_dir <- "STARsolo_STARsolo/results/STARsoloGene/raw"

# Read in the matrix
raw_counts <- readMM(file = matrix_path)
dim(raw_counts)
# Read the barcodes and gene information
barcodes <- read.table(barcodes_path, header = FALSE, stringsAsFactors = FALSE)
genes <- read.table(genes_path, header = FALSE, stringsAsFactors = FALSE)
head(genes)
# Assign barcodes to the columns of the matrix and genes to the rows
colnames(raw_counts) <- barcodes$V1
rownames(raw_counts) <- genes$V1


sce <- SingleCellExperiment(
  assays = list(counts = raw_counts),
  rowData = DataFrame(gene_id = genes$V1, gene_symbol = genes$V2),
  colData = DataFrame(cell_barcodes = barcodes$V1)
)

head(sce)
# Apply the emptyDrops method
e.out <- emptyDrops(counts(sce), lower = 100)

# Check if the output has FDR values and filter cells
if ("FDR" %in% colnames(e.out)) {
    filtered_barcodes <- which(e.out$FDR <= 0.01)
} else {
    stop("FDR values not found in the output of emptyDrops.")
}
head(filtered_barcodes)
# Subset the SCE object if any barcodes passed
if (length(filtered_barcodes) > 0) {
    sce_filtered <- sce[, filtered_barcodes]
} else {
    stop("No barcodes passed the filtering criteria.")
}
head(sce_filtered)
# Check the dimensions of the filtered matrix
dim(sce_filtered)

# Save the filtered matrix in 10X format
write10xCounts(
  path = results_dir, 
  x = counts(sce_filtered), 
  gene.id = rownames(sce_filtered),
  gene.symbol = rowData(sce_filtered)$gene_symbol,
  barcodes = colData(sce_filtered)$cell_barcodes,
  version = "3",
  overwrite = TRUE
)

####################
#matrix_check <- "../results/DropletUtils/matrix.mtx"

# Read in the matrix
#my.counts <- readMM(file = matrix_check)

set.seed(2000) #ensure reproducibility for random number generation
my.counts <- sce_filtered
my.counts <- DropletUtils:::simCounts()
br.out <- barcodeRanks(my.counts)
# Check what columns are available in the output
names(br.out)

# Open a PNG device
png(output_path, width = 800, height = 600)
# Full path to save the plot 
output_path <- file.path(results_dir, "barcode_rank_plot.png")

# Plot the barcode rank plot
plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total UMI Count", main="Barcode Rank Plot")

# Order the points to plot fitted line smoothly
o <- order(br.out$rank)

# Add the fitted line to the plot
lines(br.out$rank[o], br.out$fitted[o], col="red")

# Add horizontal lines for knee and inflection points from the metadata
abline(h = metadata(br.out)$knee, col="dodgerblue", lty=2)
abline(h = metadata(br.out)$inflection, col="forestgreen", lty=2)

# Add a legend to identify the knee and inflection points
legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), legend=c("Knee", "Inflection"))
# Close the PNG device to save the file
dev.off()

# Can use the inflation value to filter hiqh quality cells(already filtered above with 100 as threshold)
inflection_value <- metadata(br.out)$inflection
inflection_value

########################
head(br.out)
head(e.out)
print(any(!is.na(e.out$LogProb)))

# Combine the relevant columns from both DataFrames
df_combined <- data.frame(
  TotalUMI = br.out$total,
  LogProb = e.out$LogProb
)

# Filter out rows where LogProb is NA
df_filtered <- df_combined[!is.na(df_combined$LogProb), ]

# Define a threshold for high-quality cells 
threshold <- inflection_value

# Add a new column to indicate high-quality cells
df_filtered$quality <- ifelse(df_filtered$TotalUMI > threshold, "High Quality", "Low Quality")

# Create scatter plot with color based on quality
ggplot(df_filtered, aes(x = TotalUMI, y = -LogProb, color = quality)) +
  geom_point() +
  scale_color_manual(values = c("High Quality" = "red", "Low Quality" = "black")) +  # Color high-quality cells red
  labs(title = "Scatter Plot of -Log Probability vs Total UMI Counts",
       x = "Total UMI Counts",
       y = "-Log Probability",
       color = "Cell Quality") +
  theme_minimal() +
  theme(legend.position = "right")  # Remove legend if you don't want it

# Full path to save the plot 
output_path <- file.path(results_dir, "scatter_plot.png")

# Save the plot using ggsave
ggsave(output_path, width = 8, height = 6, bg = "white")

