#!/usr/bin/env Rscript

# =========================================================
# Sex-marker expression heatmap
# =========================================================

suppressMessages(library(optparse))
suppressMessages(library(tidyverse))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(circlize))
suppressMessages(library(grid))

# =========================================================
# Command-line arguments
# =========================================================

option_list <- list(
  make_option(
    "--counts",
    type = "character",
    help = "Counts CSV"
  ),
  make_option(
    "--metadata",
    type = "character",
    help = "Metadata CSV"
  ),
  make_option(
    "--annotation",
    type = "character",
    help = "Gene annotation TSV"
  ),
  make_option(
    "--out",
    type = "character",
    help = "Output PNG"
  )
)

opt <- parse_args(
  OptionParser(
    option_list = option_list
  )
)

counts_file <- opt$counts
metadata_file <- opt$metadata
annotation_file <- opt$annotation
output_file <- opt$out

# =========================================================
# Check inputs
# =========================================================

if (
  is.null(counts_file) ||
  is.null(metadata_file) ||
  is.null(annotation_file) ||
  is.null(output_file)
) {
  
  stop(
    "Usage: Rscript script.R ",
    "--counts counts.csv ",
    "--metadata metadata.csv ",
    "--annotation annotation.tsv ",
    "--out output.png"
  )
}

# =========================================================
# Read metadata
# =========================================================

message("Reading metadata...")

metadata <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  fill = TRUE
)

# ---------------------------------------------------------
# Clean metadata column names
# ---------------------------------------------------------

metadata_names <- colnames(metadata)

bad_names <- is.na(metadata_names) |
  trimws(metadata_names) == ""

if (any(bad_names)) {
  
  metadata_names[bad_names] <- paste0(
    "Unnamed_",
    seq_len(sum(bad_names))
  )
}

metadata_names <- trimws(
  metadata_names
)

# Make all names unique
metadata_names <- make.unique(
  metadata_names
)

colnames(metadata) <- metadata_names

message(
  "Metadata columns: ",
  paste(
    colnames(metadata),
    collapse = ", "
  )
)

# =========================================================
# Find Barcode column
# =========================================================

barcode_col <- which(
  tolower(trimws(colnames(metadata))) == "barcode"
)

if (length(barcode_col) == 0) {
  
  stop(
    "Metadata must contain a 'Barcode' column.\n",
    "Available columns are: ",
    paste(
      colnames(metadata),
      collapse = ", "
    )
  )
}

barcode_col <- barcode_col[1]

# Standardize Barcode column name
metadata$Barcode <- metadata[[barcode_col]]

# =========================================================
# Clean metadata
# =========================================================

metadata$Barcode <- trimws(
  as.character(metadata$Barcode)
)

# Remove empty / missing barcodes
metadata <- metadata[
  !is.na(metadata$Barcode) &
    metadata$Barcode != "",
  ,
  drop = FALSE
]

if (nrow(metadata) == 0) {
  stop(
    "No valid Barcode values were found in metadata."
  )
}

# Make duplicate barcodes unique
metadata$Barcode <- make.unique(
  metadata$Barcode
)

# Set rownames
rownames(metadata) <- metadata$Barcode

# Remove Barcode column
metadata$Barcode <- NULL

# =========================================================
# Read counts
# =========================================================

message("Reading counts...")

counts <- read.csv(
  counts_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (ncol(counts) < 2) {
  stop(
    "Counts file must contain gene IDs plus at least one sample."
  )
}

# First column = gene ID
gene_ids <- as.character(
  counts[[1]]
)

counts <- counts[
  ,
  -1,
  drop = FALSE
]

rownames(counts) <- gene_ids

# Clean sample names
colnames(counts) <- trimws(
  colnames(counts)
)

# Remove Ensembl version suffix
rownames(counts) <- sub(
  "\\..*$",
  "",
  rownames(counts)
)

# Convert counts to numeric
counts <- as.data.frame(
  lapply(
    counts,
    function(x) {
      as.numeric(
        as.character(x)
      )
    }
  ),
  check.names = FALSE
)

rownames(counts) <- gene_ids

# =========================================================
# Match samples
# =========================================================

common_samples <- intersect(
  colnames(counts),
  rownames(metadata)
)

if (length(common_samples) == 0) {
  
  stop(
    "No matching samples found between counts and metadata.\n",
    "Counts samples: ",
    paste(
      head(colnames(counts), 10),
      collapse = ", "
    ),
    "\nMetadata barcodes: ",
    paste(
      head(rownames(metadata), 10),
      collapse = ", "
    )
  )
}

# Preserve counts-file order
common_samples <- colnames(counts)[
  colnames(counts) %in% common_samples
]

counts <- counts[
  ,
  common_samples,
  drop = FALSE
]

metadata <- metadata[
  common_samples,
  ,
  drop = FALSE
]

message(
  "Matched samples: ",
  length(common_samples)
)

# =========================================================
# Read annotation
# =========================================================

message("Reading gene annotation...")

annotation <- read.delim(
  annotation_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_annotation_cols <- c(
  "ensembl_gene_id",
  "gene_symbol"
)

missing_annotation_cols <- setdiff(
  required_annotation_cols,
  colnames(annotation)
)

if (length(missing_annotation_cols) > 0) {
  
  stop(
    "Annotation is missing required columns: ",
    paste(
      missing_annotation_cols,
      collapse = ", "
    )
  )
}

annotation$ensembl_gene_id <- sub(
  "\\..*$",
  "",
  as.character(
    annotation$ensembl_gene_id
  )
)

annotation$gene_symbol <- trimws(
  as.character(
    annotation$gene_symbol
  )
)

# Fill missing symbols with Ensembl ID
missing_symbols <- is.na(
  annotation$gene_symbol
) |
  annotation$gene_symbol == ""

annotation$gene_symbol[
  missing_symbols
] <- annotation$ensembl_gene_id[
  missing_symbols
]

# Remove duplicate Ensembl IDs
annotation <- annotation[
  !duplicated(
    annotation$ensembl_gene_id
  ),
  ,
  drop = FALSE
]

gene_id_to_symbol <- setNames(
  annotation$gene_symbol,
  annotation$ensembl_gene_id
)

# =========================================================
# Detect species
# =========================================================

n_human <- sum(
  grepl(
    "^ENSG[0-9]",
    rownames(counts)
  )
)

n_mouse <- sum(
  grepl(
    "^ENSMUSG[0-9]",
    rownames(counts)
  )
)

if (
  n_mouse > n_human &&
  n_mouse > 0
) {
  
  species <- "mouse"
  
} else if (n_human > 0) {
  
  species <- "human"
  
} else {
  
  species <- "unknown"
}

message(
  "Detected species: ",
  species,
  " (ENSG: ",
  n_human,
  ", ENSMUSG: ",
  n_mouse,
  ")"
)

# =========================================================
# Define sex markers
# =========================================================

if (species == "mouse") {
  
  female_markers <- c(
    "Xist"
  )
  
  male_markers <- c(
    "Sry",
    "Ddx3y",
    "Eif2s3y",
    "Uty",
    "Zfy1",
    "Zfy2"
  )
  
} else if (species == "human") {
  
  female_markers <- c(
    "XIST"
  )
  
  male_markers <- c(
    "SRY",
    "DDX3Y",
    "EIF2S3Y",
    "UTY",
    "ZFY"
  )
  
} else {
  
  stop(
    "Could not determine species from Ensembl gene IDs."
  )
}

sex_markers <- c(
  female_markers,
  male_markers
)

# =========================================================
# Match markers to annotation
# =========================================================

marker_annotation <- annotation %>%
  filter(
    gene_symbol %in% sex_markers
  )

if (nrow(marker_annotation) == 0) {
  
  warning(
    "No sex-marker genes were found in annotation."
  )
  
  png(
    output_file,
    width = 10,
    height = 8,
    units = "in",
    res = 300
  )
  
  plot.new()
  
  text(
    0.5,
    0.5,
    "No sex-marker genes found",
    cex = 1.5
  )
  
  dev.off()
  
  quit(
    save = "no",
    status = 0
  )
}

marker_ids <- marker_annotation$ensembl_gene_id

marker_symbols <- marker_annotation$gene_symbol

message(
  "Matched markers: ",
  paste(
    sort(unique(marker_symbols)),
    collapse = ", "
  )
)

# =========================================================
# Keep markers present in counts
# =========================================================

marker_ids_present <- marker_ids[
  marker_ids %in% rownames(counts)
]

if (length(marker_ids_present) == 0) {
  
  warning(
    "Sex-marker genes were found in annotation ",
    "but none were present in counts."
  )
  
  png(
    output_file,
    width = 10,
    height = 8,
    units = "in",
    res = 300
  )
  
  plot.new()
  
  text(
    0.5,
    0.5,
    "No sex-marker genes present in counts",
    cex = 1.5
  )
  
  dev.off()
  
  quit(
    save = "no",
    status = 0
  )
}

# =========================================================
# Extract marker expression
# =========================================================

counts_sex <- counts[
  marker_ids_present,
  ,
  drop = FALSE
]

# Map IDs to symbols
symbols <- gene_id_to_symbol[
  rownames(counts_sex)
]

symbols <- as.character(
  symbols
)

missing_symbols <- is.na(symbols) |
  symbols == ""

symbols[missing_symbols] <- rownames(
  counts_sex
)[missing_symbols]

rownames(counts_sex) <- symbols

# =========================================================
# Remove genes with zero expression
# =========================================================

keep_genes <- rowSums(
  counts_sex,
  na.rm = TRUE
) > 0

counts_sex <- counts_sex[
  keep_genes,
  ,
  drop = FALSE
]

if (nrow(counts_sex) == 0) {
  
  warning(
    "All sex-marker genes have zero expression."
  )
  
  png(
    output_file,
    width = 10,
    height = 8,
    units = "in",
    res = 300
  )
  
  plot.new()
  
  text(
    0.5,
    0.5,
    "All sex-marker genes have zero expression",
    cex = 1.5
  )
  
  dev.off()
  
  quit(
    save = "no",
    status = 0
  )
}

# =========================================================
# Collapse duplicated gene symbols
# =========================================================

counts_sex <- rowsum(
  as.matrix(counts_sex),
  group = rownames(counts_sex)
)

# =========================================================
# Log transform
# =========================================================

counts_sex_log <- log10(
  counts_sex + 1
)

# =========================================================
# Find Sex column
# =========================================================

sex_col <- intersect(
  c(
    "Sex",
    "sex",
    "Sample Sex",
    "Sample sex"
  ),
  colnames(metadata)
)[1]

if (is.na(sex_col)) {
  
  warning(
    "No Sex column found in metadata. ",
    "All samples will be labelled Unknown."
  )
  
  metadata$sex <- "Unknown"
  
} else {
  
  # IMPORTANT:
  # Correct extraction of the selected metadata column
  metadata$sex <- metadata[[sex_col]]
  
  metadata$sex <- trimws(
    as.character(
      metadata$sex
    )
  )
  
  metadata$sex <- case_when(
    
    metadata$sex %in% c(
      "Male",
      "male",
      "M",
      "m"
    ) ~ "Male",
    
    metadata$sex %in% c(
      "Female",
      "female",
      "F",
      "f"
    ) ~ "Female",
    
    TRUE ~ "Unknown"
  )
}

metadata$sex[
  is.na(metadata$sex)
] <- "Unknown"

metadata$sex <- factor(
  metadata$sex,
  levels = c(
    "Male",
    "Female",
    "Unknown"
  )
)

# =========================================================
# Align sex annotation to heatmap columns
# =========================================================

sample_sex <- metadata$sex

names(sample_sex) <- rownames(metadata)

sample_sex <- sample_sex[
  colnames(counts_sex_log)
]

# =========================================================
# Sex annotation colours
# =========================================================

sex_colors <- c(
  Male = "#2c7fb8",
  Female = "#e41a1c",
  Unknown = "grey70"
)

ha <- HeatmapAnnotation(
  
  Sex = sample_sex,
  
  col = list(
    Sex = sex_colors
  ),
  
  annotation_name_gp = gpar(
    fontsize = 11,
    fontface = "bold"
  ),
  
  simple_anno_size = unit(
    5,
    "mm"
  )
)

# =========================================================
# Plot dimensions
#
# Keep width moderate so labels do not become tiny.
# =========================================================

n_samples <- ncol(
  counts_sex_log
)

n_genes <- nrow(
  counts_sex_log
)

# Moderate width
plot_width <- max(
  10,
  min(
    18,
    7 + n_samples * 0.22
  )
)

# Height
plot_height <- max(
  8,
  5 + n_genes * 0.45
)

# Readable sample-name font
column_fontsize <- if (n_samples > 50) {
  
  7
  
} else if (n_samples > 35) {
  
  8
  
} else if (n_samples > 20) {
  
  9
  
} else {
  
  10
}

# Readable gene-name font
row_fontsize <- if (n_genes > 20) {
  
  8
  
} else if (n_genes > 12) {
  
  9
  
} else {
  
  10
}

message(
  "Samples: ",
  n_samples
)

message(
  "Genes: ",
  n_genes
)

message(
  "Plot size: ",
  round(plot_width, 1),
  " x ",
  round(plot_height, 1),
  " inches"
)

message(
  "Sample label font size: ",
  column_fontsize
)

# =========================================================
# Create heatmap
# =========================================================

png(
  output_file,
  width = plot_width,
  height = plot_height,
  units = "in",
  res = 300
)

# Handle case where all values are identical
heatmap_min <- min(
  counts_sex_log,
  na.rm = TRUE
)

heatmap_max <- max(
  counts_sex_log,
  na.rm = TRUE
)

if (heatmap_min == heatmap_max) {
  
  heatmap_max <- heatmap_min + 1
}

ht <- Heatmap(
  
  counts_sex_log,
  
  name = "log10(count+1)",
  
  bottom_annotation = ha,
  
  cluster_rows = TRUE,
  
  cluster_columns = TRUE,
  
  show_row_names = TRUE,
  
  show_column_names = TRUE,
  
  column_names_rot = 45,
  
  column_names_gp = gpar(
    fontsize = column_fontsize
  ),
  
  row_names_gp = gpar(
    fontsize = row_fontsize
  ),
  
  col = colorRamp2(
    c(
      heatmap_min,
      heatmap_max
    ),
    c(
      "white",
      "red"
    )
  )
)

# =========================================================
# Draw
# =========================================================

draw(
  ht,
  padding = unit(
    c(
      8,
      8,
      50,
      8
    ),
    "mm"
  )
)

dev.off()

# =========================================================
# Finished
# =========================================================

message("")
message("==========================================")
message("Sex-marker heatmap saved:")
message(output_file)
message("==========================================")