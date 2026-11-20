# ------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------

suppressPackageStartupMessages({
    library(optparse)
    library(DESeq2)
    library(ggplot2)
    library(dplyr)
})


# ------------------------------------------------------------
# 2. Command-line arguments
# ------------------------------------------------------------

option_list <- list(

    make_option(
        "--counts",
        type = "character",
        help = "Merged count matrix CSV"
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
        help = "Output DGE CSV"
    ),

    make_option(
        "--comparison",
        type = "character",
        help = "Comparison name, e.g. KO_vs_Con"
    ),

    make_option(
        "--reference",
        type = "character",
        help = "Reference group, e.g. Con"
    ),

    make_option(
        "--dge_cat",
        type = "character",
        default = "Experimental Group",
        help = "Metadata column containing experimental groups"
    ),

    make_option(
        "--filt",
        type = "logical",
        default = TRUE,
        help = "Apply count filtering"
    ),

    make_option(
        "--lfc_threshold",
        type = "numeric",
        default = 1,
        help = "Absolute log2 fold-change threshold"
    ),

    make_option(
        "--padj_threshold",
        type = "numeric",
        default = 0.05,
        help = "Adjusted p-value threshold"
    ),

    make_option(
        "--study_col",
        type = "character",
        default = NULL,
        help = "Optional batch/study column"
    )
)


opt <- parse_args(
    OptionParser(option_list = option_list)
)


# ------------------------------------------------------------
# 3. Assign arguments
# ------------------------------------------------------------

counts_file   <- opt$counts
metadata_file <- opt$metadata
annot_file    <- opt$annotation
output_file   <- opt$out

comparison_name <- opt$comparison
reference_group <- opt$reference

dge_cat <- make.names(
    opt$dge_cat
)

study_col <- if (
    !is.null(opt$study_col) &&
    length(opt$study_col) > 0 &&
    !is.na(opt$study_col) &&
    opt$study_col != ""
) {
    make.names(opt$study_col)
} else {
    NULL
}

apply_filter <- opt$filt
lfc_thresh   <- opt$lfc_threshold
padj_thresh  <- opt$padj_threshold


# ------------------------------------------------------------
# 4. Validate required arguments
# ------------------------------------------------------------

required_args <- list(
    counts = counts_file,
    metadata = metadata_file,
    annotation = annot_file,
    output = output_file,
    comparison = comparison_name,
    reference = reference_group
)

for (nm in names(required_args)) {

    value <- required_args[[nm]]

    if (
        is.null(value) ||
        length(value) == 0 ||
        is.na(value) ||
        value == ""
    ) {

        stop(
            sprintf(
                "ERROR: required argument '%s' was not supplied.",
                nm
            )
        )
    }
}


# ------------------------------------------------------------
# 5. Infer numerator from comparison name
# ------------------------------------------------------------
#
# KO_vs_Con   -> KO
# SAK3_vs_Veh -> SAK3
#
# ------------------------------------------------------------

if (!grepl("_vs_", comparison_name, fixed = TRUE)) {

    stop(
        paste0(
            "Invalid comparison name: ",
            comparison_name,
            "\nExpected format: NUMERATOR_vs_REFERENCE"
        )
    )
}

comparison_parts <- strsplit(
    comparison_name,
    "_vs_",
    fixed = TRUE
)[[1]]

if (length(comparison_parts) != 2) {

    stop(
        paste0(
            "Could not interpret comparison: ",
            comparison_name,
            "\nExpected format: NUMERATOR_vs_REFERENCE"
        )
    )
}

numerator_group <- comparison_parts[1]
comparison_reference <- comparison_parts[2]


# ------------------------------------------------------------
# 6. Validate comparison/reference consistency
# ------------------------------------------------------------

if (
    comparison_reference != reference_group
) {

    stop(
        paste0(
            "Comparison/reference mismatch.\n",
            "Comparison: ",
            comparison_name,
            "\n",
            "Reference supplied: ",
            reference_group,
            "\n",
            "Reference implied by comparison name: ",
            comparison_reference
        )
    )
}


# ------------------------------------------------------------
# 7. Create output directory
# ------------------------------------------------------------

dir.create(
    dirname(output_file),
    recursive = TRUE,
    showWarnings = FALSE
)


# ------------------------------------------------------------
# 8. Helper for safe filenames
# ------------------------------------------------------------

safe_name <- function(x) {

    gsub(
        "[^A-Za-z0-9_-]",
        "_",
        x
    )
}


# ------------------------------------------------------------
# 9. Volcano plot function
# ------------------------------------------------------------

create_volcano_plot <- function(
    df,
    contrast_name,
    outdir,
    lfc_thresh,
    padj_thresh
) {

    safe <- safe_name(
        contrast_name
    )

    df <- df %>%
        mutate(

            neg_log10_padj = ifelse(
                !is.na(padj) &
                padj > 0,
                -log10(padj),
                NA_real_
            ),

            gene_label = ifelse(
                is.na(gene_symbol) |
                gene_symbol == "",
                ensembl_gene_id,
                gene_symbol
            ),

            direction = case_when(

                significant &
                log2FoldChange > 0 ~ "Up",

                significant &
                log2FoldChange < 0 ~ "Down",

                TRUE ~ "NS"
            )
        )


    top_up <- df %>%
        filter(
            significant,
            log2FoldChange > 0
        ) %>%
        arrange(padj) %>%
        slice_head(n = 20)


    top_down <- df %>%
        filter(
            significant,
            log2FoldChange < 0
        ) %>%
        arrange(padj) %>%
        slice_head(n = 20)


    label_df <- bind_rows(
        top_up,
        top_down
    )


    p <- ggplot(
        df,
        aes(
            x = log2FoldChange,
            y = neg_log10_padj
        )
    ) +

        geom_point(
            aes(color = direction),
            alpha = 0.6,
            size = 1.2,
            na.rm = TRUE
        ) +

        scale_color_manual(
            values = c(
                "Up" = "red",
                "Down" = "blue",
                "NS" = "grey70"
            )
        ) +

        geom_vline(
            xintercept = c(
                -lfc_thresh,
                lfc_thresh
            ),
            linetype = "dashed"
        ) +

        geom_hline(
            yintercept = -log10(padj_thresh),
            linetype = "dashed"
        ) +

        geom_text(
            data = label_df,
            aes(label = gene_label),
            size = 3,
            check_overlap = TRUE
        ) +

        labs(
            title = paste(
                "Volcano:",
                contrast_name
            ),
            x = "Log2 Fold Change",
            y = "-Log10 Adjusted P-value"
        ) +

        theme_bw()


    outfile <- file.path(
        outdir,
        paste0(
            "volcano_",
            safe,
            ".png"
        )
    )


    ggsave(
        outfile,
        p,
        width = 8,
        height = 6,
        dpi = 150
    )


    message(
        "Volcano saved: ",
        outfile
    )
}


# ------------------------------------------------------------
# 10. Read metadata
# ------------------------------------------------------------

message("Reading metadata...")

metadata <- read.csv(
    metadata_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
)


# ------------------------------------------------------------
# 11. Validate Barcode
# ------------------------------------------------------------

if (
    !"Barcode" %in% colnames(metadata)
) {

    stop(
        "Metadata must contain a column named 'Barcode'."
    )
}


metadata <- metadata[
    !is.na(metadata$Barcode) &
    trimws(metadata$Barcode) != "",
    ,
    drop = FALSE
]


metadata$Barcode <- trimws(
    metadata$Barcode
)

metadata$Barcode <- make.unique(
    metadata$Barcode
)

rownames(metadata) <- metadata$Barcode

metadata$Barcode <- NULL


# Make R-safe column names
colnames(metadata) <- make.names(
    colnames(metadata)
)


# ------------------------------------------------------------
# 12. Validate metadata columns
# ------------------------------------------------------------

required_metadata <- dge_cat

if (
    !required_metadata %in% colnames(metadata)
) {

    stop(
        paste(
            "Missing metadata column:",
            required_metadata
        )
    )
}


# ------------------------------------------------------------
# 13. Read counts
# ------------------------------------------------------------

message("Reading counts...")

counts <- read.csv(
    counts_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
)


if (ncol(counts) < 2) {

    stop(
        "Count matrix contains fewer than two columns."
    )
}


# First column = gene ID
rownames(counts) <- counts[[1]]

counts <- counts[
    ,
    -1,
    drop = FALSE
]


# Clean Ensembl version suffix
rownames(counts) <- sub(
    "\\..*$",
    "",
    rownames(counts)
)


# ------------------------------------------------------------
# 14. Match samples
# ------------------------------------------------------------

common_samples <- intersect(
    colnames(counts),
    rownames(metadata)
)


if (length(common_samples) < 2) {

    stop(
        paste0(
            "Fewer than 2 samples overlap between ",
            "counts and metadata."
        )
    )
}


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


metadata <- metadata[
    colnames(counts),
    ,
    drop = FALSE
]


if (
    !identical(
        colnames(counts),
        rownames(metadata)
    )
) {

    stop(
        "Counts and metadata sample order could not be aligned."
    )
}


message(
    "Matched samples: ",
    length(common_samples)
)


# ------------------------------------------------------------
# 15. Experimental group
# ------------------------------------------------------------

metadata[[dge_cat]] <- trimws(
    as.character(
        metadata[[dge_cat]]
    )
)


# ------------------------------------------------------------
# 16. Check requested groups
# ------------------------------------------------------------

available_groups <- unique(
    metadata[[dge_cat]]
)

message(
    "Available experimental groups: ",
    paste(
        available_groups,
        collapse = ", "
    )
)

message(
    "Numerator: ",
    numerator_group
)

message(
    "Reference: ",
    reference_group
)


required_groups <- c(
    numerator_group,
    reference_group
)

missing_groups <- setdiff(
    required_groups,
    available_groups
)


if (length(missing_groups) > 0) {

    stop(
        paste0(
            "Required group(s) not present: ",
            paste(
                missing_groups,
                collapse = ", "
            )
        )
    )
}


# ------------------------------------------------------------
# 17. Restrict to requested comparison
# ------------------------------------------------------------

keep <- metadata[[dge_cat]] %in%
    required_groups

metadata_sub <- metadata[
    keep,
    ,
    drop = FALSE
]

counts_sub <- counts[
    ,
    rownames(metadata_sub),
    drop = FALSE
]


# ------------------------------------------------------------
# 18. Check replication
# ------------------------------------------------------------

group_counts <- table(
    metadata_sub[[dge_cat]]
)

message("")
message("Sample counts:")
print(group_counts)
message("")


if (
    any(
        group_counts[
            required_groups
        ] < 2
    )
) {

    stop(
        paste0(
            "Both groups require at least 2 samples for ",
            comparison_name,
            "."
        )
    )
}


# ------------------------------------------------------------
# 19. Set factor levels
# ------------------------------------------------------------

metadata_sub[[dge_cat]] <- factor(
    metadata_sub[[dge_cat]],
    levels = c(
        reference_group,
        numerator_group
    )
)


# ------------------------------------------------------------
# 20. Optional study/batch covariate
# ------------------------------------------------------------

if (!is.null(study_col)) {

    if (
        !study_col %in% colnames(metadata_sub)
    ) {

        stop(
            paste(
                "Requested study column not found:",
                study_col
            )
        )
    }


    metadata_sub[[study_col]] <- factor(
        metadata_sub[[study_col]]
    )

    metadata_sub[[study_col]] <-
        droplevels(
            metadata_sub[[study_col]]
        )


    if (
        nlevels(
            metadata_sub[[study_col]]
        ) > 1
    ) {

        design_formula <- as.formula(
            paste(
                "~",
                study_col,
                "+",
                dge_cat
            )
        )

    } else {

        design_formula <- as.formula(
            paste(
                "~",
                dge_cat
            )
        )

        message(
            "Study has one level; not included in design."
        )
    }

} else {

    design_formula <- as.formula(
        paste(
            "~",
            dge_cat
        )
    )
}


message(
    "DESeq2 design: ",
    deparse(design_formula)
)


# ------------------------------------------------------------
# 21. Create DESeq2 object
# ------------------------------------------------------------

dds <- DESeqDataSetFromMatrix(
    countData = round(
        as.matrix(
            counts_sub
        )
    ),
    colData = metadata_sub,
    design = design_formula
)


# ------------------------------------------------------------
# 22. Count filtering
# ------------------------------------------------------------

if (apply_filter) {

    keep_genes <- rowSums(
        counts(dds) >= 10
    ) >= 3

    dds <- dds[
        keep_genes,
    ]

    message(
        "Genes retained after count filter: ",
        nrow(dds)
    )
}


if (nrow(dds) == 0) {

    stop(
        "No genes remain after filtering."
    )
}


# ------------------------------------------------------------
# 23. Run DESeq2
# ------------------------------------------------------------

message("Running DESeq2...")

dds <- DESeq(
    dds
)


# ------------------------------------------------------------
# 24. Extract requested contrast
# ------------------------------------------------------------

res <- results(
    dds,
    contrast = c(
        dge_cat,
        numerator_group,
        reference_group
    ),
    alpha = padj_thresh
)


# ------------------------------------------------------------
# 25. Convert to data frame
# ------------------------------------------------------------

res_df <- as.data.frame(
    res
)

res_df$ensembl_gene_id <-
    rownames(res_df)


# ------------------------------------------------------------
# 26. Read annotation
# ------------------------------------------------------------

message("Reading annotation...")

annot <- read.delim(
    annot_file,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
)


required_annotation <- c(
    "ensembl_gene_id",
    "gene_symbol"
)

missing_annotation <- setdiff(
    required_annotation,
    colnames(annot)
)

if (length(missing_annotation) > 0) {

    stop(
        paste(
            "Annotation missing column(s):",
            paste(
                missing_annotation,
                collapse = ", "
            )
        )
    )
}


annot <- annot %>%
    mutate(
        ensembl_gene_id = trimws(
            ensembl_gene_id
        ),
        gene_symbol = trimws(
            gene_symbol
        )
    ) %>%
    filter(
        ensembl_gene_id != ""
    ) %>%
    distinct(
        ensembl_gene_id,
        .keep_all = TRUE
    )


# ------------------------------------------------------------
# 27. Add annotation
# ------------------------------------------------------------

res_df <- res_df %>%
    left_join(
        annot %>%
            select(
                ensembl_gene_id,
                gene_symbol
            ),
        by = "ensembl_gene_id"
    )


# ------------------------------------------------------------
# 28. Significance
# ------------------------------------------------------------

res_df$significant <-
    !is.na(res_df$padj) &
    res_df$padj < padj_thresh &
    !is.na(res_df$log2FoldChange) &
    abs(res_df$log2FoldChange) >= lfc_thresh


# ------------------------------------------------------------
# 29. Direction
# ------------------------------------------------------------

res_df$direction <- "NS"

res_df$direction[
    res_df$significant &
    res_df$log2FoldChange > 0
] <- "Up"

res_df$direction[
    res_df$significant &
    res_df$log2FoldChange < 0
] <- "Down"


# ------------------------------------------------------------
# 30. Add comparison information
# ------------------------------------------------------------

res_df$comparison <- comparison_name

res_df$numerator <- numerator_group

res_df$denominator <- reference_group


# ------------------------------------------------------------
# 31. Sort results
# ------------------------------------------------------------

res_df <- res_df %>%
    arrange(padj)


# ------------------------------------------------------------
# 32. Save DGE results
# ------------------------------------------------------------

write.csv(
    res_df,
    output_file,
    row.names = FALSE
)


message(
    "DGE results saved: ",
    output_file
)


# ------------------------------------------------------------
# 33. Summary
# ------------------------------------------------------------

n_up <- sum(
    res_df$significant &
    res_df$log2FoldChange > 0,
    na.rm = TRUE
)

n_down <- sum(
    res_df$significant &
    res_df$log2FoldChange < 0,
    na.rm = TRUE
)

message(
    "Significant up: ",
    n_up
)

message(
    "Significant down: ",
    n_down
)


# ------------------------------------------------------------
# 34. VST
# ------------------------------------------------------------

message("Generating VST counts...")

vst_obj <- vst(
    dds,
    blind = FALSE
)

vst_df <- as.data.frame(
    assay(vst_obj)
)

vst_df$ensembl_gene_id <-
    rownames(vst_df)


vst_file <- file.path(
    dirname(output_file),
    "vst_normalised_counts.csv"
)


write.csv(
    vst_df,
    vst_file,
    row.names = FALSE
)


message(
    "VST saved: ",
    vst_file
)


# ------------------------------------------------------------
# 35. Volcano plot
# ------------------------------------------------------------

create_volcano_plot(
    res_df,
    comparison_name,
    dirname(output_file),
    lfc_thresh,
    padj_thresh
)


# ------------------------------------------------------------
# 36. Finished
# ------------------------------------------------------------

message("")
message(
    "================================================"
)

message(
    "DGE completed successfully"
)

message(
    "Comparison: ",
    comparison_name
)

message(
    "Numerator: ",
    numerator_group
)

message(
    "Reference: ",
    reference_group
)

message(
    "================================================"
)
