#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(optparse)
    library(dplyr)
    library(ggplot2)
    library(clusterProfiler)
    library(org.Mm.eg.db)
    library(org.Hs.eg.db)
    library(AnnotationDbi)
    library(ReactomePA)
})


# =========================================================
# COMMAND-LINE OPTIONS
# =========================================================

option_list <- list(

    make_option(
        "--dge",
        type = "character",
        help = "DGE results CSV file"
    ),

    make_option(
        "--annotation",
        type = "character",
        default = NULL,
        help = "Annotation file"
    ),

    make_option(
        "--outdir",
        type = "character",
        help = "Output directory"
    ),

    make_option(
        "--comparison",
        type = "character",
        help = "Comparison name"
    ),

    make_option(
        "--organism",
        type = "character",
        default = "mouse",
        help = "Organism: mouse or human [default: mouse]"
    ),

    make_option(
        "--databases",
        type = "character",
        default = "GO,KEGG,REACTOME",
        help = "Comma-separated databases"
    ),

    make_option(
        "--min_genes",
        type = "integer",
        default = 10,
        help = "Minimum gene-set size"
    ),

    make_option(
        "--max_genes",
        type = "integer",
        default = 500,
        help = "Maximum gene-set size"
    ),

    make_option(
        "--gsea_fdr",
        type = "double",
        default = 0.05,
        help = "FDR threshold"
    )
)


opt <- parse_args(
    OptionParser(
        option_list = option_list
    )
)


# =========================================================
# VALIDATE ARGUMENTS
# =========================================================

if (is.null(opt$dge)) {
    stop("--dge is required")
}

if (is.null(opt$outdir)) {
    stop("--outdir is required")
}

if (is.null(opt$comparison)) {
    stop("--comparison is required")
}


# =========================================================
# CREATE OUTPUT DIRECTORY
# =========================================================

if (!dir.exists(opt$outdir)) {

    dir.create(
        opt$outdir,
        recursive = TRUE,
        showWarnings = FALSE
    )
}


# =========================================================
# ORGANISM
# =========================================================

organism <- tolower(
    opt$organism
)

if (!organism %in% c("mouse", "human")) {

    stop(
        "Supported organisms: mouse, human. Got: ", organism,
        call. = FALSE
    )
}

# Organism-specific databases: gene IDs (OrgDb), KEGG code, Reactome name
org_db       <- if (organism == "mouse") org.Mm.eg.db else org.Hs.eg.db
kegg_code    <- if (organism == "mouse") "mmu" else "hsa"
reactome_org <- organism


# =========================================================
# READ DGE RESULTS
# =========================================================

message(
    "Reading DGE results: ",
    opt$dge
)

dge <- read.csv(
    opt$dge,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

if (nrow(dge) == 0) {

    stop(
        "DGE file contains no rows."
    )
}


message(
    "DGE rows: ",
    nrow(dge)
)


# =========================================================
# CHECK REQUIRED DGE COLUMNS
# =========================================================

if (!"gene_symbol" %in% colnames(dge)) {

    stop(
        "The DGE results table does not contain a 'gene_symbol' column.\n",
        "Available columns are:\n",
        paste(
            colnames(dge),
            collapse = ", "
        )
    )
}


message(
    "Using gene_symbol column from DGE results."
)


# =========================================================
# IDENTIFY STATISTIC COLUMN
# =========================================================

stat_candidates <- c(
    "stat",
    "statistic"
)

stat_matches <- stat_candidates[
    stat_candidates %in% colnames(dge)
]

if (length(stat_matches) == 0) {

    stop(
        "Could not find a DESeq2 statistic column. ",
        "Expected one of: ",
        paste(
            stat_candidates,
            collapse = ", "
        )
    )
}

stat_col <- stat_matches[1]


# =========================================================
# IDENTIFY GENE COLUMN
# =========================================================

gene_candidates <- c(
    "gene",
    "Gene",
    "gene_id",
    "GeneID",
    "geneID",
    "ensembl_gene_id",
    "ENSEMBL"
)

gene_matches <- gene_candidates[
    gene_candidates %in% colnames(dge)
]


if (length(gene_matches) == 0) {

    stop(
        "Could not find an Ensembl gene column. ",
        "Expected one of: ",
        paste(
            gene_candidates,
            collapse = ", "
        )
    )
}

gene_col <- gene_matches[1]


message(
    "Using gene column: ",
    gene_col
)

message(
    "Using statistic column: ",
    stat_col
)


# =========================================================
# PREPARE RANKING DATA
# =========================================================

ranking_df <- dge %>%

    dplyr::mutate(

        Ensembl = as.character(
            .data[[gene_col]]
        ),

        gene_symbol = as.character(
            gene_symbol
        ),

        statistic = suppressWarnings(
            as.numeric(
                .data[[stat_col]]
            )
        )
    ) %>%

    dplyr::filter(

        !is.na(Ensembl),

        Ensembl != "",

        !is.na(statistic),

        is.finite(statistic)
    ) %>%

    dplyr::mutate(

        # Remove Ensembl version suffix
        #
        # ENSMUSG00000012345.7
        #
        # becomes
        #
        # ENSMUSG00000012345

        Ensembl = sub(
            "\\..*$",
            "",
            Ensembl
        )
    ) %>%

    dplyr::filter(
        Ensembl != ""
    )


if (nrow(ranking_df) == 0) {

    stop(
        "No valid genes/statistics available for GSEA."
    )
}


# =========================================================
# REMOVE DUPLICATE ENSEMBL IDS
#
# Keep the row with the largest absolute statistic.
# =========================================================

ranking_df <- ranking_df %>%

    dplyr::group_by(
        Ensembl
    ) %>%

    dplyr::slice_max(

        order_by = abs(statistic),

        n = 1,

        with_ties = FALSE
    ) %>%

    dplyr::ungroup()


# =========================================================
# MAP ENSEMBL -> ENTREZ
#
# IMPORTANT:
#
# gene_symbol already exists in the DGE results table.
# Therefore we deliberately do NOT request SYMBOL here.
# =========================================================

message(
    "Mapping Ensembl IDs to Entrez IDs..."
)


gene_map <- AnnotationDbi::select(

    org_db,

    keys = unique(
        ranking_df$Ensembl
    ),

    keytype = "ENSEMBL",

    columns = c(
        "ENSEMBL",
        "ENTREZID"
    )
)


# =========================================================
# CLEAN ANNOTATION TABLE
# =========================================================

gene_map <- gene_map %>%

    dplyr::filter(

        !is.na(ENSEMBL),

        !is.na(ENTREZID),

        ENTREZID != ""
    ) %>%

    dplyr::distinct(

        ENSEMBL,

        ENTREZID
    )


if (nrow(gene_map) == 0) {

    stop(
        "No Ensembl IDs could be mapped to Entrez IDs."
    )
}


# =========================================================
# JOIN ANNOTATION TO RANKING
# =========================================================

ranking_df <- ranking_df %>%

    dplyr::inner_join(

        gene_map,

        by = c(
            "Ensembl" = "ENSEMBL"
        )
    )


if (nrow(ranking_df) == 0) {

    stop(
        "No genes remained after Ensembl -> Entrez mapping."
    )
}


# =========================================================
# VERIFY GENE SYMBOL SURVIVED JOIN
# =========================================================

if (!"gene_symbol" %in% colnames(ranking_df)) {

    stop(
        "The 'gene_symbol' column was lost during annotation join."
    )
}


# =========================================================
# REMOVE DUPLICATE ENTREZ IDS
#
# Keep the gene with the largest absolute statistic.
# =========================================================

ranking_df <- ranking_df %>%

    dplyr::group_by(
        ENTREZID
    ) %>%

    dplyr::slice_max(

        order_by = abs(statistic),

        n = 1,

        with_ties = FALSE
    ) %>%

    dplyr::ungroup()


# =========================================================
# CREATE ENTREZ -> GENE SYMBOL MAP
#
# This is used ONLY for displaying core_enrichment.
#
# GSEA itself continues to use Entrez IDs.
# =========================================================

entrez_to_symbol <- ranking_df %>%

    dplyr::select(
        ENTREZID,
        gene_symbol
    ) %>%

    dplyr::filter(
        !is.na(ENTREZID),
        ENTREZID != "",
        !is.na(gene_symbol),
        gene_symbol != ""
    ) %>%

    dplyr::distinct(
        ENTREZID,
        .keep_all = TRUE
    )


# =========================================================
# CREATE RANKED GENE LIST
# =========================================================

gene_list <- ranking_df$statistic

names(gene_list) <- as.character(
    ranking_df$ENTREZID
)

gene_list <- sort(
    gene_list,
    decreasing = TRUE
)


if (length(gene_list) < opt$min_genes) {

    stop(
        "Too few mapped genes for GSEA: ",
        length(gene_list)
    )
}


message(
    "Genes used for GSEA: ",
    length(gene_list)
)


# =========================================================
# SAVE GSEA GENE RANKING
# =========================================================

ranking_output <- ranking_df %>%

    dplyr::arrange(
        dplyr::desc(statistic)
    ) %>%

    dplyr::select(

        Ensembl,

        ENTREZID,

        gene_symbol,

        statistic
    )


ranking_file <- file.path(

    opt$outdir,

    "GSEA_gene_ranking.csv"
)


write.csv(

    ranking_output,

    ranking_file,

    row.names = FALSE
)


# =========================================================
# DETERMINE REQUESTED DATABASES
# =========================================================

databases <- toupper(

    trimws(

        unlist(

            strsplit(

                opt$databases,

                ","
            )
        )
    )
)


databases <- unique(

    databases[
        databases %in%
            c(
                "GO",
                "KEGG",
                "REACTOME"
            )
    ]
)


if (length(databases) == 0) {

    stop(

        "No supported databases requested. ",

        "Supported databases: GO, KEGG, REACTOME"
    )
}


message(

    "Requested databases: ",

    paste(

        databases,

        collapse = ", "
    )
)


# =========================================================
# GO GSEA
# =========================================================

go_result <- NULL


if ("GO" %in% databases) {

    message(
        "Running GO GSEA..."
    )

    go_result <- tryCatch(

        {

            clusterProfiler::gseGO(

                geneList = gene_list,

                OrgDb = org_db,

                keyType = "ENTREZID",

                ont = "BP",

                minGSSize = opt$min_genes,

                maxGSSize = opt$max_genes,

                pvalueCutoff = 1,

                pAdjustMethod = "BH",

                verbose = FALSE
            )

        },

        error = function(e) {

            message(

                "GO GSEA failed: ",

                conditionMessage(e)
            )

            NULL
        }
    )
}


# =========================================================
# KEGG GSEA
# =========================================================

kegg_result <- NULL


if ("KEGG" %in% databases) {

    message(
        "Running KEGG GSEA..."
    )

    kegg_result <- tryCatch(

        {

            clusterProfiler::gseKEGG(

                geneList = gene_list,

                organism = kegg_code,

                keyType = "ncbi-geneid",

                minGSSize = opt$min_genes,

                maxGSSize = opt$max_genes,

                pvalueCutoff = 1,

                pAdjustMethod = "BH",

                verbose = FALSE
            )

        },

        error = function(e) {

            message(

                "KEGG GSEA failed: ",

                conditionMessage(e)
            )

            NULL
        }
    )
}


# =========================================================
# REACTOME GSEA
# =========================================================

reactome_result <- NULL


if ("REACTOME" %in% databases) {

    message(
        "Running Reactome GSEA..."
    )

    reactome_result <- tryCatch(

        {

            ReactomePA::gsePathway(

                geneList = gene_list,

                organism = reactome_org,

                pvalueCutoff = 1,

                minGSSize = opt$min_genes,

                maxGSSize = opt$max_genes,

                pAdjustMethod = "BH",

                verbose = FALSE
            )

        },

        error = function(e) {

            message(

                "Reactome GSEA failed: ",

                conditionMessage(e)
            )

            NULL
        }
    )
}


# =========================================================
# CONVERT RESULTS TO DATA FRAMES
# =========================================================

go_df <- NULL

kegg_df <- NULL

reactome_df <- NULL


if (!is.null(go_result)) {

    go_df <- as.data.frame(
        go_result
    )
}


if (!is.null(kegg_result)) {

    kegg_df <- as.data.frame(
        kegg_result
    )
}


if (!is.null(reactome_result)) {

    reactome_df <- as.data.frame(
        reactome_result
    )
}


# =========================================================
# CONVERT CORE ENRICHMENT ENTREZ IDS
# TO GENE SYMBOLS
#
# clusterProfiler stores core_enrichment as:
#
#   12575/17869/22059/...
#
# We convert this to:
#
#   Cdkn1a/Myc/Trp53/...
#
# while leaving the actual GSEA calculation unchanged.
# =========================================================

convert_core_enrichment_to_symbols <- function(

    df,

    entrez_map

) {

    if (

        is.null(df) ||

        nrow(df) == 0 ||

        !"core_enrichment" %in% colnames(df)

    ) {

        return(df)
    }


    symbol_lookup <- setNames(

        entrez_map$gene_symbol,

        as.character(
            entrez_map$ENTREZID
        )
    )


    df$core_enrichment <- vapply(

        strsplit(

            as.character(
                df$core_enrichment
            ),

            "/",

            fixed = TRUE
        ),

        function(ids) {

            ids <- trimws(ids)

            symbols <- unname(
                symbol_lookup[ids]
            )

            # If an Entrez ID has no symbol in the
            # DGE-derived mapping, retain the ID rather
            # than silently dropping it.

            output <- ifelse(

                !is.na(symbols) & symbols != "",

                symbols,

                ids
            )

            paste(

                output,

                collapse = "/"
            )
        },

        FUN.VALUE = character(1)
    )


    df
}


# =========================================================
# REPLACE CORE ENRICHMENT IDS WITH GENE SYMBOLS
# =========================================================

go_df <- convert_core_enrichment_to_symbols(

    go_df,

    entrez_to_symbol
)


kegg_df <- convert_core_enrichment_to_symbols(

    kegg_df,

    entrez_to_symbol
)


reactome_df <- convert_core_enrichment_to_symbols(

    reactome_df,

    entrez_to_symbol
)


# =========================================================
# ADD SIGNIFICANCE + DIRECTION
#
# Significance:
#
#   p.adjust <= gsea_fdr
#
# Direction:
#
#   NES > 0 = UP
#   NES < 0 = DOWN
# =========================================================

add_gsea_annotations <- function(

    df,

    fdr_threshold

) {

    if (

        is.null(df) ||

        nrow(df) == 0

    ) {

        return(df)
    }


    df %>%

        dplyr::mutate(

            NES = as.numeric(
                NES
            ),

            p.adjust = as.numeric(
                p.adjust
            ),

            Significance = dplyr::case_when(

                !is.na(p.adjust) &

                    p.adjust <= fdr_threshold ~

                    "Significant",

                !is.na(p.adjust) &

                    p.adjust > fdr_threshold ~

                    "Not significant",

                TRUE ~

                    NA_character_
            ),

            Direction = dplyr::case_when(

                !is.na(NES) &

                    NES > 0 ~

                    "UP",

                !is.na(NES) &

                    NES < 0 ~

                    "DOWN",

                !is.na(NES) &

                    NES == 0 ~

                    "NEUTRAL",

                TRUE ~

                    NA_character_
            )
        )
}


go_df <- add_gsea_annotations(

    go_df,

    opt$gsea_fdr
)


kegg_df <- add_gsea_annotations(

    kegg_df,

    opt$gsea_fdr
)


reactome_df <- add_gsea_annotations(

    reactome_df,

    opt$gsea_fdr
)


# =========================================================
# SAVE GO RESULTS
# =========================================================

go_file <- file.path(

    opt$outdir,

    "GO_GSEA.csv"
)


if (!is.null(go_df)) {

    write.csv(

        go_df,

        go_file,

        row.names = FALSE
    )

} else {

    write.csv(

        data.frame(),

        go_file,

        row.names = FALSE
    )
}


# =========================================================
# SAVE KEGG RESULTS
# =========================================================

kegg_file <- file.path(

    opt$outdir,

    "KEGG_GSEA.csv"
)


if (!is.null(kegg_df)) {

    write.csv(

        kegg_df,

        kegg_file,

        row.names = FALSE
    )

} else {

    write.csv(

        data.frame(),

        kegg_file,

        row.names = FALSE
    )
}


# =========================================================
# SAVE REACTOME RESULTS
# =========================================================

reactome_file <- file.path(

    opt$outdir,

    "REACTOME_GSEA.csv"
)


if (!is.null(reactome_df)) {

    write.csv(

        reactome_df,

        reactome_file,

        row.names = FALSE
    )

} else {

    write.csv(

        data.frame(),

        reactome_file,

        row.names = FALSE
    )
}


# =========================================================
# PREPARE DOTPLOT DATA
# =========================================================

plot_list <- list()


if (

    !is.null(go_df) &&

    nrow(go_df) > 0

) {

    plot_list[[length(plot_list) + 1]] <-

        go_df %>%

        dplyr::mutate(

            Database = "GO"
        )
}


if (

    !is.null(kegg_df) &&

    nrow(kegg_df) > 0

) {

    plot_list[[length(plot_list) + 1]] <-

        kegg_df %>%

        dplyr::mutate(

            Database = "KEGG"
        )
}


if (

    !is.null(reactome_df) &&

    nrow(reactome_df) > 0

) {

    plot_list[[length(plot_list) + 1]] <-

        reactome_df %>%

        dplyr::mutate(

            Database = "REACTOME"
        )
}


plot_file <- file.path(

    opt$outdir,

    "pathway_dotplot.png"
)


# =========================================================
# NO GSEA RESULTS
# =========================================================

if (length(plot_list) == 0) {

    message(
        "No GSEA results available."
    )


    p <- ggplot() +

        annotate(

            "text",

            x = 1,

            y = 1,

            label = "No GSEA results available",

            size = 6
        ) +

        theme_void()


    ggsave(

        filename = plot_file,

        plot = p,

        width = 18,

        height = 8,

        dpi = 300
    )


} else {


    # =====================================================
    # COMBINE DATABASE RESULTS
    # =====================================================

    combined <- dplyr::bind_rows(

        plot_list
    )


    combined <- combined %>%

        dplyr::mutate(

            NES = as.numeric(
                NES
            ),

            p.adjust = as.numeric(
                p.adjust
            )
        ) %>%

        dplyr::filter(

            !is.na(NES),

            is.finite(NES),

            !is.na(p.adjust),

            is.finite(p.adjust)
        )


    # =====================================================
    # SIGNIFICANT PATHWAYS ONLY
    # =====================================================

    combined_sig <- combined %>%

        dplyr::filter(

            p.adjust <= opt$gsea_fdr
        )


    message(

        "Significant pathways at FDR <= ",

        opt$gsea_fdr,

        ": ",

        nrow(combined_sig)
    )


    # =====================================================
    # TOP 10 POSITIVE NES PER DATABASE
    # =====================================================

    up_pathways <- combined_sig %>%

        dplyr::filter(

            NES > 0
        ) %>%

        dplyr::group_by(

            Database
        ) %>%

        dplyr::arrange(

            dplyr::desc(NES),

            .by_group = TRUE
        ) %>%

        dplyr::slice_head(

            n = 10
        ) %>%

        dplyr::ungroup() %>%

        dplyr::mutate(

            Direction = "UP"
        )


    # =====================================================
    # TOP 10 NEGATIVE NES PER DATABASE
    # =====================================================

    down_pathways <- combined_sig %>%

        dplyr::filter(

            NES < 0
        ) %>%

        dplyr::group_by(

            Database
        ) %>%

        dplyr::arrange(

            NES,

            .by_group = TRUE
        ) %>%

        dplyr::slice_head(

            n = 10
        ) %>%

        dplyr::ungroup() %>%

        dplyr::mutate(

            Direction = "DOWN"
        )


    # =====================================================
    # COMBINE UP + DOWN
    # =====================================================

    combined_sig <- dplyr::bind_rows(

        up_pathways,

        down_pathways
    )


    if (nrow(combined_sig) == 0) {

        message(
            "No positive or negative significant NES pathways available."
        )


        p <- ggplot() +

            annotate(

                "text",

                x = 1,

                y = 1,

                label = paste0(

                    "No significant pathways at FDR <= ",

                    opt$gsea_fdr
                ),

                size = 6
            ) +

            theme_void()


        ggsave(

            filename = plot_file,

            plot = p,

            width = 18,

            height = 8,

            dpi = 300
        )


    } else {


        # =================================================
        # PATHWAY LABEL
        # =================================================

        combined_sig <- combined_sig %>%

            dplyr::mutate(

                Pathway = dplyr::case_when(

                    !is.na(Description) &

                        Description != "" ~

                        Description,

                    !is.na(ID) &

                        ID != "" ~

                        ID,

                    TRUE ~

                        "Unknown pathway"
                )
            )


        # =================================================
        # DOT SIZE = -LOG10(FDR)
        # =================================================

        combined_sig <- combined_sig %>%

            dplyr::mutate(

                neg_log10_fdr =

                    -log10(

                        pmax(

                            p.adjust,

                            .Machine$double.xmin
                        )
                    )
            )


        # =================================================
        # MAKE PATHWAY NAMES UNIQUE
        # =================================================

        combined_sig <- combined_sig %>%

            dplyr::mutate(

                Pathway = make.unique(

                    as.character(

                        Pathway
                    )
                )
            )


        # =================================================
        # ORDER PATHWAYS WITHIN DATABASE
        # =================================================

        combined_sig <- combined_sig %>%

            dplyr::group_by(

                Database
            ) %>%

            dplyr::arrange(

                NES,

                .by_group = TRUE
            ) %>%

            dplyr::mutate(

                Pathway = factor(

                    Pathway,

                    levels = unique(

                        Pathway
                    )
                )
            ) %>%

            dplyr::ungroup()


        # =================================================
        # DOTPLOT
        # =================================================

        p <- ggplot(

            combined_sig,

            aes(

                x = NES,

                y = Pathway,

                size = neg_log10_fdr,

                colour = Direction
            )

        ) +

            geom_vline(

                xintercept = 0,

                linetype = "dashed"
            ) +

            geom_point(

                alpha = 0.85
            ) +

            facet_wrap(

                ~Database,

                scales = "free_y",

                ncol = 1
            ) +

            labs(

                title = paste(

                    "GSEA Pathway Analysis:",

                    opt$comparison
                ),

                x = "Normalized Enrichment Score (NES)",

                y = NULL,

                size = expression(

                    -log[10]("FDR")
                ),

                colour = "Direction"
            ) +

            scale_size_continuous(

                range = c(

                    2.5,

                    10
                )
            ) +

            scale_colour_manual(

                values = c(

                    "UP" = "#D73027",

                    "DOWN" = "#4575B4"
                )
            ) +

            theme_bw() +

            theme(

                plot.title =

                    element_text(

                        face = "bold",

                        size = 14
                    ),

                axis.text.y =

                    element_text(

                        size = 9
                    ),

                axis.text.x =

                    element_text(

                        size = 10
                    ),

                axis.title.x =

                    element_text(

                        size = 11
                    ),

                strip.text =

                    element_text(

                        face = "bold",

                        size = 12
                    ),

                legend.title =

                    element_text(

                        face = "bold"
                    ),

                plot.margin =

                    margin(

                        10,

                        20,

                        10,

                        30
                    ),

                panel.spacing =

                    unit(

                        1.5,

                        "lines"
                    )
            )


        # =================================================
        # SAVE DOTPLOT
        # =================================================

        ggsave(

            filename = plot_file,

            plot = p,

            width = 18,

            height = 14,

            dpi = 300
        )
    }
}


# =========================================================
# FINAL OUTPUT MESSAGE
# =========================================================

message("")

message(
    "GSEA analysis completed successfully."
)

message(
    "Output directory: ",
    opt$outdir
)

message(
    "Gene ranking: ",
    ranking_file
)

message(
    "GO results: ",
    go_file
)

message(
    "KEGG results: ",
    kegg_file
)

message(
    "Reactome results: ",
    reactome_file
)

message(
    "Dotplot: ",
    plot_file
)