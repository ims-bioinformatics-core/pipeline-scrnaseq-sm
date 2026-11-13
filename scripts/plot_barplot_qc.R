#!/usr/bin/env Rscript

# Per-sample FastQC metrics from the MultiQC data folder, coloured by metadata columns.
# FASTQ pairs of one sample (lanes / runs) are combined: reads are summed, percentages are
# weighted by reads, worst-case metrics (lowest quality, highest adapter) take the worst lane.
#
# Outputs (in --outdir):
#   qc_metrics.csv            one row per sample: metadata + QC metrics
#   barplot_qc_<column>.png   one figure per metadata column in --group_cols
#   qc_metrics_per_lane.csv   one row per FASTQ pair: sequencing run, lane + QC metrics
#   barplot_qc_run_lane.png   per-lane values coloured by run + lane, beside the --group_cols of each sample
#   stripplot_qc_groups.png   read pairs and R2 % GC across each column in --group_cols

suppressPackageStartupMessages({
    library(optparse)
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(patchwork)
})

# ---------------------------------------------------------
# Arguments
# ---------------------------------------------------------

option_list <- list(
    make_option(
        "--multiqc_dir",
        type = "character",
        help = "MultiQC data folder (multiqc_report_data)"
    ),
    make_option(
        "--filepaths",
        type = "character",
        help = "filepaths.csv (sample,fastq_1,fastq_2)"
    ),
    make_option(
        "--metadata",
        type = "character",
        help = "Metadata CSV with a Barcode column"
    ),
    make_option(
        "--group_cols",
        type = "character",
        default = "",
        help = "Comma-separated metadata columns to colour samples by"
    ),
    make_option(
        "--outdir",
        type = "character",
        help = "Output directory"
    )
)

opt <- parse_args(
    OptionParser(
        option_list = option_list
    )
)

group_cols <- trimws(strsplit(opt$group_cols, ",", fixed = TRUE)[[1]])
group_cols <- group_cols[group_cols != ""]

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

cat("MultiQC data: ", opt$multiqc_dir, "\n", sep = "")
cat("Group columns: ", paste(group_cols, collapse = ", "), "\n", sep = "")

# ---------------------------------------------------------
# FASTQ -> sample / read
# ---------------------------------------------------------

filepaths <- read.csv(opt$filepaths, check.names = FALSE, stringsAsFactors = FALSE)

fastq_base <- function(x) {
    sub("\\.(fq|fastq)(\\.gz)?$", "", basename(trimws(x)))
}

files <- bind_rows(
    data.frame(pair = seq_len(nrow(filepaths)), sample = trimws(filepaths[[1]]), file = fastq_base(filepaths[[2]]), read = "R1"),
    data.frame(pair = seq_len(nrow(filepaths)), sample = trimws(filepaths[[1]]), file = fastq_base(filepaths[[3]]), read = "R2")
)

# Run and lane from <SLX>.<sample>.<flowcell>.s_<lane>.r_<read>; other names -> "Unknown"
name_pattern <- "^.*\\.([^.]+)\\.(s_[0-9]+)\\.r_[12]$"
files$run <- ifelse(grepl(name_pattern, files$file), sub(name_pattern, "\\1", files$file), "Unknown")
files$lane <- ifelse(grepl(name_pattern, files$file), sub(name_pattern, "\\2", files$file), "Unknown")

# ---------------------------------------------------------
# MultiQC tables
# ---------------------------------------------------------

fastqc <- read.delim(
    file.path(opt$multiqc_dir, "multiqc_fastqc.txt"),
    check.names = FALSE,
    stringsAsFactors = FALSE
)

missing_files <- setdiff(files$file, fastqc$Sample)
if (length(missing_files) > 0) {
    stop(
        length(missing_files), " FASTQ files from filepaths.csv are not in MultiQC, e.g. ",
        paste(head(missing_files, 3), collapse = ", ")
    )
}

qc <- files %>%
    inner_join(fastqc, by = c("file" = "Sample")) %>%
    transmute(
        pair, sample, file, read, run, lane,
        total = `Total Sequences`,
        gc = `%GC`,
        unique_pct = total_deduplicated_percentage
    )

# "(position, value)" cells of a MultiQC line-plot table -> file, position, value
read_plot_table <- function(path) {
    if (!file.exists(path)) {
        return(NULL)
    }
    lines <- readLines(path)[-1]
    bind_rows(lapply(strsplit(lines, "\t", fixed = TRUE), function(cells) {
        pairs <- cells[-1][cells[-1] != ""]
        data.frame(
            file = cells[1],
            value = as.numeric(sub("^\\(.*, *(.*)\\)$", "\\1", pairs))
        )
    }))
}

# Lowest per-position mean quality of each file
quality <- read_plot_table(
    file.path(opt$multiqc_dir, "fastqc_per_base_sequence_quality_plot.txt")
)
if (!is.null(quality)) {
    quality <- quality %>%
        group_by(file) %>%
        summarise(min_q = min(value, na.rm = TRUE), .groups = "drop")
    qc <- left_join(qc, quality, by = "file")
}

# Highest adapter % at any position, any adapter type (rows are "<file> - <adapter>")
adapter <- read_plot_table(
    file.path(opt$multiqc_dir, "fastqc_adapter_content_plot.txt")
)
if (!is.null(adapter)) {
    adapter <- adapter %>%
        mutate(file = sub(" - [^-]*$", "", file)) %>%
        group_by(file) %>%
        summarise(max_adapter = max(value, na.rm = TRUE), .groups = "drop")
    qc <- left_join(qc, adapter, by = "file")
}

# Overrepresented sequences (absent when no file has any)
overrep_path <- file.path(opt$multiqc_dir, "fastqc_overrepresented_sequences_plot.txt")
if (file.exists(overrep_path)) {
    overrep <- read.delim(overrep_path, check.names = FALSE, stringsAsFactors = FALSE)
    overrep$overrep_pct <- rowSums(overrep[, -1, drop = FALSE], na.rm = TRUE)
    qc <- left_join(qc, overrep[, c("Sample", "overrep_pct")], by = c("file" = "Sample"))
}

for (column in c("min_q", "max_adapter", "overrep_pct")) {
    if (!column %in% colnames(qc)) {
        qc[[column]] <- NA_real_
    }
}

# ---------------------------------------------------------
# One row per sample and read
# ---------------------------------------------------------

na_if_empty <- function(x, value) {
    if (all(is.na(x))) NA_real_ else value
}

per_read <- qc %>%
    group_by(sample, read) %>%
    summarise(
        n_files = n(),
        reads = sum(total),
        gc = weighted.mean(gc, total, na.rm = TRUE),
        unique_pct = weighted.mean(unique_pct, total, na.rm = TRUE),
        overrep_pct = na_if_empty(overrep_pct, weighted.mean(overrep_pct, total, na.rm = TRUE)),
        min_q = na_if_empty(min_q, min(min_q, na.rm = TRUE)),
        max_adapter = na_if_empty(max_adapter, max(max_adapter, na.rm = TRUE)),
        .groups = "drop"
    )

wide <- per_read %>%
    pivot_wider(
        id_cols = sample,
        names_from = read,
        values_from = c(n_files, reads, gc, unique_pct, overrep_pct, min_q, max_adapter)
    ) %>%
    mutate(`Read pairs (M)` = reads_R2 / 1e6) %>%
    transmute(
        sample,
        `FASTQ pairs` = n_files_R2,
        `Read pairs (M)`,
        `R1 % GC` = gc_R1,
        `R2 % GC` = gc_R2,
        `R1 % unique reads` = unique_pct_R1,
        `R2 % unique reads` = unique_pct_R2,
        `R1 % overrepresented` = overrep_pct_R1,
        `R2 % overrepresented` = overrep_pct_R2,
        `R1 max % adapter` = max_adapter_R1,
        `R2 max % adapter` = max_adapter_R2,
        `R1 lowest mean Q` = min_q_R1,
        `R2 lowest mean Q` = min_q_R2
    )

metrics <- setdiff(colnames(wide), c("sample", "FASTQ pairs"))

# ---------------------------------------------------------
# Metadata
# ---------------------------------------------------------

metadata <- read.csv(
    opt$metadata,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    colClasses = "character",
    fileEncoding = "UTF-8-BOM"
)
colnames(metadata) <- trimws(colnames(metadata))
metadata[] <- lapply(metadata, trimws)

if (!"Barcode" %in% colnames(metadata)) {
    stop("Metadata must contain a 'Barcode' column. Found: ", paste(colnames(metadata), collapse = ", "))
}

missing_cols <- setdiff(group_cols, colnames(metadata))
if (length(missing_cols) > 0) {
    stop(
        "Columns not found in metadata: ", paste(missing_cols, collapse = ", "),
        "\nAvailable columns: ", paste(colnames(metadata), collapse = ", ")
    )
}

merged <- metadata %>%
    rename(sample = Barcode) %>%
    right_join(wide, by = "sample") %>%
    arrange(sample)

write.csv(
    rename(merged, Barcode = sample),
    file.path(opt$outdir, "qc_metrics.csv"),
    row.names = FALSE
)

merged$label <- merged$sample
if ("Sample name" %in% colnames(merged)) {
    has_name <- !is.na(merged$`Sample name`) & merged$`Sample name` != ""
    merged$label[has_name] <- paste0(merged$sample[has_name], " (", merged$`Sample name`[has_name], ")")
}

# ---------------------------------------------------------
# Plots
# ---------------------------------------------------------

okabe_ito <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00", "#F0E442", "#999999")

plot_group <- function(group_col) {
    values <- merged[[group_col]]
    values[is.na(values) | values == ""] <- "Unknown"

    # Numeric columns such as Batch: 1, 2, ... 10 rather than 1, 10, 2
    known <- setdiff(unique(values), "Unknown")
    levels_group <- if (all(!is.na(suppressWarnings(as.numeric(known))))) {
        known[order(as.numeric(known))]
    } else {
        sort(known)
    }
    if ("Unknown" %in% values) {
        levels_group <- c(levels_group, "Unknown")
    }

    colours <- if (length(levels_group) <= length(okabe_ito)) {
        okabe_ito[seq_along(levels_group)]
    } else {
        hcl.colors(length(levels_group), "Dynamic")
    }
    names(colours) <- levels_group

    # Samples ordered by group, then by depth; first sample on top
    plot_df <- merged %>%
        mutate(group = factor(values, levels = levels_group)) %>%
        arrange(group, desc(`Read pairs (M)`))
    plot_df$label <- factor(plot_df$label, levels = rev(plot_df$label))

    long <- plot_df %>%
        select(label, group, all_of(metrics)) %>%
        pivot_longer(all_of(metrics), names_to = "metric", values_to = "value") %>%
        filter(!is.na(value)) %>%
        mutate(metric = factor(metric, levels = metrics))

    n_panels <- length(unique(long$metric))
    n_rows <- ceiling(n_panels / 4)

    p <- ggplot(long, aes(x = value, y = label, fill = group)) +
        geom_col(width = 0.8) +
        facet_wrap(~metric, ncol = 4, scales = "free_x") +
        scale_fill_manual(values = colours, name = group_col) +
        labs(x = NULL, y = NULL, title = paste0("FastQC metrics per sample, by ", group_col)) +
        theme_bw(base_size = 10) +
        theme(
            panel.grid.major.y = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.y = element_text(size = 7),
            strip.background = element_rect(fill = "grey95", colour = NA),
            legend.position = "top"
        )

    output <- file.path(opt$outdir, paste0("barplot_qc_", gsub(" ", "_", group_col), ".png"))
    ggsave(
        output,
        p,
        width = 14,
        height = max(5, n_rows * (nrow(plot_df) * 0.16 + 1) + 1),
        dpi = 200,
        limitsize = FALSE
    )
    cat("Saved: ", output, "\n", sep = "")
}

invisible(lapply(group_cols, plot_group))

# ---------------------------------------------------------
# Per FASTQ pair: sequencing run and lane
# ---------------------------------------------------------

per_lane <- qc %>%
    select(-file) %>%
    pivot_wider(
        id_cols = c(pair, sample, run, lane),
        names_from = read,
        values_from = c(total, gc, unique_pct, overrep_pct, min_q, max_adapter)
    ) %>%
    transmute(
        sample,
        `Sequencing run` = run,
        Lane = lane,
        `Read pairs (M)` = total_R2 / 1e6,
        `R1 % GC` = gc_R1,
        `R2 % GC` = gc_R2,
        `R1 % unique reads` = unique_pct_R1,
        `R2 % unique reads` = unique_pct_R2,
        `R1 % overrepresented` = overrep_pct_R1,
        `R2 % overrepresented` = overrep_pct_R2,
        `R1 max % adapter` = max_adapter_R1,
        `R2 max % adapter` = max_adapter_R2,
        `R1 lowest mean Q` = min_q_R1,
        `R2 lowest mean Q` = min_q_R2
    ) %>%
    arrange(sample, `Sequencing run`, Lane)

write.csv(
    rename(per_lane, Barcode = sample),
    file.path(opt$outdir, "qc_metrics_per_lane.csv"),
    row.names = FALSE
)

clean_group <- function(x) {
    x[is.na(x) | x == ""] <- "Unknown"
    x
}

level_order <- function(values) {
    known <- setdiff(unique(values), "Unknown")
    numbers <- suppressWarnings(as.numeric(known))
    ordered <- if (all(!is.na(numbers))) known[order(numbers)] else sort(known)
    c(ordered, intersect("Unknown", values))
}

mix <- function(colour, with, amount) {
    colorRampPalette(c(colour, with))(101)[round(amount * 100) + 1]
}

# Run = hue, lane = shade within the hue
run_lane_keys <- per_lane %>%
    distinct(`Sequencing run`, Lane) %>%
    mutate(lane_number = suppressWarnings(as.numeric(sub("^s_", "", Lane)))) %>%
    arrange(`Sequencing run`, lane_number) %>%
    mutate(key = ifelse(Lane == "Unknown", `Sequencing run`, paste(`Sequencing run`, Lane)))

runs <- unique(run_lane_keys$`Sequencing run`)
run_colours <- if (length(runs) <= length(okabe_ito)) {
    okabe_ito[seq_along(runs)]
} else {
    hcl.colors(length(runs), "Dynamic")
}
names(run_colours) <- runs

run_lane_colours <- unlist(lapply(runs, function(run) {
    n <- sum(run_lane_keys$`Sequencing run` == run)
    base <- run_colours[[run]]
    if (n == 1) base else colorRampPalette(c(mix(base, "white", 0.6), base, mix(base, "black", 0.45)))(n)
}))
names(run_lane_colours) <- run_lane_keys$key

# Read pairs add up across lanes -> stacked bars; other metrics -> one point per lane.
# Group columns are shown as labelled tiles next to each sample.
plot_run_lane <- function() {
    primary <- if ("Batch" %in% group_cols) "Batch" else group_cols[1]
    order_df <- merged
    order_df$primary <- if (length(group_cols) > 0) clean_group(order_df[[primary]]) else "All"
    order_df$primary <- factor(order_df$primary, levels = level_order(order_df$primary))
    order_df <- arrange(order_df, primary, desc(`Read pairs (M)`))
    sample_levels <- rev(order_df$label)

    # Lines between the groups of the sorting column
    group_ends <- cumsum(rle(as.character(order_df$primary))$lengths)
    separators <- nrow(order_df) - group_ends[-length(group_ends)] + 0.5

    panels <- c("Read pairs (M)", "R2 % GC", "R2 % unique reads", "R2 max % adapter")
    long <- per_lane %>%
        left_join(select(merged, sample, label), by = "sample") %>%
        left_join(select(run_lane_keys, `Sequencing run`, Lane, key), by = c("Sequencing run", "Lane")) %>%
        mutate(
            key = factor(key, levels = run_lane_keys$key),
            label = factor(label, levels = sample_levels)
        ) %>%
        select(label, key, all_of(panels)) %>%
        pivot_longer(all_of(panels), names_to = "metric", values_to = "value") %>%
        filter(!is.na(value)) %>%
        mutate(metric = factor(metric, levels = panels))

    main <- ggplot(mapping = aes(x = value, y = label)) +
        geom_hline(yintercept = separators, colour = "grey75", linewidth = 0.3) +
        geom_col(
            data = filter(long, metric == "Read pairs (M)"),
            aes(fill = key),
            width = 0.8,
            colour = "white",
            linewidth = 0.15
        ) +
        geom_point(
            data = filter(long, metric != "Read pairs (M)"),
            aes(colour = key),
            size = 1.3,
            position = position_jitter(width = 0, height = 0.22, seed = 1)
        ) +
        facet_wrap(~metric, nrow = 1, scales = "free_x") +
        scale_y_discrete(limits = sample_levels) +
        scale_fill_manual(values = run_lane_colours, name = "Run and lane", drop = FALSE) +
        scale_colour_manual(values = run_lane_colours, guide = "none", drop = FALSE) +
        guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
        labs(x = NULL, y = NULL) +
        theme_bw(base_size = 10) +
        theme(
            panel.grid.major.y = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            strip.background = element_rect(fill = "grey95", colour = NA),
            legend.position = "bottom"
        )

    tiles <- bind_rows(lapply(group_cols, function(column) {
        values <- clean_group(merged[[column]])
        levels_column <- level_order(values)
        colours <- if (length(levels_column) <= length(okabe_ito)) {
            okabe_ito[seq_along(levels_column)]
        } else {
            hcl.colors(length(levels_column), "Dynamic")
        }
        data.frame(
            label = merged$label,
            column = column,
            value = values,
            fill = vapply(colours[match(values, levels_column)], mix, character(1), with = "white", amount = 0.65)
        )
    }))

    title <- "FastQC metrics per lane: each sample with its sequencing runs, lanes and groups"

    figure <- if (nrow(tiles) > 0) {
        tiles$column <- factor(tiles$column, levels = group_cols)
        tile_plot <- ggplot(tiles, aes(x = column, y = label)) +
            geom_tile(aes(fill = fill), colour = "white", linewidth = 0.4) +
            geom_text(aes(label = value), size = 2.2, colour = "grey15") +
            geom_hline(yintercept = separators, colour = "grey75", linewidth = 0.3) +
            scale_fill_identity() +
            scale_x_discrete(position = "top", expand = c(0, 0)) +
            scale_y_discrete(limits = sample_levels, expand = c(0, 0.6)) +
            labs(x = NULL, y = NULL) +
            theme_minimal(base_size = 10) +
            theme(
                panel.grid = element_blank(),
                axis.text.y = element_text(size = 7),
                axis.text.x = element_text(size = 8)
            )
        tile_plot + main + plot_layout(widths = c(length(group_cols) * 0.9, 9)) +
            plot_annotation(title = title)
    } else {
        main + plot_annotation(title = title)
    }

    output <- file.path(opt$outdir, "barplot_qc_run_lane.png")
    ggsave(
        output,
        figure,
        width = 15,
        height = max(5, nrow(merged) * 0.2 + 2),
        dpi = 200,
        limitsize = FALSE
    )
    cat("Saved: ", output, "\n", sep = "")
}

# One point per sample; line = group median; samples far from the rest are named
plot_groups <- function() {
    panels <- c("Read pairs (M)", "R2 % GC")
    all_levels <- unique(unlist(lapply(group_cols, function(column) {
        level_order(clean_group(merged[[column]]))
    })))

    long <- bind_rows(lapply(group_cols, function(column) {
        merged %>%
            transmute(sample, column = column, level = clean_group(.data[[column]]), across(all_of(panels)))
    })) %>%
        pivot_longer(all_of(panels), names_to = "metric", values_to = "value") %>%
        filter(!is.na(value)) %>%
        mutate(
            level = factor(level, levels = all_levels),
            column = factor(column, levels = group_cols),
            metric = factor(metric, levels = panels)
        ) %>%
        group_by(metric) %>%
        mutate(
            spread = mad(value),
            far = spread > 0 & abs(value - median(value)) / spread > 3.5
        ) %>%
        ungroup()

    tests <- long %>%
        group_by(column, metric) %>%
        filter(n_distinct(level) > 1) %>%
        summarise(p = kruskal.test(value ~ droplevels(level))$p.value, .groups = "drop") %>%
        mutate(text = paste0("Kruskal-Wallis p = ", signif(p, 2)))

    figure <- ggplot(long, aes(x = level, y = value)) +
        geom_point(
            colour = "grey35",
            size = 1.6,
            alpha = 0.75,
            position = position_jitter(width = 0.14, height = 0, seed = 1)
        ) +
        stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", width = 0.55, linewidth = 0.3) +
        geom_text(
            data = filter(long, far),
            aes(label = sample),
            hjust = 0,
            position = position_nudge(x = 0.2),
            size = 2.3,
            colour = "grey15"
        ) +
        geom_text(
            data = tests,
            aes(x = -Inf, y = Inf, label = text),
            hjust = -0.05,
            vjust = 1.5,
            size = 2.6,
            colour = "grey30"
        ) +
        facet_grid(metric ~ column, scales = "free", space = "free_x", switch = "y") +
        scale_x_discrete(expand = expansion(add = 0.8)) +
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
        labs(x = NULL, y = NULL, title = "Read pairs and R2 GC across sample groups") +
        theme_bw(base_size = 10) +
        theme(
            panel.grid.major.x = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_rect(fill = "grey95", colour = NA),
            strip.placement = "outside"
        )

    output <- file.path(opt$outdir, "stripplot_qc_groups.png")
    ggsave(output, figure, width = 14, height = 7, dpi = 200)
    cat("Saved: ", output, "\n", sep = "")
}

plot_run_lane()
if (length(group_cols) > 0) {
    plot_groups()
}

cat("QC metrics for", nrow(merged), " samples -> ", opt$outdir, "\n", sep = "")
