#!/usr/bin/env Rscript

# Per-sample FastQC metrics from the MultiQC data folder, coloured by metadata columns.
# FASTQ pairs of one sample (lanes / runs) are combined: reads are summed, percentages are
# weighted by reads, worst-case metrics (lowest quality, highest adapter) take the worst lane.
#
# Outputs (in --outdir):
#   qc_metrics.csv            one row per sample: metadata + QC metrics
#   qc_metrics_per_lane.csv   one row per FASTQ pair: sequencing run, lane + QC metrics
#   barplot_qc_run_lane.png   per-lane values coloured by run + lane, beside the --group_cols of each sample
#   stripplot_qc_groups.png   all metrics, R1 and R2 side by side, across each column in --group_cols

suppressPackageStartupMessages({
    library(optparse)
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(patchwork)
    library(ggrepel)
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

# Far from the other samples: more than 3.5 median absolute deviations from the median
is_outlier <- function(x) {
    spread <- mad(x, na.rm = TRUE)
    !is.na(x) & spread > 0 & abs(x - median(x, na.rm = TRUE)) / spread > 3.5
}

outlier_colour <- "#D55E00"

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

# Read pairs add up across lanes -> stacked bars; other metrics -> mean over lanes + lane range.
# Group columns are shown as labelled tiles next to each sample.
plot_run_lane <- function() {
    # Samples sorted by the --group_cols in the order given, then by read pairs (deepest first);
    # lines separate the groups of the first two columns (e.g. condition x sex)
    sort_keys <- lapply(group_cols, function(group_col) {
        values <- clean_group(merged[[group_col]])
        as.integer(factor(values, levels = level_order(values)))
    })
    order_df <- merged[do.call(order, c(sort_keys, list(-merged$`Read pairs (M)`))), ]
    sample_levels <- rev(order_df$label)

    separators <- numeric(0)
    if (length(group_cols) > 0) {
        block <- do.call(paste, c(lapply(head(group_cols, 2), function(group_col) clean_group(order_df[[group_col]])), sep = " | "))
        block_ends <- cumsum(rle(block)$lengths)
        separators <- nrow(order_df) - block_ends[-length(block_ends)] + 0.5
    }

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

    # Other metrics: one bar per sample = read-weighted mean over all runs and lanes,
    # line = lowest to highest lane, outlier samples coloured and labelled
    lane_panels <- setdiff(panels, "Read pairs (M)")
    lane_summary <- per_lane %>%
        left_join(select(merged, sample, label), by = "sample") %>%
        select(label, weight = `Read pairs (M)`, all_of(lane_panels)) %>%
        pivot_longer(all_of(lane_panels), names_to = "metric", values_to = "value") %>%
        filter(!is.na(value)) %>%
        group_by(label, metric) %>%
        summarise(
            mean = weighted.mean(value, weight),
            low = min(value),
            high = max(value),
            .groups = "drop"
        ) %>%
        group_by(metric) %>%
        mutate(outlier = is_outlier(mean), median = median(mean)) %>%
        ungroup() %>%
        mutate(
            label = factor(label, levels = sample_levels),
            metric = factor(metric, levels = panels)
        )
    medians <- distinct(lane_summary, metric, median)

    main <- ggplot(mapping = aes(x = value, y = label)) +
        geom_hline(yintercept = separators, colour = "grey45", linewidth = 0.4) +
        geom_col(
            data = filter(long, metric == "Read pairs (M)"),
            aes(fill = key),
            width = 0.8,
            colour = "white",
            linewidth = 0.15
        ) +
        geom_col(
            data = lane_summary,
            aes(x = mean),
            fill = "grey72",
            width = 0.75
        ) +
        geom_col(
            data = filter(lane_summary, outlier),
            aes(x = mean),
            fill = outlier_colour,
            width = 0.75
        ) +
        geom_vline(data = medians, aes(xintercept = median), linetype = "dashed", colour = "grey35", linewidth = 0.35) +
        geom_linerange(
            data = lane_summary,
            aes(x = mean, xmin = low, xmax = high),
            colour = "grey15",
            linewidth = 0.45
        ) +
        geom_text(
            data = filter(lane_summary, outlier),
            aes(x = high, label = sprintf("%.1f", mean)),
            hjust = -0.2,
            size = 2.4,
            colour = "grey15"
        ) +
        facet_wrap(~metric, nrow = 1, scales = "free_x") +
        scale_x_continuous(expand = expansion(mult = c(0, 0.16))) +
        scale_y_discrete(limits = sample_levels) +
        scale_fill_manual(values = run_lane_colours, name = "Run and lane", drop = FALSE) +
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

    title <- paste0("FastQC metrics per sample, sorted by ", paste(c(group_cols, "read pairs"), collapse = ", "))
    subtitle <- "Read pairs split by sequencing run and lane; other metrics averaged across runs and lanes"
    caption <- "Grey bar: read-weighted mean over all runs and lanes; black line: lowest to highest lane; dashed line: median of all samples; orange-red bar + value: more than 3.5 MADs from the median"

    figure <- if (nrow(tiles) > 0) {
        tiles$column <- factor(tiles$column, levels = group_cols)
        tile_plot <- ggplot(tiles, aes(x = column, y = label)) +
            geom_tile(aes(fill = fill), colour = "white", linewidth = 0.4) +
            geom_text(aes(label = value), size = 2.2, colour = "grey15") +
            geom_hline(yintercept = separators, colour = "grey45", linewidth = 0.4) +
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
            plot_annotation(title = title, subtitle = subtitle, caption = caption)
    } else {
        main + plot_annotation(title = title, subtitle = subtitle, caption = caption)
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

# All FastQC metrics across each group column: one point per sample, R1 and R2 side by side,
# line = median per group and read, samples far from the rest are named
plot_groups <- function() {
    read_metrics <- c("% GC", "% unique reads", "% overrepresented", "max % adapter", "lowest mean Q")
    panels <- c("Read pairs (M)", read_metrics)
    read_names <- c(R1 = "R1 (barcode + UMI)", R2 = "R2 (cDNA)")
    read_colours <- c("R1 (barcode + UMI)" = "#0072B2", "R2 (cDNA)" = "#E69F00", "Read pair" = "grey35")

    per_sample <- bind_rows(
        transmute(merged, sample, metric = "Read pairs (M)", read = "Read pair", value = `Read pairs (M)`),
        bind_rows(lapply(read_metrics, function(metric_name) {
            bind_rows(lapply(names(read_names), function(read_id) {
                source_col <- paste(read_id, metric_name)
                if (!source_col %in% colnames(merged)) {
                    return(NULL)
                }
                transmute(merged, sample, metric = metric_name, read = read_names[[read_id]], value = .data[[source_col]])
            }))
        }))
    ) %>%
        filter(!is.na(value)) %>%
        group_by(metric, read) %>%
        mutate(far = is_outlier(value)) %>%
        ungroup()

    all_levels <- unique(unlist(lapply(group_cols, function(group_col) {
        level_order(clean_group(merged[[group_col]]))
    })))

    long <- bind_rows(lapply(group_cols, function(group_col) {
        groups <- data.frame(sample = merged$sample, level = clean_group(merged[[group_col]]))
        per_sample %>%
            left_join(groups, by = "sample") %>%
            mutate(column = group_col)
    })) %>%
        mutate(
            level = factor(level, levels = all_levels),
            column = factor(column, levels = group_cols),
            metric = factor(metric, levels = panels),
            read = factor(read, levels = names(read_colours))
        )

    tests <- long %>%
        group_by(column, metric, read) %>%
        filter(n_distinct(level) > 1, n_distinct(value) > 1) %>%
        summarise(p = kruskal.test(value ~ droplevels(level))$p.value, .groups = "drop") %>%
        mutate(prefix = ifelse(read == "Read pair", "", paste0(sub(" .*", "", as.character(read)), " "))) %>%
        group_by(column, metric) %>%
        summarise(text = paste0(prefix, "p = ", signif(p, 2), collapse = "   "), .groups = "drop")

    dodge <- position_dodge(width = 0.7)

    figure <- ggplot(long, aes(x = level, y = value, colour = read)) +
        geom_point(
            size = 1.5,
            alpha = 0.8,
            position = position_jitterdodge(jitter.width = 0.12, jitter.height = 0, dodge.width = 0.7, seed = 1)
        ) +
        stat_summary(
            aes(group = read),
            fun = median,
            fun.min = median,
            fun.max = median,
            geom = "crossbar",
            width = 0.32,
            linewidth = 0.35,
            position = dodge,
            show.legend = FALSE
        ) +
        geom_text_repel(
            data = filter(long, far),
            aes(label = sample, group = read),
            position = dodge,
            size = 2.1,
            colour = "grey15",
            box.padding = 0.35,
            min.segment.length = 0,
            segment.size = 0.2,
            segment.colour = "grey50",
            max.overlaps = Inf,
            seed = 1,
            show.legend = FALSE
        ) +
        geom_text(
            data = tests,
            aes(x = -Inf, y = Inf, label = text),
            inherit.aes = FALSE,
            hjust = -0.04,
            vjust = 1.4,
            size = 2.4,
            colour = "grey30"
        ) +
        facet_grid(metric ~ column, scales = "free", space = "free_x", switch = "y") +
        scale_colour_manual(values = read_colours, name = NULL) +
        scale_x_discrete(expand = expansion(add = 0.7)) +
        scale_y_continuous(expand = expansion(mult = c(0.14, 0.22))) +
        guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
        labs(
            x = NULL,
            y = NULL,
            title = "FastQC metrics across sample groups: R1 and R2 side by side",
            caption = paste0(
                "Read pairs (M): sequenced fragments, each = one R1 read + one R2 read. ",
                "One point per sample (all runs and lanes combined); line = median.\n",
                "Named = more than 3.5 MADs from the median of all samples; Kruskal-Wallis p, unadjusted; ",
                "no point = no FastQC value (e.g. no adapter found in R1)."
            )
        ) +
        theme_bw(base_size = 10) +
        theme(
            panel.grid.major.x = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_rect(fill = "grey95", colour = NA),
            strip.placement = "outside",
            legend.position = "top",
            plot.caption = element_text(hjust = 0)
        )

    output <- file.path(opt$outdir, "stripplot_qc_groups.png")
    ggsave(output, figure, width = 14, height = 2.1 * length(unique(long$metric)) + 1.8, dpi = 200)
    cat("Saved: ", output, "\n", sep = "")
}

plot_run_lane()
if (length(group_cols) > 0) {
    plot_groups()
}

cat("QC metrics for", nrow(merged), " samples -> ", opt$outdir, "\n", sep = "")
