#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  out <- list()
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop(paste("Invalid argument:", key))
    }
    if (i == length(args)) {
      stop(paste("Missing value for argument:", key))
    }
    out[[sub("^--", "", key)]] <- args[[i + 1]]
    i <- i + 2
  }
  out
}

opt <- parse_args(args)
required <- c(
  "synonymous", "lof", "missense", "others",
  "sample_tsv", "batch_color_tsv", "out_prefix"
)
missing <- setdiff(required, names(opt))
if (length(missing) > 0) {
  stop(paste("Missing required arguments:", paste(missing, collapse = ", ")))
}

bins <- if (!is.null(opt$bins)) as.integer(opt$bins) else 120L

sample_df <- read.delim(opt$sample_tsv, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("entity.sample_id", "Batch") %in% colnames(sample_df))) {
  stop("sample.tsv must contain columns: entity.sample_id and Batch")
}
sample_df <- sample_df[, c("entity.sample_id", "Batch")]
colnames(sample_df) <- c("sample_id", "Batch")

batch_color_df <- read.delim(opt$batch_color_tsv, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("batch", "color") %in% colnames(batch_color_df))) {
  stop("CDP.batch_color.tsv must contain columns: batch and color")
}
batch_colors <- setNames(batch_color_df$color, batch_color_df$batch)

plot_one_class <- function(stat_file, class_name, out_png) {
  stat_df <- read.delim(stat_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  req <- c("sample_id", "het_count", "hom_count", "all_count")
  if (!all(req %in% colnames(stat_df))) {
    stop(paste("Stat file missing required columns:", stat_file))
  }

  stat_df$het_count <- as.numeric(stat_df$het_count)
  stat_df$hom_count <- as.numeric(stat_df$hom_count)
  stat_df$all_count <- as.numeric(stat_df$all_count)

  merged <- merge(stat_df, sample_df, by = "sample_id", all.x = TRUE, sort = FALSE)
  merged$Batch[is.na(merged$Batch) | merged$Batch == ""] <- "Unknown"

  palette <- batch_colors
  if (!("Unknown" %in% names(palette))) {
    palette <- c(palette, Unknown = "#9e9e9e")
  }

  long_df <- rbind(
    data.frame(sample_id = merged$sample_id, Batch = merged$Batch, metric = "All", count = merged$all_count, stringsAsFactors = FALSE),
    data.frame(sample_id = merged$sample_id, Batch = merged$Batch, metric = "Het", count = merged$het_count, stringsAsFactors = FALSE),
    data.frame(sample_id = merged$sample_id, Batch = merged$Batch, metric = "Hom", count = merged$hom_count, stringsAsFactors = FALSE)
  )
  long_df$metric <- factor(long_df$metric, levels = c("All", "Het", "Hom"))

  facet_layer <- facet_wrap(~metric, ncol = 1, scales = "free_y")
  if (utils::packageVersion("ggplot2") >= "3.5.0") {
    facet_layer <- facet_wrap(
      ~metric,
      ncol = 1,
      scales = "free_y",
      axes = "all_x",
      axis.labels = "all_x"
    )
  }

  p <- ggplot(long_df, aes(x = count, fill = Batch)) +
    geom_histogram(bins = bins, color = "grey20", linewidth = 0.15, position = "stack") +
    facet_layer +
    scale_fill_manual(values = palette, drop = TRUE) +
    labs(
      title = paste0(class_name, " variant count per sample"),
      x = "Variant count per sample",
      y = "Number of samples",
      fill = "Batch"
    ) +
    theme_bw() +
    guides(fill = guide_legend(ncol = 1, byrow = TRUE)) +
    theme(
      legend.position = c(0.01, 0.99),
      legend.justification = c(0, 1),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 6),
      legend.key.size = grid::unit(0.25, "cm"),
      legend.background = element_rect(fill = "white", color = "grey70"),
      legend.margin = margin(2, 2, 2, 2),
      plot.title = element_text(face = "bold")
    )

  ggsave(out_png, p, width = 9, height = 11, dpi = 200)
}

plot_one_class(opt$synonymous, "synonymous", paste0(opt$out_prefix, ".synonymous.hist.png"))
plot_one_class(opt$lof, "lof", paste0(opt$out_prefix, ".lof.hist.png"))
plot_one_class(opt$missense, "missense", paste0(opt$out_prefix, ".missense.hist.png"))
plot_one_class(opt$others, "others", paste0(opt$out_prefix, ".others.hist.png"))
