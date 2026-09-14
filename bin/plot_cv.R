#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("The R package 'ggplot2' is required")
}
if (length(args) < 3) {
  stop("Usage: plot_cv.R OUTPUT_TSV OUTPUT_PNG LOG [LOG ...]")
}

output_tsv <- args[1]
output_png <- args[2]
logs <- args[-c(1, 2)]

extract_cv <- function(path) {
  lines <- readLines(path, warn = FALSE)
  hit <- grep("CV error", lines, value = TRUE)
  if (!length(hit)) return(NULL)
  line <- tail(hit, 1)
  k <- as.integer(sub(".*K[[:space:]]*=[[:space:]]*([0-9]+).*", "\\1", line))
  err <- as.numeric(sub(".*:[[:space:]]*([0-9.eE+-]+)[[:space:]]*$", "\\1", line))
  if (is.na(k) || is.na(err)) return(NULL)
  data.frame(K = k, CV_error = err)
}

rows <- Filter(Negate(is.null), lapply(logs, extract_cv))
if (!length(rows)) {
  stop("No ADMIXTURE CV errors were found. Set admixture_cv to a value greater than zero.")
}

dat <- do.call(rbind, rows)
dat <- dat[order(dat$K), , drop = FALSE]
write.table(dat, output_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

p <- ggplot2::ggplot(dat, ggplot2::aes(x = K, y = CV_error)) +
  ggplot2::geom_line(linewidth = 0.8, colour = "#2F4858") +
  ggplot2::geom_point(size = 2.8, colour = "#2F4858") +
  ggplot2::scale_x_continuous(breaks = dat$K) +
  ggplot2::labs(
    title = "ADMIXTURE cross-validation",
    x = "K",
    y = "Cross-validation error"
  ) +
  ggplot2::theme_bw(base_size = 13) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5)
  )

ggplot2::ggsave(
  filename = output_png,
  plot = p,
  width = 8,
  height = 5.5,
  units = "in",
  dpi = 300,
  bg = "white"
)
