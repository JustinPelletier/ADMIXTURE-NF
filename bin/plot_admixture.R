#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("The R package 'ggplot2' is required")
}

if (!requireNamespace("patchwork", quietly = TRUE)) {
  stop("The R package 'patchwork' is required")
}

if (length(args) != 9) {
  stop(paste(
    "Usage: plot_admixture.R Q FAM K CAG_METADATA REF_METADATA",
    "CAG_KEEP REF_KEEP OUTPUT_PNG OUTPUT_ORDER"
  ))
}


# -------------------------------------------------------------------------
# Arguments
# -------------------------------------------------------------------------

q_path <- args[1]
fam_path <- args[2]
k <- as.integer(args[3])
cag_meta_path <- args[4]
ref_meta_path <- args[5]
cag_keep_path <- args[6]
ref_keep_path <- args[7]
output_png <- args[8]
output_order <- args[9]


# -------------------------------------------------------------------------
# Input functions
# -------------------------------------------------------------------------

read_keep <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    return(character())
  }

  x <- readLines(path, warn = FALSE)
  x <- trimws(sub("\t.*$", "", x))
  x <- x[nzchar(x) & !grepl("^#", x)]

  # Remove a possible header.
  x[!tolower(x) %in% c("sample_id", "iid", "id")]
}


read_metadata <- function(path, cohort_name) {
  empty_metadata <- data.frame(
    sample_id = character(),
    population = character(),
    cohort = character(),
    metadata_order = integer(),
    stringsAsFactors = FALSE
  )

  if (!file.exists(path) || file.info(path)$size == 0) {
    return(empty_metadata)
  }

  x <- tryCatch(
    read.delim(
      path,
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      comment.char = ""
    ),
    error = function(e) {
      stop(
        "Could not read metadata file ",
        path,
        ": ",
        conditionMessage(e)
      )
    }
  )

  if (nrow(x) == 0) {
    return(empty_metadata)
  }

  required <- c("sample_id", "population")

  if (!all(required %in% names(x))) {
    stop(path, " must contain columns: sample_id and population")
  }

  x <- x[, required, drop = FALSE]

  x$sample_id <- trimws(as.character(x$sample_id))
  x$population <- trimws(as.character(x$population))

  x <- x[
    nzchar(x$sample_id) &
      nzchar(x$population),
    ,
    drop = FALSE
  ]

  if (anyDuplicated(x$sample_id)) {
    stop("Duplicate sample_id in ", path)
  }

  x$cohort <- cohort_name
  x$metadata_order <- seq_len(nrow(x))

  x
}


# Numerically order labels such as:
# Cluster_1, Cluster_2, ..., Cluster_10
order_cluster_levels <- function(x) {
  cluster_idx <- grepl(
    "^Cluster[_ ]?[0-9]+$",
    x,
    ignore.case = TRUE
  )

  if (sum(cluster_idx) > 1) {
    cluster_labels <- x[cluster_idx]

    cluster_numbers <- as.integer(
      sub(
        ".*?([0-9]+)$",
        "\\1",
        cluster_labels,
        perl = TRUE
      )
    )

    x[cluster_idx] <- cluster_labels[
      order(cluster_numbers, cluster_labels)
    ]
  }

  x
}


# -------------------------------------------------------------------------
# Read ADMIXTURE Q file
# -------------------------------------------------------------------------

q <- as.matrix(
  read.table(
    q_path,
    header = FALSE,
    check.names = FALSE
  )
)

storage.mode(q) <- "numeric"


# -------------------------------------------------------------------------
# Read PLINK FAM file
# -------------------------------------------------------------------------

fam <- read.table(
  fam_path,
  header = FALSE,
  stringsAsFactors = FALSE
)


# -------------------------------------------------------------------------
# Validate ADMIXTURE and FAM input
# -------------------------------------------------------------------------

if (is.na(k) || k < 1) {
  stop("K must be an integer greater than or equal to 1")
}

if (ncol(fam) < 2) {
  stop("The FAM file has fewer than two columns")
}

if (nrow(q) != nrow(fam)) {
  stop(
    "Q and FAM row counts differ: ",
    nrow(q),
    " rows in Q and ",
    nrow(fam),
    " rows in FAM"
  )
}

if (ncol(q) != k) {
  stop(
    "The number of Q columns does not equal K: ",
    ncol(q),
    " columns in Q and K = ",
    k
  )
}

if (any(!is.finite(q))) {
  stop("The Q file contains non-numeric or non-finite values")
}

if (any(q < 0)) {
  stop("The Q file contains negative ancestry proportions")
}


# -------------------------------------------------------------------------
# Normalize Q rows exactly to 1
# -------------------------------------------------------------------------

q_row_sums <- rowSums(q)

if (any(!is.finite(q_row_sums))) {
  stop("At least one Q-file row has a non-finite sum")
}

if (any(q_row_sums <= 0)) {
  stop("At least one Q-file row has a sum less than or equal to zero")
}

# This removes minor ADMIXTURE rounding differences such as
# 0.999999 or 1.000001.
q <- q / q_row_sums


# -------------------------------------------------------------------------
# Extract sample IDs
# -------------------------------------------------------------------------

samples <- as.character(fam[[2]])

if (anyDuplicated(samples)) {
  stop("The FAM file contains duplicate sample IDs")
}


# -------------------------------------------------------------------------
# Read keep lists and metadata
# -------------------------------------------------------------------------

cag_keep <- read_keep(cag_keep_path)
ref_keep <- read_keep(ref_keep_path)

cag_meta <- read_metadata(
  cag_meta_path,
  "CaG"
)

ref_meta <- read_metadata(
  ref_meta_path,
  "Reference"
)

meta <- rbind(
  ref_meta,
  cag_meta
)

if (anyDuplicated(meta$sample_id)) {
  stop("At least one sample_id occurs in both metadata files")
}

overlapping_keep_samples <- intersect(
  cag_keep,
  ref_keep
)

if (length(overlapping_keep_samples) > 0) {
  stop(
    length(overlapping_keep_samples),
    " sample IDs occur in both keep lists. ",
    "The first overlapping sample is: ",
    overlapping_keep_samples[1]
  )
}


# -------------------------------------------------------------------------
# Assign cohort and population labels
# -------------------------------------------------------------------------

cohort <- rep(
  NA_character_,
  length(samples)
)

population <- rep(
  NA_character_,
  length(samples)
)

names(cohort) <- samples
names(population) <- samples


# Assign labels using metadata first.
if (nrow(meta) > 0) {
  matched <- match(
    samples,
    meta$sample_id
  )

  has_metadata <- !is.na(matched)

  cohort[has_metadata] <- meta$cohort[
    matched[has_metadata]
  ]

  population[has_metadata] <- meta$population[
    matched[has_metadata]
  ]
}


# Assign cohort using keep lists for samples without metadata.
cohort[
  is.na(cohort) &
    samples %in% ref_keep
] <- "Reference"

cohort[
  is.na(cohort) &
    samples %in% cag_keep
] <- "CaG"


# Retain samples that cannot be assigned to either cohort.
cohort[is.na(cohort)] <- "Unlabelled"


# Assign general population labels when detailed metadata is absent.
population[
  is.na(population) &
    cohort == "Reference"
] <- "Reference"

population[
  is.na(population) &
    cohort == "CaG"
] <- "CaG"

population[is.na(population)] <- "Unlabelled"


# -------------------------------------------------------------------------
# Define cohort ordering
# -------------------------------------------------------------------------

cohort_levels <- c(
  "Reference",
  "CaG",
  "Unlabelled"
)

cohort_levels <- cohort_levels[
  cohort_levels %in% cohort
]


# -------------------------------------------------------------------------
# Define population ordering within each cohort
# -------------------------------------------------------------------------

reference_population_levels <- unique(
  c(
    ref_meta$population,
    if (
      any(
        cohort == "Reference" &
          population == "Reference"
      )
    ) {
      "Reference"
    }
  )
)

reference_population_levels <- reference_population_levels[
  reference_population_levels %in%
    population[cohort == "Reference"]
]


cag_population_levels <- unique(
  c(
    cag_meta$population,
    if (
      any(
        cohort == "CaG" &
          population == "CaG"
      )
    ) {
      "CaG"
    }
  )
)

cag_population_levels <- cag_population_levels[
  cag_population_levels %in%
    population[cohort == "CaG"]
]

cag_population_levels <- order_cluster_levels(
  cag_population_levels
)


unlabelled_population_levels <- unique(
  population[cohort == "Unlabelled"]
)


population_levels <- unique(
  c(
    reference_population_levels,
    cag_population_levels,
    unlabelled_population_levels
  )
)


# -------------------------------------------------------------------------
# Determine sample ordering
# -------------------------------------------------------------------------

metadata_rank <- rep(
  Inf,
  length(samples)
)

if (nrow(meta) > 0) {
  idx <- match(
    samples,
    meta$sample_id
  )

  has_metadata <- !is.na(idx)

  metadata_rank[has_metadata] <- meta$metadata_order[
    idx[has_metadata]
  ]
}


cohort_rank <- match(
  cohort,
  cohort_levels
)

population_rank <- match(
  population,
  population_levels
)


ord <- order(
  cohort_rank,
  population_rank,
  metadata_rank,
  seq_along(samples)
)


# -------------------------------------------------------------------------
# Write plotted sample order
# -------------------------------------------------------------------------

plot_order <- data.frame(
  plot_position = seq_along(ord),
  sample_id = samples[ord],
  cohort = cohort[ord],
  population = population[ord],
  stringsAsFactors = FALSE
)


# Reset sample position inside each cohort/population facet.
facet_group <- interaction(
  plot_order$cohort,
  plot_order$population,
  drop = TRUE
)

plot_order$facet_position <- ave(
  plot_order$plot_position,
  facet_group,
  FUN = seq_along
)


write.table(
  plot_order,
  output_order,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# -------------------------------------------------------------------------
# Prepare ancestry data
# -------------------------------------------------------------------------

q_ordered <- q[
  ord,
  ,
  drop = FALSE
]

k_component_levels <- paste0(
  "K",
  seq_len(k)
)


plot_data <- data.frame(
  sample_id = rep(
    plot_order$sample_id,
    times = k
  ),
  cohort = rep(
    plot_order$cohort,
    times = k
  ),
  population = rep(
    plot_order$population,
    times = k
  ),
  facet_position = rep(
    plot_order$facet_position,
    times = k
  ),
  k_component = factor(
    rep(
      k_component_levels,
      each = nrow(q_ordered)
    ),
    levels = k_component_levels
  ),
  proportion = as.vector(q_ordered),
  stringsAsFactors = FALSE
)


plot_data$cohort <- factor(
  plot_data$cohort,
  levels = cohort_levels
)


# -------------------------------------------------------------------------
# Validate plotting data
# -------------------------------------------------------------------------

if (any(is.na(plot_data$k_component))) {
  stop("At least one ancestry component could not be assigned a K label")
}

if (any(!is.finite(plot_data$proportion))) {
  stop("The plotting data contain non-finite ancestry proportions")
}


# -------------------------------------------------------------------------
# Fixed ADMIXTURE color palette
# -------------------------------------------------------------------------

admixture_palette <- c(
  "#d00000",
  "#ffba08",
  "#cbff8c",
  "#8fe388",
  "#1b998b",
  "#3185fc",
  "#5d2e8c",
  "#46237a",
  "#ff7b9c",
  "#ff9b85"
)

if (k > length(admixture_palette)) {
  stop(
    "K = ",
    k,
    " exceeds the number of available colors (",
    length(admixture_palette),
    ")."
  )
}

colors <- stats::setNames(
  admixture_palette[seq_len(k)],
  k_component_levels
)


# -------------------------------------------------------------------------
# Function to create one cohort panel
# -------------------------------------------------------------------------

make_cohort_plot <- function(
  data,
  cohort_name,
  population_order
) {
  if (nrow(data) == 0) {
    return(NULL)
  }

  present_populations <- population_order[
    population_order %in% data$population
  ]

  if (length(present_populations) == 0) {
    stop(
      "No population levels were found for cohort: ",
      cohort_name
    )
  }

  data$population <- factor(
    as.character(data$population),
    levels = present_populations
  )

  if (any(is.na(data$population))) {
    stop(
      "At least one population in cohort ",
      cohort_name,
      " does not have a corresponding facet level"
    )
  }

  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = facet_position,
      y = proportion,
      fill = k_component
    )
  ) +
    ggplot2::geom_col(
      width = 1.01,
      colour = NA,
      linewidth = 0
    ) +
    ggplot2::facet_grid(
      cols = ggplot2::vars(population),
      scales = "free_x",
      space = "free_x",
      switch = "x"
    ) +
    ggplot2::scale_fill_manual(
      values = colors,
      breaks = names(colors),
      drop = FALSE,

      # If a component ever fails to match the palette, it will be
      # shown in black instead of silently appearing white.
      na.value = "#000000"
    ) +
    ggplot2::scale_y_continuous(
      expand = c(0, 0),
      breaks = c(0, 0.5, 1)
    ) +
    ggplot2::scale_x_continuous(
      expand = c(0, 0)
    ) +
    ggplot2::coord_cartesian(
      ylim = c(0, 1),
      expand = FALSE,
      clip = "on"
    ) +
    ggplot2::labs(
      title = cohort_name,
      x = NULL,
      y = "Ancestry proportion"
    ) +
    ggplot2::theme_bw(
      base_size = 12
    ) +
    ggplot2::theme(
      panel.spacing.x = grid::unit(
        0.12,
        "lines"
      ),

      panel.grid = ggplot2::element_blank(),

      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),

      strip.placement = "outside",

      strip.background = ggplot2::element_rect(
        fill = "grey95",
        colour = "grey50"
      ),

      strip.text.x.bottom = ggplot2::element_text(
        angle = 45,
        hjust = 1,
        vjust = 1
      ),

      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0,
        size = 12
      ),

      legend.position = "none"
    )
}


# -------------------------------------------------------------------------
# Create separate cohort panels
# -------------------------------------------------------------------------

plot_list <- list()


if ("Reference" %in% cohort_levels) {
  reference_data <- plot_data[
    as.character(plot_data$cohort) == "Reference",
    ,
    drop = FALSE
  ]

  plot_list[[length(plot_list) + 1]] <- make_cohort_plot(
    data = reference_data,
    cohort_name = "Reference",
    population_order = reference_population_levels
  )
}


if ("CaG" %in% cohort_levels) {
  cag_data <- plot_data[
    as.character(plot_data$cohort) == "CaG",
    ,
    drop = FALSE
  ]

  plot_list[[length(plot_list) + 1]] <- make_cohort_plot(
    data = cag_data,
    cohort_name = "CaG",
    population_order = cag_population_levels
  )
}


if ("Unlabelled" %in% cohort_levels) {
  unlabelled_data <- plot_data[
    as.character(plot_data$cohort) == "Unlabelled",
    ,
    drop = FALSE
  ]

  plot_list[[length(plot_list) + 1]] <- make_cohort_plot(
    data = unlabelled_data,
    cohort_name = "Unlabelled",
    population_order = unlabelled_population_levels
  )
}


# Remove any empty plots.
plot_list <- Filter(
  Negate(is.null),
  plot_list
)

if (length(plot_list) == 0) {
  stop("No samples were available for plotting")
}


# -------------------------------------------------------------------------
# Combine cohort panels
# -------------------------------------------------------------------------

p <- patchwork::wrap_plots(
  plot_list,
  ncol = 1,
  heights = rep(
    1,
    length(plot_list)
  )
) +
  patchwork::plot_annotation(
    title = paste(
      "ADMIXTURE K =",
      k
    ),
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0.5,
        size = 14
      )
    )
  )


# -------------------------------------------------------------------------
# Determine output dimensions
# -------------------------------------------------------------------------

samples_per_cohort <- table(
  factor(
    cohort,
    levels = cohort_levels
  )
)

maximum_row_samples <- max(
  samples_per_cohort
)

plot_width <- max(
  12,
  min(
    40,
    8 + maximum_row_samples / 350
  )
)

plot_height <- max(
  6,
  3.75 * length(plot_list)
)


# -------------------------------------------------------------------------
# Save PNG
# -------------------------------------------------------------------------

ggplot2::ggsave(
  filename = output_png,
  plot = p,
  width = plot_width,
  height = plot_height,
  units = "in",
  dpi = 300,
  limitsize = FALSE,
  bg = "white"
)