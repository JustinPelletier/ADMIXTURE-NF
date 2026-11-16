#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

required_packages <- c(
  "ggplot2",
  "ggh4x",
  "patchwork"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    "The following R packages are required: ",
    paste(missing_packages, collapse = ", ")
  )
}

if (length(args) != 9) {
  stop(paste(
    "Usage: plot_admixture.R Q FAM K DATASET_METADATA REF_METADATA",
    "DATASET_KEEP REF_KEEP OUTPUT_PNG OUTPUT_ORDER"
  ))
}


# -------------------------------------------------------------------------
# Arguments
# -------------------------------------------------------------------------

q_path <- args[1]
fam_path <- args[2]
k <- as.integer(args[3])

dataset_meta_path <- args[4]
reference_meta_path <- args[5]

dataset_keep_path <- args[6]
reference_keep_path <- args[7]

output_png <- args[8]
output_order <- args[9]

output_pdf <- sub(
  "\\.png$",
  ".pdf",
  output_png,
  ignore.case = TRUE
)


# -------------------------------------------------------------------------
# Functions for reading input files
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


read_metadata <- function(
  path,
  cohort_name,
  require_super_population = FALSE
) {
  empty_metadata <- data.frame(
    sample_id = character(),
    population = character(),
    super_population = character(),
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

  required_columns <- c(
    "sample_id",
    "population"
  )

  if (require_super_population) {
    required_columns <- c(
      required_columns,
      "super_population"
    )
  }

  if (!all(required_columns %in% names(x))) {
    stop(
      path,
      " must contain the columns: ",
      paste(required_columns, collapse = ", ")
    )
  }

  # The DataSet metadata does not require a super_population column.
  if (!"super_population" %in% names(x)) {
    x$super_population <- cohort_name
  }

  x <- x[
    ,
    c(
      "sample_id",
      "population",
      "super_population"
    ),
    drop = FALSE
  ]

  x$sample_id <- trimws(
    as.character(x$sample_id)
  )

  x$population <- trimws(
    as.character(x$population)
  )

  x$super_population <- trimws(
    as.character(x$super_population)
  )

  x <- x[
    nzchar(x$sample_id) &
      nzchar(x$population) &
      nzchar(x$super_population),
    ,
    drop = FALSE
  ]

  if (anyDuplicated(x$sample_id)) {
    stop(
      "Duplicate sample_id values were found in ",
      path
    )
  }

  x$cohort <- cohort_name
  x$metadata_order <- seq_len(nrow(x))

  x
}


# Numerically order Cluster_1, Cluster_2, ..., Cluster_10.
order_cluster_levels <- function(x) {
  is_cluster <- grepl(
    "^Cluster[_ ]?[0-9]+$",
    x,
    ignore.case = TRUE
  )

  cluster_labels <- x[is_cluster]
  other_labels <- x[!is_cluster]

  if (length(cluster_labels) > 0) {
    cluster_numbers <- as.integer(
      sub(
        ".*?([0-9]+)$",
        "\\1",
        cluster_labels,
        perl = TRUE
      )
    )

    cluster_labels <- cluster_labels[
      order(cluster_numbers, cluster_labels)
    ]
  }

  c(
    cluster_labels,
    other_labels
  )
}


# -------------------------------------------------------------------------
# Read ADMIXTURE and PLINK files
# -------------------------------------------------------------------------

q <- as.matrix(
  read.table(
    q_path,
    header = FALSE,
    check.names = FALSE
  )
)

storage.mode(q) <- "numeric"

fam <- read.table(
  fam_path,
  header = FALSE,
  stringsAsFactors = FALSE
)


# -------------------------------------------------------------------------
# Validate ADMIXTURE input
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
    "The Q file contains ",
    ncol(q),
    " columns, but K = ",
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
# Normalize every Q-file row to exactly 1
# -------------------------------------------------------------------------

q_row_sums <- rowSums(q)

if (any(!is.finite(q_row_sums))) {
  stop("At least one Q-file row has a non-finite sum")
}

if (any(q_row_sums <= 0)) {
  stop("At least one Q-file row has a sum <= 0")
}

q <- q / q_row_sums


# -------------------------------------------------------------------------
# Read sample information
# -------------------------------------------------------------------------

samples <- as.character(fam[[2]])

if (anyDuplicated(samples)) {
  stop("The FAM file contains duplicate sample IDs")
}

dataset_keep <- read_keep(dataset_keep_path)
reference_keep <- read_keep(reference_keep_path)

dataset_meta <- read_metadata(
  dataset_meta_path,
  cohort_name = "DataSet",
  require_super_population = FALSE
)

reference_meta <- read_metadata(
  reference_meta_path,
  cohort_name = "Reference",
  require_super_population = TRUE
)

meta <- rbind(
  reference_meta,
  dataset_meta
)

if (anyDuplicated(meta$sample_id)) {
  stop("At least one sample_id occurs in both metadata files")
}

overlapping_keep_samples <- intersect(
  dataset_keep,
  reference_keep
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
# Assign cohort, population and super-population
# -------------------------------------------------------------------------

cohort <- rep(
  NA_character_,
  length(samples)
)

population <- rep(
  NA_character_,
  length(samples)
)

super_population <- rep(
  NA_character_,
  length(samples)
)

metadata_rank <- rep(
  Inf,
  length(samples)
)

names(cohort) <- samples
names(population) <- samples
names(super_population) <- samples


# Assign metadata labels first.
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

  super_population[has_metadata] <- meta$super_population[
    matched[has_metadata]
  ]

  metadata_rank[has_metadata] <- meta$metadata_order[
    matched[has_metadata]
  ]
}


# Assign remaining samples using the keep lists.
cohort[
  is.na(cohort) &
    samples %in% reference_keep
] <- "Reference"

cohort[
  is.na(cohort) &
    samples %in% dataset_keep
] <- "DataSet"

cohort[is.na(cohort)] <- "Unlabelled"


# Assign fallback population labels.
population[
  is.na(population) &
    cohort == "Reference"
] <- "Reference"

population[
  is.na(population) &
    cohort == "DataSet"
] <- "DataSet"

population[
  is.na(population)
] <- "Unlabelled"


# Assign fallback super-population labels.
super_population[
  is.na(super_population) &
    cohort == "Reference"
] <- "Reference"

super_population[
  is.na(super_population) &
    cohort == "DataSet"
] <- "DataSet"

super_population[
  is.na(super_population)
] <- "Unlabelled"


# -------------------------------------------------------------------------
# Define cohort order
# -------------------------------------------------------------------------

cohort_levels <- c(
  "Reference",
  "DataSet",
  "Unlabelled"
)

cohort_levels <- cohort_levels[
  cohort_levels %in% cohort
]


# -------------------------------------------------------------------------
# Define reference super-population order
# -------------------------------------------------------------------------

preferred_super_population_order <- c(
  "EUR",
  "AFR",
  "MID",
  "CSA",
  "EAS",
  "AMR",
  "OCE"
)

observed_reference_super_populations <- unique(
  super_population[
    cohort == "Reference"
  ]
)

reference_super_population_levels <- c(
  preferred_super_population_order[
    preferred_super_population_order %in%
      observed_reference_super_populations
  ],
  observed_reference_super_populations[
    !observed_reference_super_populations %in%
      preferred_super_population_order
  ]
)

reference_super_population_levels <- unique(
  reference_super_population_levels
)


# -------------------------------------------------------------------------
# Define reference population order within super-populations
# -------------------------------------------------------------------------

reference_population_levels <- unlist(
  lapply(
    reference_super_population_levels,
    function(super_pop) {
      unique(
        reference_meta$population[
          reference_meta$super_population == super_pop
        ]
      )
    }
  ),
  use.names = FALSE
)

reference_population_levels <- unique(
  c(
    reference_population_levels,
    population[
      cohort == "Reference" &
        !population %in% reference_population_levels
    ]
  )
)


# -------------------------------------------------------------------------
# Define DataSet population order
# -------------------------------------------------------------------------

dataset_population_levels <- unique(
  c(
    dataset_meta$population,
    population[
      cohort == "DataSet" &
        !population %in% dataset_meta$population
    ]
  )
)

dataset_population_levels <- order_cluster_levels(
  dataset_population_levels
)


unlabelled_population_levels <- unique(
  population[
    cohort == "Unlabelled"
  ]
)


# -------------------------------------------------------------------------
# Order samples
# -------------------------------------------------------------------------

cohort_rank <- match(
  cohort,
  cohort_levels
)

reference_super_rank <- match(
  super_population,
  reference_super_population_levels
)

reference_population_rank <- match(
  population,
  reference_population_levels
)

dataset_population_rank <- match(
  population,
  dataset_population_levels
)

unlabelled_population_rank <- match(
  population,
  unlabelled_population_levels
)


population_rank <- rep(
  Inf,
  length(samples)
)

super_population_rank <- rep(
  Inf,
  length(samples)
)


reference_idx <- cohort == "Reference"
dataset_idx <- cohort == "DataSet"
unlabelled_idx <- cohort == "Unlabelled"

super_population_rank[reference_idx] <-
  reference_super_rank[reference_idx]

population_rank[reference_idx] <-
  reference_population_rank[reference_idx]

super_population_rank[dataset_idx] <- 1

population_rank[dataset_idx] <-
  dataset_population_rank[dataset_idx]

super_population_rank[unlabelled_idx] <- 1

population_rank[unlabelled_idx] <-
  unlabelled_population_rank[unlabelled_idx]


ord <- order(
  cohort_rank,
  super_population_rank,
  population_rank,
  metadata_rank,
  seq_along(samples)
)


# -------------------------------------------------------------------------
# Create and write the plotted sample order
# -------------------------------------------------------------------------

plot_order <- data.frame(
  plot_position = seq_along(ord),
  sample_id = samples[ord],
  cohort = cohort[ord],
  super_population = super_population[ord],
  population = population[ord],
  stringsAsFactors = FALSE
)


facet_group <- interaction(
  plot_order$cohort,
  plot_order$super_population,
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
# Prepare long-format ancestry data
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
  super_population = rep(
    plot_order$super_population,
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
# Fixed ADMIXTURE palette
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
    ")"
  )
}

colors <- stats::setNames(
  admixture_palette[seq_len(k)],
  k_component_levels
)


# -------------------------------------------------------------------------
# Shared plot layers
# -------------------------------------------------------------------------

common_plot_layers <- list(
  ggplot2::geom_col(
    width = 1.01,
    colour = NA,
    linewidth = 0
  ),

  ggplot2::scale_fill_manual(
    values = colors,
    breaks = names(colors),
    drop = FALSE,
    na.value = "#000000"
  ),

  ggplot2::scale_y_continuous(
    expand = c(0, 0),
    breaks = c(0, 0.5, 1)
  ),

  ggplot2::scale_x_continuous(
    expand = c(0, 0)
  ),

  ggplot2::coord_cartesian(
    ylim = c(0, 1),
    expand = FALSE,
    clip = "on"
  ),

  ggplot2::theme_bw(
    base_size = 12
  ),

  ggplot2::theme(
    panel.spacing.x = grid::unit(
      0.08,
      "lines"
    ),

    panel.grid = ggplot2::element_blank(),

    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),

    strip.placement = "outside",

    plot.title = ggplot2::element_text(
      face = "bold",
      hjust = 0,
      size = 12
    ),

    legend.position = "none"
  )
)


# -------------------------------------------------------------------------
# Reference plot with nested super-population strips
# -------------------------------------------------------------------------

make_reference_plot <- function(data) {
  if (nrow(data) == 0) {
    return(NULL)
  }

  present_super_populations <-
    reference_super_population_levels[
      reference_super_population_levels %in%
        data$super_population
    ]

  present_populations <-
    reference_population_levels[
      reference_population_levels %in%
        data$population
    ]

  data$super_population <- factor(
    as.character(data$super_population),
    levels = present_super_populations
  )

  data$population <- factor(
    as.character(data$population),
    levels = present_populations
  )

  if (
    any(is.na(data$super_population)) ||
      any(is.na(data$population))
  ) {
    stop(
      "At least one Reference sample could not be assigned ",
      "to a super-population or population facet"
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
    common_plot_layers +

    ggh4x::facet_nested(
      cols = ggplot2::vars(
        super_population,
        population
      ),
      scales = "free_x",
      space = "free_x",
      switch = "x",

      nest_line = ggplot2::element_line(
        colour = "grey45",
        linewidth = 0.35
      ),

      strip = ggh4x::strip_nested(
        bleed = FALSE,

        # Variable strip sizes allow the super-population strip to
        # remain compact while population labels remain angled.
        size = "variable",

        by_layer_x = TRUE,

        text_x = list(
          # Outermost strip: super-population.
          ggplot2::element_text(
            angle = 0,
            face = "bold",
            size = 8,
            hjust = 0.5,
            vjust = 0.5,
            margin = ggplot2::margin(
              t = 1.5,
              r = 2,
              b = 1.5,
              l = 2,
              unit = "pt"
            )
          ),

          # Inner strip: population.
          ggplot2::element_text(
            angle = 90,
            face = "plain",
            size = 7,
            hjust = 1,
            vjust = 1,
            margin = ggplot2::margin(
              t = 2,
              r = 2,
              b = 2,
              l = 2,
              unit = "pt"
            )
          )
        ),

        background_x = list(
          # Super-population background.
          ggplot2::element_rect(
            fill = "grey85",
            colour = "grey40",
            linewidth = 0.35
          ),

          # Population background.
          ggplot2::element_rect(
            fill = "grey95",
            colour = "grey60",
            linewidth = 0.25
          )
        )
      )
    ) +

    ggplot2::labs(
      title = "Reference",
      x = NULL,
      y = "Ancestry proportion"
    )
}


# -------------------------------------------------------------------------
# DataSet or unlabelled plot
# -------------------------------------------------------------------------

make_standard_plot <- function(
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

  data$population <- factor(
    as.character(data$population),
    levels = present_populations
  )

  if (any(is.na(data$population))) {
    stop(
      "At least one population in ",
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
    common_plot_layers +

    ggplot2::facet_grid(
      cols = ggplot2::vars(population),
      scales = "free_x",
      space = "free_x",
      switch = "x"
    ) +

    ggplot2::labs(
      title = cohort_name,
      x = NULL,
      y = "Ancestry proportion"
    ) +

    ggplot2::theme(
      strip.background = ggplot2::element_rect(
        fill = "grey95",
        colour = "grey50"
      ),

      strip.text.x.bottom = ggplot2::element_text(
        angle = 90,
        hjust = 1,
        vjust = 1,
        size = 8
      )
    )
}


# -------------------------------------------------------------------------
# Create plot panels
# -------------------------------------------------------------------------

plot_list <- list()


if ("Reference" %in% cohort_levels) {
  reference_data <- plot_data[
    as.character(plot_data$cohort) == "Reference",
    ,
    drop = FALSE
  ]

  plot_list[[length(plot_list) + 1]] <-
    make_reference_plot(reference_data)
}


if ("DataSet" %in% cohort_levels) {
  dataset_data <- plot_data[
    as.character(plot_data$cohort) == "DataSet",
    ,
    drop = FALSE
  ]

  plot_list[[length(plot_list) + 1]] <-
    make_standard_plot(
      data = dataset_data,
      cohort_name = "DataSet",
      population_order = dataset_population_levels
    )
}


if ("Unlabelled" %in% cohort_levels) {
  unlabelled_data <- plot_data[
    as.character(plot_data$cohort) == "Unlabelled",
    ,
    drop = FALSE
  ]

  plot_list[[length(plot_list) + 1]] <-
    make_standard_plot(
      data = unlabelled_data,
      cohort_name = "Unlabelled",
      population_order = unlabelled_population_levels
    )
}


plot_list <- Filter(
  Negate(is.null),
  plot_list
)

if (length(plot_list) == 0) {
  stop("No plot panels could be generated")
}


# -------------------------------------------------------------------------
# Combine panels
# -------------------------------------------------------------------------

combined_plot <- patchwork::wrap_plots(
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
  4.0 * length(plot_list)
)


# -------------------------------------------------------------------------
# Save PNG
# -------------------------------------------------------------------------

ggplot2::ggsave(
  filename = output_png,
  plot = combined_plot,
  width = plot_width,
  height = plot_height,
  units = "in",
  dpi = 600,
  limitsize = FALSE,
  bg = "white"
)

ggplot2::ggsave(
  filename = output_pdf,
  plot = combined_plot,
  width = plot_width,
  height = plot_height,
  units = "in",
  device = "pdf",
  limitsize = FALSE,
  bg = "white"
)
