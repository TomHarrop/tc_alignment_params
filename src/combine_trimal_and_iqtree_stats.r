#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)

FileNameToParamString <- function(x) {
    return(
        unlist(
            strsplit(x, "/", fixed = TRUE)
        )[[3]]
    )
}

ParseParamString <- function(x) {
    split_string <- unlist(strsplit(x, ".", fixed = TRUE))
    gap <- sub("gap", "", paste(split_string[4:5], collapse = "."))
    cov <- sub("cov", "", paste(split_string[6:7], collapse = "."))
    return(
        list(
            marker = split_string[1],
            marker_format = split_string[2],
            align_method = split_string[3],
            gap = as.numeric(gap),
            cov = as.numeric(cov)
        )
    )
}


ReadStatsFiles <- function(path, pattern) {
    stats_files <- list.files(
        path,
        pattern = pattern,
        full.names = TRUE,
        recursive = TRUE,
    )

    names(stats_files) <- sapply(
        stats_files, FileNameToParamString
    )

    my_stats <- rbindlist(
        lapply(stats_files, fread),
        idcol = "param_string"
    )

    return(my_stats)
}

# globals
bin_breaks <- seq(0, 1, by = 0.1)

# TODO! it's actually NUC.NT.mafft_auto.gap0.9.cov0.4 but we haven't run that
# yet.
default_param_string <- "NUC.NT.mafft_auto.gap0.8.cov0.4"


# read the trimal stats
trimal_stats <- ReadStatsFiles(
    "output/023_trimal-stats",
    "stats.csv"
)

# read the iqtree stats
site_info_stats <- ReadStatsFiles(
    "output/033_iqtree-stats",
    "site_info.csv"
)
iqtree_metrics <- ReadStatsFiles(
    "output/033_iqtree-stats",
    "metrics.csv"
)

# combine the per-site statistics
per_locus_stats <- merge(
    trimal_stats,
    site_info_stats,
    by.x = c("param_string", "locus"),
    by.y = c("param_string", "Name"),
    all = TRUE
)

per_locus_stats[locus == "4471" & param_string == "NUC.NT.muscle_super5.gap0.8.cov0.0", ]

# Summarise the per-site statistics. There are lots more, check the individual
# stats files.
summarised_alignment_stats <- per_locus_stats[, .(
    MedianGapScore_median = median(MedianGapScore),
    sum_sites = sum(Site, na.rm = TRUE),
    sum_Infor = sum(Infor, na.rm = TRUE)
), by = .(param_string)]

# merge the metrics
all_metrics <- merge(
    summarised_alignment_stats,
    dcast(iqtree_metrics, param_string ~ Metric, value.var = "Value"),
    by = "param_string",
    all = TRUE
)

# parse the parameters
parsed_param_strings <- all_metrics[,
    ParseParamString(param_string),
    by = param_string
]

all_metrics_with_params <- merge(
    all_metrics,
    parsed_param_strings,
    by = "param_string",
    all = TRUE
)

# some sort of combined metric
all_metrics_with_params[
    ,
    normalised_total_tree_length := `Total tree length (sum of branch lengths)` / max(`Total tree length (sum of branch lengths)`, na.rm = TRUE)
]
all_metrics_with_params[
    ,
    normalised_sum_of_informative_sites := sum_Infor / max(sum_Infor, na.rm = TRUE)
]

all_metrics_with_params[
    ,
    naive_tree_score := ((2 * MedianGapScore_median) + normalised_sum_of_informative_sites - (2 * normalised_total_tree_length)) / 5
]


# the default results
default_results <- all_metrics_with_params[param_string == default_param_string]
best_results <- all_metrics_with_params[which.max(naive_tree_score)]

# plots
gp <- ggplot(
    all_metrics_with_params,
    aes(
        x = MedianGapScore_median,
        y = normalised_sum_of_informative_sites,
        colour = gap,
        shape = cov
    )
) +
    scale_colour_viridis_b(breaks = bin_breaks) +
    scale_shape_binned() +
    facet_wrap(~align_method) +
    geom_point() +
    geom_point(
        data = default_results,
        colour = "black",
        shape = 1,
        size = 5,
    ) +
    geom_point(
        data = best_results,
        colour = "red",
        shape = 1,
        size = 5,
    ) +
    geom_vline(
        xintercept = default_results$MedianGapScore_median,
        colour = "black",
        linetype = "dashed",
        alpha = 0.5
    ) +
    geom_hline(
        yintercept = default_results$normalised_sum_of_informative_sites,
        colour = "black",
        linetype = "dashed",
        alpha = 0.5
    )

ggsave("test/total_informative_sites.pdf",
    gp,
    width = 10,
    height = 7.5,
    units = "in",
    device = cairo_pdf
)

gp2 <- ggplot(
    all_metrics_with_params,
    aes(
        x = MedianGapScore_median,
        y = normalised_total_tree_length,
        colour = gap,
        shape = cov
    )
) +
    scale_colour_viridis_b(breaks = bin_breaks) +
    scale_shape_binned() +
    facet_wrap(~align_method) +
    geom_point(size = 4) +
    geom_point(
        data = default_results,
        colour = "black",
        shape = 1,
        size = 5,
    ) +
    geom_point(
        data = best_results,
        colour = "red",
        shape = 1,
        size = 5,
    ) +
    geom_vline(
        xintercept = default_results$MedianGapScore_median,
        colour = "black",
        linetype = "dashed",
        alpha = 0.5
    ) +
    geom_hline(
        yintercept = default_results$normalised_total_tree_length,
        colour = "black",
        linetype = "dashed",
        alpha = 0.5
    )

ggsave("test/total_tree_length.pdf",
    gp2,
    width = 10,
    height = 7.5,
    units = "in",
    device = cairo_pdf
)

gp3 <- ggplot(
    all_metrics_with_params,
    aes(
        x = normalised_sum_of_informative_sites,
        y = normalised_total_tree_length,
        colour = gap,
        shape = cov
    )
) +
    scale_colour_viridis_b(breaks = bin_breaks) +
    scale_shape_binned() +
    facet_wrap(~align_method) +
    geom_point(size = 4) +
    geom_point(
        data = default_results,
        colour = "black",
        shape = 1,
        size = 6,
    ) +
    geom_point(
        data = best_results,
        colour = "red",
        shape = 1,
        size = 6,
    ) +
    geom_vline(
        xintercept = default_results$normalised_sum_of_informative_sites,
        colour = "black",
        linetype = "dashed",
        alpha = 0.5
    ) +
    geom_hline(
        yintercept = default_results$normalised_total_tree_length,
        colour = "black",
        linetype = "dashed",
        alpha = 0.5
    )


gp3 <- ggplot(
    all_metrics_with_params,
    aes(
        x = cov,
        y = naive_tree_score,
        colour = as.factor(gap),
        group = as.factor(gap)
    )
) +
    scale_colour_viridis_d() +
    facet_wrap(~align_method) +
    geom_point(shape = 16,size = 1, alpha = 0.8) +
    geom_path() +
    geom_point(
        data = default_results,
        colour = "black",
        shape = 1,
        size = 6,
    ) +
    geom_point(
        data = best_results,
        colour = "red",
        shape = 1,
        size = 6,
    ) +
    geom_hline(
        yintercept = default_results$naive_tree_score,
        colour = "black",
        linetype = "dashed",
        alpha = 0.5
    )

ggsave("test/naive_metric.pdf",
    gp3,
    width = 10,
    height = 7.5,
    units = "in",
    device = cairo_pdf
)
