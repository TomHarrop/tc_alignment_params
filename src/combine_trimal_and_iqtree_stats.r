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
    wscore <- sub("wscore", "", paste(split_string[8:9], collapse = "."))
    return(
        list(
            marker = split_string[1],
            marker_format = split_string[2],
            align_method = split_string[3],
            gap = as.numeric(gap),
            cov = as.numeric(cov),
            wscore = as.numeric(wscore)
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

default_param_string <- "NUC.NT.mafft_auto.gap0.9.cov0.4.wscore0.0"

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
all_metrics_with_params[naive_tree_score == best_results$naive_tree_score]

# plots
naive_tree_score_plot <- ggplot(
    all_metrics_with_params,
    aes(
        x = cov,
        y = naive_tree_score,
        colour = as.factor(gap),
        group = as.factor(gap)
    )
) +
    scale_colour_viridis_d(guide = guide_legend(title = "--clipkit_gaps")) +
    facet_grid(align_method ~ as.factor(wscore)) +
    xlab("--min_coverage") +
    geom_point(shape = 16, size = 1, alpha = 0.8) +
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
    ylab("Naïve tree score") +
    geom_hline(
        yintercept = default_results$naive_tree_score,
        colour = "black",
        linetype = "dashed",
        alpha = 0.5
    )

ggsave("test/naive_tree_score.pdf",
    naive_tree_score_plot,
    width = 210,
    height = 297,
    units = "mm",
    device = cairo_pdf
)


normalised_branch_length_plot <- ggplot(
    all_metrics_with_params,
    aes(
        x = cov,
        y = normalised_total_tree_length,
        colour = as.factor(gap),
        group = as.factor(gap)
    )
) +
    scale_colour_viridis_d(guide = guide_legend(title = "--clipkit_gaps")) +
    facet_grid(align_method ~ as.factor(wscore)) +
    xlab("--min_coverage") +
    geom_point(shape = 16, size = 1, alpha = 0.8) +
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
    ylab("Normalised total tree length") +
    geom_hline(
        yintercept = default_results$normalised_total_tree_length,
        colour = "black",
        linetype = "dashed",
        alpha = 0.5
    )

ggsave("test/normalised_branch_length.pdf",
    normalised_branch_length_plot,
    width = 210,
    height = 297,
    units = "mm",
    device = cairo_pdf
)


normalised_informative_site_plot <- ggplot(
    all_metrics_with_params,
    aes(
        x = cov,
        y = normalised_sum_of_informative_sites,
        colour = as.factor(gap),
        group = as.factor(gap)
    )
) +
    scale_colour_viridis_d(guide = guide_legend(title = "--clipkit_gaps")) +
    facet_grid(align_method ~ as.factor(wscore)) +
    xlab("--min_coverage") +
    geom_point(shape = 16, size = 1, alpha = 0.8) +
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
    ylab("Normalised sum of informative sites") +
    geom_hline(
        yintercept = default_results$normalised_sum_of_informative_sites,
        colour = "black",
        linetype = "dashed",
        alpha = 0.5
    )

ggsave("test/normalised_sum_of_informative_sites.pdf",
    normalised_informative_site_plot,
    width = 210,
    height = 297,
    units = "mm",
    device = cairo_pdf
)

gap_score_plot <- ggplot(
    all_metrics_with_params,
    aes(
        x = cov,
        y = MedianGapScore_median,
        colour = as.factor(gap),
        group = as.factor(gap)
    )
) +
    scale_colour_viridis_d(guide = guide_legend(title = "--clipkit_gaps")) +
    facet_grid(align_method ~ as.factor(wscore)) +
    xlab("--min_coverage") +
    geom_point(shape = 16, size = 1, alpha = 0.8) +
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
    ylab("Median of median gap score") +
    geom_hline(
        yintercept = default_results$MedianGapScore_median,
        colour = "black",
        linetype = "dashed",
        alpha = 0.5
    )

ggsave("test/MedianGapScore_median.pdf",
    gap_score_plot,
    width = 210,
    height = 297,
    units = "mm",
    device = cairo_pdf
)
