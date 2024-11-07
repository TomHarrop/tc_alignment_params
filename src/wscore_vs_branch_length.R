#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)


this_marker_type <- "NUC"
stats_file <- "data/03_extractions/captus-assembly_extract.stats.tsv"
branch_length_file <- "output/033_iqtree-stats/NUC.NT.muscle_align.gap0.2.cov0.9.wscore0.2/branch_lengths.csv"

branch_lengths <- fread(branch_length_file)

stats <- fread(stats_file)

median_wscores <- stats[
    sample_name %in% branch_lengths$leaf_name
][
    hit == 0 & marker_type == this_marker_type,
    .(median_wscore = median(wscore, na.rm = TRUE)),
    by = sample_name
]

pd <- merge(
    branch_lengths,
    median_wscores,
    by.x = "leaf_name", by.y = "sample_name"
)

test1 <- lm(
    terminal_branch_length ~ 1 + median_wscore,
    data = pd
)
summary(test1)
cor.test(
  pd$median_wscore,
  pd$terminal_branch_length,
  method = "spearman"
)



dummy_df <- pd[, .(median_wscore)]
pred <- cbind(
    dummy_df,
    predict(test1, newdata = dummy_df, interval = "confidence")
)

ggplot(pd, aes(y = terminal_branch_length, x = median_wscore)) +
    geom_point() +
    geom_line(
        mapping = aes(y = fit),
        data = pred,
        colour = "red"
    ) +
    labs(
        y = "Terminal branch length",
        x = "Median wscore"
    ) +
    theme_minimal()
