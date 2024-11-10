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
    .(
        median_wscore = median(wscore, na.rm = TRUE),
        median_pct_recovered = median(pct_recovered, na.rm = TRUE)
    ),
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
pred1 <- cbind(
    pd[, .(median_wscore)],
    predict(test1, newdata = pd[, .(median_wscore)], interval = "confidence")
)

test2 <- lm(
    terminal_branch_length ~ 1 + median_pct_recovered,
    data = pd
)
summary(test2)
cor.test(
    pd$median_pct_recovered,
    pd$terminal_branch_length,
    method = "spearman"
)
pred2 <- cbind(
    pd[, .(median_pct_recovered)],
    predict(
        test2,
        newdata = pd[, .(median_pct_recovered)],
        interval = "confidence"
    )
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





ggplot(pd, aes(y = terminal_branch_length, x = median_pct_recovered)) +
    geom_point() +
    geom_line(
        mapping = aes(y = fit),
        data = pred2,
        colour = "red"
    ) +
    labs(
        y = "Terminal branch length",
        x = "Median % recovered"
    ) +
    theme_minimal()
