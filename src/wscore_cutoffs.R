#!/usr/bin/env Rscript

if (exists("snakemake")) {
    log <- file(snakemake@log[["log"]], open = "wt")
    sink(log, type = "message")
    sink(log, append = TRUE, type = "output")

    stats_file <- snakemake@params[["stats_file"]]
    discarded_samples_file <- snakemake@output[["discarded_samples"]]
    discarded_loci_file <- snakemake@output[["discarded_loci"]]

    sample_wscore_cutoff <- snakemake@wildcards[["sample_wscore_cutoff"]]
    message("locus_wscore_cutoff not implemented, setting to sample_wscore_cutoff")
    locus_wscore_cutoff <- sample_wscore_cutoff
    this_marker_type <- snakemake@wildcards[["marker"]]
} else {
    locus_wscore_cutoff <- 0.25
    sample_wscore_cutoff <- 0.25
    this_marker_type <- "NUC"
    stats_file <- "data/03_extractions/captus-assembly_extract.stats.tsv"
    discarded_samples_file <- "test/discarded_samples.csv"
    discarded_loci_file <- "test/discarded_loci.csv"
}

library(data.table)

stats <- fread(stats_file)

kept_samples <- stats[
    hit == 0 & marker_type == this_marker_type,
    median(wscore, na.rm = TRUE),
    by = sample_name
][V1 >= sample_wscore_cutoff, unique(sample_name)]

discarded_samples <- stats[
    !sample_name %in% kept_samples,
    unique(sample_name)
]



kept_loci <- stats[
    hit == 0 &
        sample_name %in% kept_samples &
        marker_type == this_marker_type,
    median(wscore, na.rm = TRUE),
    by = locus
][V1 >= locus_wscore_cutoff, unique(locus)]
discarded_loci <- stats[
    marker_type == this_marker_type & !locus %in% kept_loci,
    unique(locus)
]


fwrite(list(discarded_samples), discarded_samples_file)
fwrite(list(discarded_loci), discarded_loci_file)

sessionInfo()
