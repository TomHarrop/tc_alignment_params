#!/usr/bin/env python3

from functools import cache
from pathlib import Path
import numpy as np
import random
import shutil
import tarfile

#############
# FUNCTIONS #
#############


@cache
def get_module_snakefile(module, tag):
    return github(
        "tomharrop/smk-modules",
        path=f"modules/{module}/Snakefile",
        tag=tag,
    )


# just get the alphabetically last directory
def find_alignment_directory(wildcards, input):
    captus_tarfile = Path(input.tarfile)
    if not captus_tarfile.exists():
        return [""]
    with tarfile.open(captus_tarfile, "r") as tar:
        fasta_files = [
            Path(x.name)
            for x in tar.getmembers()
            if x.isfile() and x.name.endswith(".fna")
        ]
    parent_dirs = sorted(set(x.parent for x in fasta_files))
    return parent_dirs[-1]


def ruleio_expand(ruleio):
    return expand(
        [str(x) for x in ruleio],
        align_method=align_methods,
        clipkit_gap=clipkit_gaps,
        min_coverage=min_coverages,
        marker=["NUC"],
        marker_format=["NT"],
        sample_wscore_cutoff=sample_wscore_cutoffs,
    )


###########
# GLOBALS #
###########

# Paths
data = Path("data")
outdir = Path("output")
logdir = Path(outdir, "logs")
benchdir = Path(logdir, "benchmarks")

captus_extractions_directory = Path(data, "03_extractions")

# Params
align_methods = [
    "mafft_auto",
    "mafft_genafpair",
    "mafft_localpair",
    "mafft_globalpair",
    "mafft_retree1",
    "mafft_retree2",
    "muscle_align",
    "muscle_super5",
]
clipkit_gaps = [round(float(x), 1) for x in np.linspace(0, 1, 10)]
min_coverages = [round(float(x), 1) for x in np.linspace(0, 1, 10)]
markers = ["NUC", "PTD", "MIT", "DNA", "CLR", "ALL"]
formats = ["AA", "NT", "GE", "GF", "MA", "MF", "ALL"]
# Note, even at wscore >= 0.0, some samples get discarded. This happens if no
# markes of the correct type were assembled. It's correct to remove these
# samples.
sample_wscore_cutoffs = [
    round(float(x), 1) for x in np.linspace(0, 0.6, num=4)
]

param_string = (
    "{marker}."
    "{marker_format}."
    "{align_method}."
    "gap{clipkit_gap}."
    "cov{min_coverage}."
    "wscore{sample_wscore_cutoff}"
)

# how many samples do we want? it's quite slow.
n_samples = 30

# containers
biopython = "docker://quay.io/biocontainers/biopython:1.81"
captus = "docker://quay.io/biocontainers/captus:1.0.1--pyhdfd78af_2"
ete3 = "docker://quay.io/biocontainers/ete3:3.1.1--py36_0"
pandas = "docker://quay.io/biocontainers/pandas:2.2.1"
r = "docker://ghcr.io/tomharrop/r-containers:r2u_24.04_cv1"
trimal = "docker://quay.io/biocontainers/trimal:1.5.0--h4ac6f70_0"

# modules
module_tag = "0.7.1"

iqtree_snakefile = get_module_snakefile("iqtree", module_tag)
trimal_snakefile = get_module_snakefile("trimal", module_tag)


#########
# RULES #
#########


wildcard_constraints:
    align_method="|".join(align_methods),
    clipkit_gap="|".join(map(str, clipkit_gaps)),
    min_coverage="|".join(map(str, min_coverages)),
    marker="|".join(markers),
    marker_format="|".join(formats),
    sample_wscore_cutoff="|".join(map(str, sample_wscore_cutoffs)),


rule get_tip_length:
    input:
        treefile=Path(
            outdir,
            "030_iqtree",
            param_string,
            "tree.treefile",
        ),
    output:
        branch_length_csv=Path(
            outdir,
            "033_iqtree-stats",
            param_string,
            "branch_lengths.csv",
        ),
    log:
        Path(
            logdir,
            "get_tip_length",
            param_string + ".log",
        ),
    resources:
        time=lambda wildcards, attempt: 10 * attempt,
    container:
        ete3
    shell:
        "python3 src/get_tip_length.py "
        "{input.treefile} "
        "> {output.branch_length_csv} "
        "2> {log}"


rule parse_iqtree_file:
    input:
        iqtree_file=Path(
            outdir,
            "030_iqtree",
            param_string,
            "tree.iqtree",
        ),
    output:
        summary_metrics=Path(
            outdir,
            "033_iqtree-stats",
            param_string,
            "metrics.csv",
        ),
        site_data=Path(
            outdir,
            "033_iqtree-stats",
            param_string,
            "site_info.csv",
        ),
    log:
        Path(
            logdir,
            "parse_iqtree_file",
            param_string + ".log",
        ),
    resources:
        time=lambda wildcards, attempt: 10 * attempt,
    container:
        pandas
    script:
        "src/parse_iqtree_file.py"


module iqtree:
    snakefile:
        iqtree_snakefile
    config:
        {
            "alignment_directory": Path(
                outdir,
                "027_iqtree-input",
                param_string,
            ),
            "outdir": Path(
                outdir,
                "030_iqtree",
                param_string,
            ),
        }


use rule * from iqtree as iqtree_*


use rule iqtree from iqtree as iqtree_iqtree with:
    resources:
        mem_mb=lambda wildcards, attempt: int(32e3 * attempt),
        time=lambda wildcards, attempt: int(180 * attempt),
    threads: 48


rule setup_iqtree_input:
    input:
        tarfile=Path(outdir, "025_trimal-processed", param_string, "kept.tar"),
    output:
        outdir=temp(
            directory(
                Path(
                    outdir,
                    "027_iqtree-input",
                    param_string,
                )
            )
        ),
    container:
        trimal
    shadow:
        "minimal"
    shell:
        "tmpdir=$(mktemp -d) && "
        "tar -xf {input.tarfile} -C ${{tmpdir}} && "
        "mv ${{tmpdir}} {output.outdir}"


rule process_trimal_files:
    input:
        tarfile=Path(outdir, "020_trimal", param_string, "trimal.tar"),
    output:
        kept=Path(outdir, "025_trimal-processed", param_string, "kept.tar"),
        discarded=Path(
            outdir, "025_trimal-processed", param_string, "discarded.tar"
        ),
        empty=Path(outdir, "025_trimal-processed", param_string, "empty.tar"),
    log:
        Path(
            logdir,
            "process_trimal_files",
            param_string + ".log",
        ),
    benchmark:
        Path(
            benchdir,
            "process_trimal_files",
            param_string + ".txt",
        )
    threads: 12
    resources:
        time=lambda wildcards, attempt: 10 * attempt,
    container:
        biopython
    script:
        "src/process_trimal_files.py"


rule parse_trimal_stats:
    input:
        stats_tarfile=Path(
            outdir, "020_trimal", param_string, "trimal.stats.tar"
        ),
    output:
        trimal_stats=Path(
            outdir, "023_trimal-stats", param_string, "stats.csv"
        ),
    log:
        Path(
            logdir,
            "parse_trimal_stats",
            param_string + ".log",
        ),
    resources:
        time=lambda wildcards, attempt: 10 * attempt,
    container:
        pandas
    script:
        "src/parse_trimal_stats.py"


rule trimal:
    input:
        tarfile=Path(outdir, "010_captus-align", param_string + ".tar"),
    output:
        tarfile=Path(outdir, "020_trimal", param_string, "trimal.tar"),
        stats=Path(outdir, "020_trimal", param_string, "trimal.stats.tar"),
        htmlout=Path(outdir, "020_trimal", param_string, "trimal.htmlout.tar"),
    params:
        alignment_directory=find_alignment_directory,
    log:
        Path(
            logdir,
            "trimal",
            param_string + ".log",
        ),
    threads: 1
    resources:
        time=lambda wildcards, attempt: 10 * attempt,
    container:
        trimal
    shadow:
        "minimal"
    shell:
        "tar -xf {input.tarfile} "
        "./{params.alignment_directory} && "
        "mkdir trimmed stats htmlout && "
        "for input in {params.alignment_directory}/*.fna ; do "
        "output=$(basename ${{input}}) ; "
        "trimal "
        "-in ${{input}} "
        "-out trimmed/${{output}} "
        "-automated1 "
        "-htmlout htmlout/${{output}}.html "
        "-sgc "
        "-sgt "
        "-sident "
        "-soverlap "
        "-ssc "
        "-sst "
        "> stats/${{output}}.stats "
        "2>> {log} ; "
        "done && "
        "tar -cf {output.tarfile} --directory trimmed . && "
        "tar -cf {output.stats} --directory stats . && "
        "tar -cf {output.htmlout} --directory htmlout ."


rule captus_align:
    input:
        extraction_dir=Path(outdir, "005_selected-samples"),
        discarded_samples=Path(
            outdir,
            "007_wscore-filtering",
            "{marker}.wscore{sample_wscore_cutoff}",
            "discarded_samples.txt",
        ),
    output:
        tarfile=Path(outdir, "010_captus-align", param_string + ".tar"),
    log:
        Path(
            logdir,
            "align",
            param_string + ".log",
        ),
    benchmark:
        Path(
            benchdir,
            "align",
            param_string + ".txt",
        )
    threads: lambda wildcards, attempt: 32 * attempt
    resources:
        time=lambda wildcards, attempt: 120 * attempt,
        mem_mb=lambda wildcards, attempt: int(16e3 * attempt),
    container:
        captus
    shell:
        "tmp_indir=$(mktemp -d) && tmp_outdir=$(mktemp -d) ; "
        "cp -rL {input.extraction_dir} ${{tmp_indir}}/align_input ; "
        "while read sample; do "
        "if [ -d ${{tmp_indir}}/align_input/${{sample}}__captus-ext ]; then "
        'rm -r "${{tmp_indir}}/align_input/${{sample}}__captus-ext" ; '
        "fi ; "
        "done < {input.discarded_samples} && "
        "captus_assembly align "
        "--captus_extractions_dir ${{tmp_indir}}/align_input "
        "--out ${{tmp_outdir}} "
        '--ram "$(( {resources.mem_mb}/1000 ))" '
        "--threads {threads} "
        "--align_method {wildcards.align_method} "
        "--clipkit_gaps {wildcards.clipkit_gap} "
        "--min_coverage {wildcards.min_coverage} "
        "--markers {wildcards.marker} "
        "--format {wildcards.marker_format} "
        "&> {log} ; "
        "tar -cvf {output.tarfile} --directory ${{tmp_outdir}} ."


rule wscore_cutoffs:
    input:
        extraction_dir=Path(outdir, "005_selected-samples"),
    output:
        discarded_samples=Path(
            outdir,
            "007_wscore-filtering",
            "{marker}.wscore{sample_wscore_cutoff}",
            "discarded_samples.txt",
        ),
        discarded_loci=Path(
            outdir,
            "007_wscore-filtering",
            "{marker}.wscore{sample_wscore_cutoff}",
            "discarded_loci.txt",
        ),
    params:
        stats_file=lambda wildcards, input: Path(
            input.extraction_dir, "captus-assembly_extract.stats.tsv"
        ),
    log:
        log=Path(
            logdir,
            "wscore_cutoffs",
            "{marker}.wscore{sample_wscore_cutoff}" + ".log",
        ),
    container:
        r
    script:
        "src/wscore_cutoffs.R"


# This is very slow. Subset the extractions to speed it up.
# this is on scratch. it's probably better to copy once than symlink.
rule set_up_selected_samples:
    input:
        extraction_dir=captus_extractions_directory,
    output:
        outdir=directory(Path(outdir, "005_selected-samples")),
    params:
        n_samples=n_samples,
    run:
        captus_extractions_directory = Path(input.extraction_dir)
        n_samples = params.n_samples
        Path(output.outdir).mkdir(exist_ok=True)

        all_sample_folders = captus_extractions_directory.glob("*__captus-ext")
        auxfiles = captus_extractions_directory.glob(
            "captus-assembly_extract*"
        )

        all_sample_names = sorted(
            set(
                x.name.split("__")[0] for x in all_sample_folders if x.is_dir()
            )
        )

        random.seed(14)
        selected_samples = sorted(
            set(random.sample(all_sample_names, n_samples))
        )

        for sample in selected_samples:
            abs_from = Path(
                input.extraction_dir, f"{sample}__captus-ext"
            ).resolve()
            abs_to = Path(output.outdir, f"{sample}__captus-ext").resolve()
            shutil.copytree(abs_from, abs_to)

        for auxfile in auxfiles:
            shutil.copy(auxfile, output.outdir)


###########
# TARGETS #
###########


# Aggregate plots. Run this after the workflow finishes. There is no input for
# this rule, because some parameter combinations fail. The script lists files
# in the directories and plots whatever is there.
rule combine_trimal_and_iqtree_stats:
    output:
        combined_stats=Path(
            outdir,
            "035_combined-stats",
            "combined_stats.csv",
        ),
        naive_tree_score_plot=Path(
            outdir,
            "035_combined-stats",
            "naive_tree_score_plot.pdf",
        ),
        normalised_branch_length_plot=Path(
            outdir,
            "035_combined-stats",
            "normalised_branch_length_plot.pdf",
        ),
        gap_score_plot=Path(
            outdir,
            "035_combined-stats",
            "gap_score_plot.pdf",
        ),
        normalised_informative_site_plot=Path(
            outdir,
            "035_combined-stats",
            "normalised_informative_site_plot.pdf",
        ),
    params:
        trimal_stats_path=Path(outdir, "023_trimal-stats"),
        iqtree_stats_path=Path(outdir, "033_iqtree-stats"),
    log:
        log=Path(
            logdir,
            "combine_trimal_and_iqtree_stats.log",
        ),
    container:
        r
    script:
        "src/combine_trimal_and_iqtree_stats.R"


# actual targets
rule target:
    default_target: True
    input:
        ruleio_expand(rules.iqtree_target.input),
        ruleio_expand(rules.parse_iqtree_file.output),
        ruleio_expand(rules.parse_trimal_stats.output),
        ruleio_expand(rules.get_tip_length.output),
