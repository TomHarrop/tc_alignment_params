#!/usr/bin/env python3

from functools import cache
from pathlib import Path
import numpy as np
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


###########
# GLOBALS #
###########

# Paths
data = Path("data")
outdir = Path("output", "090_rerun-full-dataset")
logdir = Path(outdir, "logs")
benchdir = Path(logdir, "benchmarks")

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
markers = ["NUC", "PTD", "MIT", "DNA", "CLR", "ALL"]
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

# containers
ete3 = "docker://quay.io/biocontainers/ete3:3.1.1--py36_0"
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
    clipkit_gap="|".join(
        map(str, [round(float(x), 1) for x in np.arange(0, 1, 0.1)])
    ),
    min_coverage="|".join(
        map(str, [round(float(x), 1) for x in np.arange(0, 1, 0.1)])
    ),
    marker="|".join(markers),
    marker_format="|".join(["AA", "NT", "GE", "GF", "MA", "MF", "ALL"]),


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
        mem_mb=int(384e3),
        time="7-00",
    threads: 64
    retries: 0


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
    threads: 24
    resources:
        time=lambda wildcards, attempt: 60 * attempt,
    container:
        "docker://quay.io/biocontainers/biopython:1.81"
    script:
        "src/process_trimal_files.py"


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
        extraction_dir=Path(data, "03_extractions"),
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
        time=lambda wildcards, attempt: 480 * attempt,
        mem_mb=lambda wildcards, attempt: int(16e3 * attempt),
    container:
        "docker://quay.io/biocontainers/captus:1.0.1--pyhdfd78af_2"
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
        extraction_dir=Path(data, "03_extractions"),
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


###########
# TARGETS #
###########


rule target:
    default_target: True
    input:
        # defaults:
        # NUC.NT.mafft_auto.gap0.9.cov0.4.wscore0.0
        # good according to naive metric (before wscore cutoff):
        # NUC.NT.muscle_align.gap0.8.cov0.8.wscore0.0
        # best according to full param explore:
        # NUC.NT.muscle_align.gap0.8.cov0.8.wscore0.4
        expand(
            [
                str(x)
                for x in rules.iqtree_target.input
                + rules.get_tip_length.output
            ],
            zip,
            marker=["NUC", "NUC", "NUC"],
            marker_format=["NT", "NT", "NT"],
            align_method=["mafft_auto", "muscle_align", "muscle_align"],
            clipkit_gap=[0.9, 0.8, 0.8],
            min_coverage=[0.4, 0.8, 0.8],
            sample_wscore_cutoff=[0.0, 0.0, 0.4],
        ),
