#!/usr/bin/env python3

from functools import cache
from pathlib import Path
import numpy as np

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
clipkit_gaps = np.arange(0, 1, 0.1)
min_coverages = np.arange(0, 1, 0.1)
markers = ["NUC", "PTD", "MIT", "DNA", "CLR", "ALL"]
formats = ["AA", "NT", "GE", "GF", "MA", "MF", "ALL"]

param_string = (
    "{marker}."
    "{marker_format}."
    "{align_method}."
    "gap{clipkit_gap}."
    "cov{min_coverage}"
)

# containers
captus = "docker://quay.io/biocontainers/captus:1.0.1--pyhdfd78af_2"

# modules
module_tag = "0.7.0"

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


module iqtree:
    snakefile:
        iqtree_snakefile
    config:
        {
            "alignment_directory": Path(
                outdir,
                "020_trimal",
                param_string,
                "trimmed",
            ),
            "outdir": Path(
                outdir,
                "030_iqtree",
                param_string,
            ),
        }


use rule * from iqtree as iqtree_*


module trimal:
    snakefile:
        trimal_snakefile
    config:
        {
            "alignment_directory": Path(
                outdir,
                "020_trimal",
                param_string,
                "input",
            ),
            "outdir": Path(
                outdir,
                "020_trimal",
                param_string,
            ),
        }


use rule * from trimal as trimal_*


rule collect_captus_alignment_directories:
    input:
        alignment_directory=Path(
            outdir,
            "010_captus-align",
            param_string,
        ),
    output:
        directory(
            Path(
                outdir,
                "020_trimal",
                param_string,
                "input",
            )
        ),
    params:
        alignment_directory=lambda wildcards, input: Path(
            input.alignment_directory,
            "04_alignments",
        ),
    shell:
        "ln -sf "
        "$( readlink -f {params.alignment_directory} ) "
        "$( readlink -f {output} )"


rule captus_align:
    input:
        extraction_dir=captus_extractions_directory,
    output:
        outdir=directory(
            Path(
                outdir,
                "010_captus-align",
                param_string,
            )
        ),
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
    threads: lambda wildcards, attempt: 16 * attempt
    resources:
        time=lambda wildcards, attempt: 120 * attempt,
        mem_mb=lambda wildcards, attempt: int(16e3 * attempt),
    container:
        captus
    shell:
        "tmp_indir=$(mktemp -d) ; "
        "echo $tmp_indir ; "
        "exit 1 ; "
        "tmp_outdir=$(mktemp -d) ; "
        "cp -r {input.extraction_dir} ${{tmp_indir}}/align_input ; "
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
        "mv ${{tmp_outdir}} {output.outdir}"



###########
# TARGETS #
###########


rule target:
    default_target: True
    input:
        expand(
            [str(x) for x in rules.iqtree_target.input],
            align_method=align_methods,
            clipkit_gap=clipkit_gaps,
            min_coverage=min_coverages,
            marker=["NUC"],             # just for now
            marker_format=["NT"],       # ditto
        )[0],
