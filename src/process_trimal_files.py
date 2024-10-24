#!/usr/bin/env python3

from Bio import AlignIO
from pathlib import Path
from snakemake import logger
import logging


def check_for_informative_sites(trimal_file):
    with open(trimal_file, "r") as f:
        alignment = AlignIO.read(f, "fasta")

    substitutions = alignment.substitutions.select("ACTG")

    informative_sites = False
    for i, row in enumerate(substitutions):
        for j, value in enumerate(row):
            if i != j and value > 0:
                informative_sites = True
                break
        if informative_sites:
            break
    return alignment, informative_sites


def main():
    trimal_files = sorted(
        set(x for x in trimal_path.glob("*.fasta") if not x.name.startswith("."))
    )

    Path(output_path, "kept").mkdir(parents=True, exist_ok=True)
    Path(output_path, "discarded").mkdir(parents=True, exist_ok=True)
    Path(output_path, "empty").mkdir(parents=True, exist_ok=True)

    for trimal_file in trimal_files:
        try:
            alignment, informative_sites = check_for_informative_sites(trimal_file)
        except ValueError:
            # this is an empty file
            output_file = Path(output_path, "empty", trimal_file.name)
            output_file.touch()
            continue
        if informative_sites:
            output_file = Path(output_path, "kept", trimal_file.name)
        else:
            output_file = Path(output_path, "discarded", trimal_file.name)
        with open(output_file, "w") as f:
            AlignIO.write(alignment, f, "fasta")


if __name__ == "__main__":

    trimal_path = Path(snakemake.input["trimal_path"])
    output_path = Path(snakemake.params["output_path"])

    logfile = Path(snakemake.log[0])
    file_handler = logging.FileHandler(logfile)
    logger.logfile_handler = file_handler
    logger.logger.addHandler(logger.logfile_handler)

    try:
        main()
    except Exception as e:
        logger.error(e)
        raise e
