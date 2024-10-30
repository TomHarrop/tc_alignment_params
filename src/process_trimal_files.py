#!/usr/bin/env python3

from Bio import AlignIO
from pathlib import Path
from snakemake import logger
import logging
import tarfile
import tempfile


def check_for_informative_sites(trimal_file):
    with open(trimal_file, "r") as f:
        alignment = AlignIO.read(f, "fasta")
    alphabet = alignment.substitutions.alphabet
    substitutions = alignment.substitutions.select(alphabet)

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
    trimal_path = Path(tempfile.mkdtemp())

    with tarfile.open(trimal_tarfile, "r") as tar:
        tar.extractall(trimal_path)

    trimal_files = sorted(
        set(x for x in trimal_path.glob("*.fna") if not x.name.startswith("."))
    )

    # open the tarfiles for writing
    kept_tar = tarfile.open(kept_tarfile, "w")
    discarded_tar = tarfile.open(discarded_tarfile, "w")
    empty_tar = tarfile.open(empty_tarfile, "w")

    for trimal_file in trimal_files:
        try:
            alignment, informative_sites = check_for_informative_sites(trimal_file)
        except ValueError:
            # this is an empty file
            empty_tar.add(trimal_file, arcname=trimal_file.name)
            continue
        if informative_sites:
            kept_tar.add(trimal_file, arcname=trimal_file.name)
        else:
            discarded_tar.add(trimal_file, arcname=trimal_file.name)

    kept_tar.close()
    discarded_tar.close()
    empty_tar.close()


# testing
# trimal_tarfile = Path("test", "trimal.tar")
# kept_tarfile = Path("test", "kept.tar")
# discarded_tarfile = Path("test", "discarded.tar")
# empty_tarfile = Path("test", "empty.tar")
# main()


if __name__ == "__main__":

    trimal_tarfile = Path(snakemake.input["tarfile"])
    kept_tarfile = Path(snakemake.output["kept"])
    discarded_tarfile = Path(snakemake.output["discarded"])
    empty_tarfile = Path(snakemake.output["empty"])

    logfile = Path(snakemake.log[0])
    file_handler = logging.FileHandler(logfile)
    logger.logfile_handler = file_handler
    logger.logger.addHandler(logger.logfile_handler)

    try:
        main()
    except Exception as e:
        logger.error(e)
        raise e
