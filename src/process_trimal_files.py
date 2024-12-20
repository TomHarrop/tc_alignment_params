#!/usr/bin/env python3

from Bio import AlignIO
from pathlib import Path
from snakemake import logger
import logging
import multiprocessing as mp
import shutil
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
    return informative_sites


# This function monitors the queue and writes new items to the tarfile
def write_to_tarfile(queue, output_tarfile):
    i = 0
    logger.debug(f"This is the queue for {queue}, writing to {output_tarfile}")
    with tarfile.open(output_tarfile, "w") as tar:
        while True:
            item = queue.get()
            logger.debug(f"Got {item}")
            if item is None:
                logger.info(f"Sorted a total of {i} alignments to {output_tarfile}")
                break
            tar.add(item, arcname=item.name)
            i += 1
            if i % 10 == 0:
                logger.info(f"{i} alignments sorted to {output_tarfile}...")


# This function sends the trimal file to the appropriate queue
def process_trimal_file(trimal_file, write_queues):
    try:
        informative_sites = check_for_informative_sites(trimal_file)
        if informative_sites:
            write_queues["kept"].put(trimal_file)
        else:
            write_queues["discarded"].put(trimal_file)
    except ValueError:
        # this is an empty file
        write_queues["empty"].put(trimal_file)


def main():
    output_streams = ["kept", "discarded", "empty"]

    # set up pool
    logger.info(f"Processing trimal output with {threads} threads")
    pool = mp.Pool(threads + len(output_streams))

    # set up listeners
    manager = mp.Manager()
    write_queues = {k: manager.Queue() for k in output_streams}

    temp_tarfiles = {k: Path(tempfile.mktemp(suffix=".tar")) for k in output_streams}

    for k, v in temp_tarfiles.items():
        logger.info(f"Sorting {k} alignments to {v}")

    # start the listeners
    watchers = {
        k: pool.apply_async(write_to_tarfile, args=(write_queues[k], temp_tarfiles[k]))
        for k in output_streams
    }
    logger.debug(watchers)

    # get the input alignments
    trimal_path = Path(tempfile.mkdtemp())
    logger.info(f"Extracting trimal files to {trimal_path}")
    with tarfile.open(trimal_tarfile, "r") as tar:
        tar.extractall(trimal_path)

    trimal_files = sorted(
        set(x for x in trimal_path.glob("*.fna") if not x.name.startswith("."))
    )
    logger.info(f"Found {len(trimal_files)} alignment files")

    # launch a job for each alignment
    logger.info("Starting alignment processing")
    jobs = []
    for trimal_file in trimal_files:
        job = pool.apply_async(process_trimal_file, args=(trimal_file, write_queues))
        jobs.append(job)

    logger.debug(f"Queued {len(jobs)} jobs")

    # wait for all jobs to finish
    for job in jobs:
        job.get()

    logger.debug("All jobs finished")

    # close the listeners
    logger.debug("Closing write queues")
    for k, v in write_queues.items():
        v.put(None)

    logger.debug("Closing pools")
    pool.close()
    pool.join()

    # move the temp tarfiles to the final tarfiles
    logger.info("Writing output tarfiles")
    for k, v in temp_tarfiles.items():
        logger.info(f"Writing {k} alignments to {output_tarfiles[k]}")
        shutil.move(v, output_tarfiles[k])

    logger.info("Done")

# testing
# trimal_tarfile = Path("test", "trimal.tar")
# output_tarfiles = {
#     "kept": Path("test", "kept.tar"),
#     "discarded": Path("test", "discarded.tar"),
#     "empty": Path("test", "empty.tar"),
# }
# threads = 11
# main()


if __name__ == "__main__":

    trimal_tarfile = Path(snakemake.input["tarfile"])
    output_tarfiles = {
        "kept": Path(snakemake.output["kept"]),
        "discarded": Path(snakemake.output["discarded"]),
        "empty": Path(snakemake.output["empty"]),
    }

    threads = snakemake.threads

    logfile = Path(snakemake.log[0])
    file_handler = logging.FileHandler(logfile)
    logger.logfile_handler = file_handler
    logger.logger.addHandler(logger.logfile_handler)

    try:
        main()
    except Exception as e:
        logger.error(e)
        raise e
