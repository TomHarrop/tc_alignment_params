#!/usr/bin/env python3

from pathlib import Path
from snakemake import logger
import logging
import pandas as pd


def parse_site_data(site_data):
    table_data = []
    for line in site_data:
        try:
            parts = line.strip().split()
            if len(parts) == 9:
                entry = {
                    "Type": str(parts[1]),
                    "Seq": int(parts[2]),
                    "Site": int(parts[3]),
                    "Unique": int(parts[4]),
                    "Infor": int(parts[5]),
                    "Invar": int(parts[6]),
                    "Const": int(parts[7]),
                    "Name": str(parts[8]),
                }
                table_data.append(entry)
        except Exception as e:
            print(f"Error parsing line: {line}")
            raise e
    return pd.DataFrame(table_data)


def read_iqtree_file(file_path):
    site_data = []
    summary_metrics = {}
    table_started = False
    table_finished = False
    in_stats_region = False
    stats_region_finished = False

    with open(file_path, "r") as file:
        for line in file:
            if "ID\tType\tSeq\tSite\tUnique\tInfor\tInvar\tConst\tName" in line:
                table_started = True
                continue
            if table_started and not table_finished:
                if not line == "\n":
                    site_data.append(line)
                else:
                    table_finished = True
            else:
                if line == "MAXIMUM LIKELIHOOD TREE\n":
                    in_stats_region = True
                    continue
                if (
                    in_stats_region
                    and not stats_region_finished
                    and not line.startswith("WARNING")
                ):
                    parts = line.strip().split(":")
                    if len(parts) == 2:
                        try:
                            parsed_value = float(parts[1].split()[0].strip())
                            summary_metrics[parts[0]] = parsed_value
                        except ValueError as e:
                            logger.info(f"{line} is not a metric")
                    elif line.startswith("+"):
                        stats_region_finished = True

    parsed_site_data = parse_site_data(site_data)

    return parsed_site_data, summary_metrics


def main():
    parsed_data, summary_metrics = read_iqtree_file(iqtree_file)
    summary_metrics_df = pd.DataFrame(
        summary_metrics.items(), columns=["Metric", "Value"]
    )
    parsed_data.to_csv(site_data_file, index=False)
    summary_metrics_df.to_csv(summary_metrics_file, index=False)


# Example usage
# iqtree_file = Path("test/iqtree_files/tree.iqtree")
# site_data_file = Path("test", "site_data.csv")
# summary_metrics_file = Path("test", "summary_metrics.csv")
# main()

if __name__ == "__main__":
    iqtree_file = Path(snakemake.input["iqtree_file"])
    summary_metrics_file = Path(snakemake.output["summary_metrics"])
    site_data_file = Path(snakemake.output["site_data"])

    logfile = Path(snakemake.log[0])
    file_handler = logging.FileHandler(logfile)
    logger.logfile_handler = file_handler
    logger.logger.addHandler(logger.logfile_handler)

    try:
        main()
    except Exception as e:
        logger.error(e)
        raise e
