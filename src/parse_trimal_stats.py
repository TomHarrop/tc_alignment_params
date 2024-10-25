#!/usr/bin/env python3

# TODO. This script parses the trimal stats output. It needs to go in the
# trimal module.

import pandas as pd
import statistics
from pathlib import Path


def read_trimal_file(file_path):
    gap_table_data = []
    summary_metrics = {}
    table_started = False
    table_finished = False
    with open(file_path, "r") as file:
        for line in file:
            if "|	Number of Gaps" in line:
                table_started = True
                continue
            if table_started and not table_finished:
                # Skip the separator line
                if line.startswith("| Residues") or line.startswith("+"):
                    continue
                elif not line.startswith("|"):
                    gap_table_data.append(line)
                elif line.startswith("|"):
                    table_finished = True
            else:
                if line.startswith("##"):
                    parts = line.strip().split()
                    if len(parts) == 3:
                        summary_metrics[parts[1]] = float(parts[2])

    parsed_gap_table = parse_table_data(gap_table_data)
    MedianGapScore, AverageGapScore = summarise_gap_scores(parsed_gap_table)

    summary_metrics["MedianGapScore"] = MedianGapScore
    summary_metrics["AverageGapScore"] = AverageGapScore

    return summary_metrics


def parse_table_data(gap_table_data):
    data = []
    for line in gap_table_data:
        parts = line.strip().split()
        if len(parts) == 7:
            entry = {
                "Number_of_Residues": int(parts[0]),
                "Pct_Length": float(parts[1]),
                "Cumulative_Number_of_Residues": int(parts[2]),
                "Pct_Cumulative_Length": float(parts[3]),
                "Number_of_Gaps_per_Column": int(parts[4]),
                "Pct_Gaps_per_Column": float(parts[5]),
                "Gap_Score_per_Column": float(parts[6]),
            }
            data.append(entry)
    return pd.DataFrame(data)


def summarise_gap_scores(df):
    decoded_sequence = []
    for _, row in df.iterrows():
        decoded_sequence.extend(
            [row["Gap_Score_per_Column"]] * int(row["Number_of_Residues"])
        )
    MedianGapScore = float(statistics.median(decoded_sequence))
    AverageGapScore = float(statistics.mean(decoded_sequence))
    return MedianGapScore, AverageGapScore


trimal_log_dir = Path("output/020_trimal/NUC.NT.muscle_super5.gap0.8.cov0.4/logs")
trimal_stats_files = trimal_log_dir.glob("trimal.*.stats.txt")

[print(file) for file in trimal_stats_files]

summary_metrics = {}
for file_path in trimal_stats_files:
    id = file_path.name.split(".")[1]
    summary_metrics[id] = read_trimal_file(file_path)

summary_metrics_df = pd.DataFrame(summary_metrics).T