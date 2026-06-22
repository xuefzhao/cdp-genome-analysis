#!/usr/bin/env python3
import argparse
import csv
from collections import defaultdict
from typing import Dict, List


def merge_tables(input_tables: List[str]) -> Dict[str, Dict[str, int]]:
    merged: Dict[str, Dict[str, int]] = defaultdict(
        lambda: {"het_count": 0, "hom_count": 0, "all_count": 0}
    )

    required = {
        "sample_id",
        "het_count",
        "hom_count",
        "all_count",
    }

    for table in input_tables:
        with open(table, "rt", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            missing = required - set(reader.fieldnames or [])
            if missing:
                raise ValueError(
                    f"Missing required columns in {table}: {', '.join(sorted(missing))}"
                )

            for row in reader:
                sample = row["sample_id"]
                if not sample:
                    continue
                merged[sample]["het_count"] += int(row["het_count"])
                merged[sample]["hom_count"] += int(row["hom_count"])
                merged[sample]["all_count"] += int(row["all_count"])

    return merged


def write_merged(merged: Dict[str, Dict[str, int]], output_table: str) -> None:
    with open(output_table, "wt", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample_id", "het_count", "hom_count", "all_count"])
        for sample in sorted(merged):
            row = merged[sample]
            writer.writerow(
                [
                    sample,
                    row["het_count"],
                    row["hom_count"],
                    row["all_count"],
                ]
            )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Merge per-contig class-specific count tables into one per-sample table."
    )
    parser.add_argument(
        "--tables",
        required=True,
        nargs="+",
        help="Input count tables (TSV) to merge",
    )
    parser.add_argument("--out", required=True, help="Output merged TSV table")
    args = parser.parse_args()

    merged = merge_tables(args.tables)
    write_merged(merged, args.out)


if __name__ == "__main__":
    main()
