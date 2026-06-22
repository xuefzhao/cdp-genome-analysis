#!/usr/bin/env python3
import argparse
import csv
from collections import defaultdict
from pathlib import Path
from typing import Dict, Set


VARIANT_CLASSES = ("synonymous", "lof", "missense", "others")


def parse_carriers(raw: str) -> Set[str]:
    if not raw or raw == ".":
        return set()
    return {sample.strip() for sample in raw.split(",") if sample.strip() and sample != "."}


def canonical_class(raw: str) -> str:
    value = (raw or "").strip().lower()
    if value == "other":
        return "others"
    if value in VARIANT_CLASSES:
        return value
    return ""


def init_counts() -> Dict[str, int]:
    return {"het_count": 0, "hom_count": 0, "all_count": 0}


def count_per_class(bed_tsv: str) -> Dict[str, Dict[str, Dict[str, int]]]:
    per_class: Dict[str, Dict[str, Dict[str, int]]] = {
        cls: defaultdict(init_counts) for cls in VARIANT_CLASSES
    }

    with open(bed_tsv, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"impact_class", "het_carriers", "hom_carriers"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing required columns in BED TSV: {', '.join(sorted(missing))}")

        for row in reader:
            vclass = canonical_class(row["impact_class"])
            if not vclass:
                continue

            het_samples = parse_carriers(row["het_carriers"])
            hom_samples = parse_carriers(row["hom_carriers"])

            for sample in het_samples:
                per_class[vclass][sample]["het_count"] += 1
                per_class[vclass][sample]["all_count"] += 1

            for sample in hom_samples:
                per_class[vclass][sample]["hom_count"] += 1
                per_class[vclass][sample]["all_count"] += 1

    return per_class


def write_table(counts: Dict[str, Dict[str, int]], out_tsv: str) -> None:
    with open(out_tsv, "wt", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample_id", "het_count", "hom_count", "all_count"])
        for sample_id in sorted(counts):
            row = counts[sample_id]
            writer.writerow([sample_id, row["het_count"], row["hom_count"], row["all_count"]])


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Count per-sample variant counts by class (synonymous/lof/missense/others) "
            "from BED TSV output by vcf_to_bed_with_carriers.py."
        )
    )
    parser.add_argument("--bed", required=True, help="Input BED TSV file")
    parser.add_argument(
        "--out-prefix",
        required=True,
        help=(
            "Output prefix. Produces: "
            "<prefix>.synonymous.stat.tsv, <prefix>.lof.stat.tsv, "
            "<prefix>.missense.stat.tsv, <prefix>.others.stat.tsv"
        ),
    )
    args = parser.parse_args()

    counts = count_per_class(args.bed)
    prefix = Path(args.out_prefix)
    for cls in VARIANT_CLASSES:
        write_table(counts[cls], str(prefix) + f".{cls}.stat.tsv")


if __name__ == "__main__":
    main()
