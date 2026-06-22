#!/usr/bin/env python3
import argparse
import gzip
from typing import Dict, List, Optional, Tuple


LOF_CONSEQUENCES = {
    "transcript_ablation",
    "exon_loss_variant",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "start_lost",
    "stop_lost",
}

MISSENSE_CONSEQUENCES = {
    "missense_variant",
    "inframe_insertion",
    "inframe_deletion",
    "protein_altering_variant",
    "coding_sequence_variant",
}

SYNONYMOUS_CONSEQUENCES = {
    "synonymous_variant",
    "stop_retained_variant",
}

IMPACT_RANK = {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3}


def open_text(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def parse_info(info_str: str) -> Dict[str, str]:
    info: Dict[str, str] = {}
    for item in info_str.split(";"):
        if not item:
            continue
        if "=" in item:
            key, val = item.split("=", 1)
            info[key] = val
        else:
            info[item] = "1"
    return info


def allele_specific(value: str, allele_index: int) -> str:
    if not value or value == ".":
        return "."
    parts = value.split(",")
    if len(parts) == 1:
        return parts[0] if parts[0] else "."
    idx = allele_index - 1
    if 0 <= idx < len(parts):
        return parts[idx] if parts[idx] else "."
    return "."


def parse_vep_fields(header_line: str) -> List[str]:
    marker = "Format: "
    if marker not in header_line:
        return []
    fmt = header_line.split(marker, 1)[1].rstrip('">')
    return fmt.split("|")


def classify_impact(consequence: str) -> str:
    terms = set()
    for token in consequence.split("&"):
        token = token.strip()
        if token:
            terms.add(token)
    if terms & LOF_CONSEQUENCES:
        return "lof"
    if terms & MISSENSE_CONSEQUENCES:
        return "missense"
    if terms & SYNONYMOUS_CONSEQUENCES:
        return "synonymous"
    return "others"


def select_vep_annotation(
    vep_value: str, alt: str, allele_index: int, vep_fields: List[str]
) -> Tuple[str, str, str, str]:
    if not vep_value or vep_value == "." or not vep_fields:
        return ".", ".", ".", "others"

    idx_map = {field: i for i, field in enumerate(vep_fields)}
    allele_i = idx_map.get("Allele")
    impact_i = idx_map.get("IMPACT")
    symbol_i = idx_map.get("SYMBOL")
    consequence_i = idx_map.get("Consequence")
    allele_num_i = idx_map.get("ALLELE_NUM")
    pick_i = idx_map.get("PICK")
    canonical_i = idx_map.get("CANONICAL")

    best = None
    best_score = (99, 1, 1)  # impact rank, pick penalty, canonical penalty

    for ann in vep_value.split(","):
        fields = ann.split("|")
        if len(fields) < len(vep_fields):
            fields += [""] * (len(vep_fields) - len(fields))

        ann_allele = fields[allele_i] if allele_i is not None else ""
        ann_allele_num = fields[allele_num_i] if allele_num_i is not None else ""

        allele_match = ann_allele == alt
        if not allele_match and ann_allele_num.isdigit():
            allele_match = int(ann_allele_num) == allele_index
        if not allele_match:
            continue

        impact = fields[impact_i] if impact_i is not None and fields[impact_i] else "."
        consequence = (
            fields[consequence_i]
            if consequence_i is not None and fields[consequence_i]
            else "."
        )
        symbol = fields[symbol_i] if symbol_i is not None and fields[symbol_i] else "."
        pick_penalty = 0 if (pick_i is not None and fields[pick_i] == "1") else 1
        canonical_penalty = (
            0 if (canonical_i is not None and fields[canonical_i] == "YES") else 1
        )
        score = (IMPACT_RANK.get(impact, 99), pick_penalty, canonical_penalty)
        if score < best_score:
            best_score = score
            best = (impact, consequence, symbol, classify_impact(consequence))

    if best is None:
        return ".", ".", ".", "others"
    return best


def split_carriers(
    sample_names: List[str], format_str: str, sample_values: List[str], allele_index: int
) -> Tuple[str, str]:
    format_keys = format_str.split(":")
    gt_idx = format_keys.index("GT") if "GT" in format_keys else -1
    if gt_idx == -1:
        return ".", "."

    het: List[str] = []
    hom: List[str] = []
    for sample, value in zip(sample_names, sample_values):
        parts = value.split(":")
        gt = parts[gt_idx] if gt_idx < len(parts) else "."
        if not gt or gt == ".":
            continue

        alleles = []
        for token in gt.replace("|", "/").split("/"):
            if token.isdigit():
                alleles.append(int(token))
        if not alleles:
            continue

        target_count = sum(1 for a in alleles if a == allele_index)
        if target_count == 0:
            continue
        if target_count == len(alleles):
            hom.append(sample)
        else:
            het.append(sample)

    return ",".join(het) if het else ".", ",".join(hom) if hom else "."


def convert_vcf_to_bed(vcf_path: str, out_path: str) -> None:
    vep_fields: List[str] = []
    sample_names: List[str] = []

    header = [
        "chr",
        "pos",
        "end",
        "variantID",
        "ref",
        "alt",
        "size",
        "FILTER",
        "AF",
        "AC",
        "AN",
        "gnomAD_srWGS_match_type",
        "gnomAD_srWGS_match_source",
        "gnomAD_srWGS_match_filter",
        "gnomAD_srWGS_match_ID",
        "gnomAD_srWGS_match_AF",
        "gnomAD_srWES_match_type",
        "gnomAD_srWES_match_source",
        "gnomAD_srWES_match_filter",
        "gnomAD_srWES_match_ID",
        "gnomAD_srWES_match_AF",
        "IMPACT",
        "Consequence",
        "SYMBOL",
        "impact_class",
        "het_carriers",
        "hom_carriers",
    ]

    with open_text(vcf_path) as fin, open(out_path, "wt") as fout:
        fout.write("\t".join(header) + "\n")

        for line in fin:
            if line.startswith("##INFO=<ID=vep,"):
                vep_fields = parse_vep_fields(line.strip())
                continue
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                cols = line.rstrip("\n").split("\t")
                sample_names = cols[9:]
                continue

            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue

            chrom, pos, variant_id, ref, alt_str, _qual, filt, info_str = cols[:8]
            format_str = cols[8] if len(cols) > 8 else ""
            sample_values = cols[9:] if len(cols) > 9 else []
            info = parse_info(info_str)

            alts = alt_str.split(",")
            pos_int = int(pos)
            start = pos_int - 1
            end = start + len(ref)

            for i, alt in enumerate(alts, start=1):
                impact, consequence, symbol, impact_class = select_vep_annotation(
                    info.get("vep", "."), alt, i, vep_fields
                )
                het, hom = split_carriers(sample_names, format_str, sample_values, i)

                out_row = [
                    chrom,
                    str(start),
                    str(end),
                    variant_id if variant_id else ".",
                    ref,
                    alt,
                    str(len(alt) - len(ref)),
                    filt if filt else ".",
                    allele_specific(info.get("AF", "."), i),
                    allele_specific(info.get("AC", "."), i),
                    info.get("AN", "."),
                    info.get("gnomAD_srWGS_match_type", "."),
                    info.get("gnomAD_srWGS_match_source", "."),
                    info.get("gnomAD_srWGS_match_filter", "."),
                    info.get("gnomAD_srWGS_match_ID", "."),
                    allele_specific(info.get("gnomAD_srWGS_match_AF", "."), i),
                    info.get("gnomAD_srWES_match_type", "."),
                    info.get("gnomAD_srWES_match_source", "."),
                    info.get("gnomAD_srWES_match_filter", "."),
                    info.get("gnomAD_srWES_match_ID", "."),
                    allele_specific(info.get("gnomAD_srWES_match_AF", "."), i),
                    impact,
                    consequence,
                    symbol,
                    impact_class,
                    het,
                    hom,
                ]
                fout.write("\t".join(out_row) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert an annotated VCF(.gz) to BED-like TSV with VEP and carrier columns."
    )
    parser.add_argument("--vcf", required=True, help="Input VCF or VCF.GZ path")
    parser.add_argument("--out", required=True, help="Output BED/TSV path")
    args = parser.parse_args()

    convert_vcf_to_bed(args.vcf, args.out)


if __name__ == "__main__":
    main()
