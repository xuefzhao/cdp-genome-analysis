#!/usr/bin/env python3
"""
classify_variants.py
────────────────────
Classify VCF variants relative to gene/exon features from a GTF(.gz) file.

Classes assigned (in priority order):
  5UTR          – overlaps a 5' UTR feature
  3UTR          – overlaps a 3' UTR feature
  coding        – overlaps a CDS feature
  intronic      – falls inside a gene but outside any exon/UTR/CDS
  intergenic_upstream   – upstream of the nearest gene on the same strand
  intergenic_downstream – downstream of the nearest gene on the same strand

For intergenic variants:
  dist_closest  – distance to the nearest gene boundary
  dist_adjacent – total span between the two flanking genes
                  (0 if only one gene exists on that chromosome)

For intronic variants:
  dist_closest  – distance to the nearest exon boundary
  dist_adjacent – span between the two flanking exons
                  (0 if variant is before the first or after the last exon)

For 5UTR / 3UTR variants:
  dist_closest  – distance to the nearest CDS exon boundary
  dist_adjacent – '.' (not applicable)

Output: BED-like TSV with columns:
  chrom  start(0-based)  end  ID  REF  ALT  FILTER
  class  dist_closest  dist_adjacent

Usage:
  python classify_variants.py \\
      --vcf   input.vcf[.gz] \\
      --gtf   annotation.gtf[.gz] \\
      --out   output.bed \\
      [--upstream-bp   2000]   # max bp upstream to call 'intergenic_upstream'
      [--downstream-bp 2000]   # max bp downstream to call 'intergenic_downstream'

Dependencies: Python ≥ 3.8, stdlib only (gzip, bisect, argparse, re, csv)
"""

import argparse
import bisect
import csv
import gzip
import re
import sys
from collections import defaultdict

# ── helpers ──────────────────────────────────────────────────────────

def open_file(path):
    """Open plain or gzip-compressed file."""
    if path.endswith('.gz'):
        return gzip.open(path, 'rt', encoding='utf-8')
    return open(path, 'r', encoding='utf-8')


def parse_gtf_attribute(attr_str, key):
    """Extract value for a GTF attribute key."""
    m = re.search(rf'{key}\s+"([^"]+)"', attr_str)
    return m.group(1) if m else '.'


# ── GTF loader ────────────────────────────────────────────────────────

def load_gtf(gtf_path):
    """
    Parse GTF and build per-chromosome sorted interval structures.

    Returns
    -------
    genes[chrom]  : sorted list of (start, end, gene_id, strand)  [0-based, half-open)
    exons[chrom]  : sorted list of (start, end, gene_id, feature)
                    feature ∈ {'exon', 'CDS', '5UTR', '3UTR'}
    gene_starts[chrom] : sorted list of gene start positions  (for bisect)
    gene_ends[chrom]   : sorted list of gene end positions    (for bisect)
    exon_starts[chrom] : sorted list of exon/CDS/UTR starts
    exon_ends[chrom]   : sorted list of exon/CDS/UTR ends
    """
    genes    = defaultdict(list)
    exons    = defaultdict(list)
    feat_set = {'CDS', 'exon', '5UTR', '3UTR',
                'five_prime_utr', 'three_prime_utr',
                'UTR'}          # Ensembl / GENCODE variations

    print(f"Loading GTF: {gtf_path}", file=sys.stderr)
    n_genes, n_exons = 0, 0

    with open_file(gtf_path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            chrom, source, feature, start_s, end_s, score, strand, frame, attrs = parts
            start = int(start_s) - 1   # GTF is 1-based; convert to 0-based
            end   = int(end_s)         # end stays as-is → half-open [start, end)
            gene_id = parse_gtf_attribute(attrs, 'gene_id')

            if feature == 'gene':
                genes[chrom].append((start, end, gene_id, strand))
                n_genes += 1
            elif feature in feat_set:
                # Normalise UTR names
                f = feature
                if f in ('five_prime_utr',  'UTR'):   # Ensembl uses 'UTR' for both
                    f = '5UTR'
                elif f == 'three_prime_utr':
                    f = '3UTR'
                exons[chrom].append((start, end, gene_id, f))
                n_exons += 1

    # Sort all lists by start
    for chrom in genes:
        genes[chrom].sort()
    for chrom in exons:
        exons[chrom].sort()

    # Pre-build sorted position arrays for bisect
    gene_starts = {c: [g[0] for g in v] for c, v in genes.items()}
    gene_ends   = {c: [g[1] for g in v] for c, v in genes.items()}
    exon_starts = {c: [e[0] for e in v] for c, v in exons.items()}
    exon_ends   = {c: [e[1] for e in v] for c, v in exons.items()}

    print(f"  Loaded {n_genes} gene features, {n_exons} exon/CDS/UTR features "
          f"across {len(genes)} chromosomes.", file=sys.stderr)
    return genes, exons, gene_starts, gene_ends, exon_starts, exon_ends


# ── interval query helpers ─────────────────────────────────────────────

def overlapping_genes(pos, chrom, genes, gene_starts, gene_ends):
    """Return list of gene records that contain pos (0-based)."""
    if chrom not in genes:
        return []
    starts = gene_starts[chrom]
    ends   = gene_ends[chrom]
    gs     = genes[chrom]

    # Candidates: genes whose start <= pos
    right = bisect.bisect_right(starts, pos)
    result = []
    for i in range(right - 1, -1, -1):
        g = gs[i]
        if g[1] > pos:   # end > pos  → overlaps
            result.append(g)
        elif g[0] < pos - 5_000_000:  # early exit: genes are far behind
            break
    return result


def overlapping_exons(pos, chrom, exons, exon_starts, exon_ends):
    """Return list of exon records that contain pos."""
    if chrom not in exons:
        return []
    starts = exon_starts[chrom]
    es     = exons[chrom]
    right  = bisect.bisect_right(starts, pos)
    result = []
    for i in range(right - 1, -1, -1):
        e = es[i]
        if e[1] > pos:
            result.append(e)
        elif e[0] < pos - 5_000_000:
            break
    return result


def flanking_genes(pos, chrom, genes, gene_starts, gene_ends):
    """
    Return (upstream_gene, downstream_gene) closest to pos.
    Each is a gene record tuple or None.
    upstream   = rightmost gene whose end  <= pos
    downstream = leftmost  gene whose start > pos
    """
    if chrom not in genes:
        return None, None
    gs     = genes[chrom]
    starts = gene_starts[chrom]
    ends   = gene_ends[chrom]

    # downstream: first gene whose start > pos
    idx_down = bisect.bisect_right(starts, pos)
    downstream = gs[idx_down] if idx_down < len(gs) else None

    # upstream: last gene whose end <= pos
    idx_up = bisect.bisect_right(ends, pos) - 1
    upstream = gs[idx_up] if idx_up >= 0 else None

    return upstream, downstream


def flanking_exons(pos, chrom, exons, exon_starts, exon_ends):
    """
    Return (prev_exon, next_exon) flanking pos within an intronic region.
    prev_exon: rightmost exon whose end  <= pos
    next_exon: leftmost  exon whose start > pos
    """
    if chrom not in exons:
        return None, None
    es     = exons[chrom]
    starts = exon_starts[chrom]
    ends   = exon_ends[chrom]

    idx_next = bisect.bisect_right(starts, pos)
    next_exon = es[idx_next] if idx_next < len(es) else None

    idx_prev = bisect.bisect_right(ends, pos) - 1
    prev_exon = es[idx_prev] if idx_prev >= 0 else None

    return prev_exon, next_exon


# ── variant classification ────────────────────────────────────────────

def dist_to_nearest_cds(pos, chrom, exons, exon_starts, exon_ends):
    """
    For a position inside a UTR, return the distance to the nearest CDS
    exon boundary (start or end).  Searches both the flanking exon before
    and after pos, filtering to CDS features only.

    Returns int distance, or '.' if no CDS features exist on this chrom.
    """
    if chrom not in exons:
        return '.'

    es     = exons[chrom]
    starts = exon_starts[chrom]
    ends   = exon_ends[chrom]

    # Candidate: last CDS whose start <= pos
    idx_right = bisect.bisect_right(starts, pos)

    distances = []

    # Look backwards for a CDS that might overlap or be just behind pos
    for i in range(idx_right - 1, max(idx_right - 200, -1), -1):
        e = es[i]
        if e[3] != 'CDS':
            continue
        # Distance from pos to end of this CDS (pos is past the CDS end)
        if e[1] <= pos:
            distances.append(pos - e[1])
        else:
            # CDS straddles pos (shouldn't happen for UTR variants, but handle it)
            distances.append(0)
        break   # nearest one behind is enough

    # Look forwards for the next CDS
    for i in range(idx_right, min(idx_right + 200, len(es))):
        e = es[i]
        if e[3] != 'CDS':
            continue
        distances.append(e[0] - pos)
        break

    return min(distances) if distances else '.'


def classify_variant(pos, chrom,
                     genes, exons,
                     gene_starts, gene_ends,
                     exon_starts, exon_ends,
                     upstream_bp=2000, downstream_bp=2000):
    """
    Classify a single variant at 0-based `pos` on `chrom`.

    Returns
    -------
    (classification, dist_closest, dist_adjacent)
    classification : str
    dist_closest   : int or '.'
                     – intronic/UTR: distance to nearest exon / CDS boundary
                     – intergenic:   distance to nearest gene boundary
                     – coding:       '.'
    dist_adjacent  : int or '.'
                     – intronic:    span between the two flanking exons
                     – intergenic:  span between the two flanking genes
                     – UTR/coding:  '.'
    """
    # 1. Check exon-level features first (UTR / CDS / exon)
    ov_exons = overlapping_exons(pos, chrom, exons, exon_starts, exon_ends)
    if ov_exons:
        # Priority: CDS > 5UTR > 3UTR > exon (exon without CDS = non-coding exon)
        feats = [e[3] for e in ov_exons]
        if 'CDS' in feats:
            return 'coding', '.', '.'
        if '5UTR' in feats:
            dist_c = dist_to_nearest_cds(pos, chrom, exons, exon_starts, exon_ends)
            return '5UTR', dist_c, '.'
        if '3UTR' in feats:
            dist_c = dist_to_nearest_cds(pos, chrom, exons, exon_starts, exon_ends)
            return '3UTR', dist_c, '.'
        return 'coding', '.', '.'   # exon overlap, no CDS annotation → treat as coding

    # 2. Check gene overlap → intronic
    ov_genes = overlapping_genes(pos, chrom, genes, gene_starts, gene_ends)
    if ov_genes:
        prev_ex, next_ex = flanking_exons(pos, chrom, exons, exon_starts, exon_ends)
        dist_closest = '.'
        dist_adjacent = '.'

        distances = []
        if prev_ex is not None:
            distances.append(pos - prev_ex[1])   # distance from pos to exon end
        if next_ex is not None:
            distances.append(next_ex[0] - pos)   # distance from pos to exon start
        if distances:
            dist_closest = min(distances)

        if prev_ex is not None and next_ex is not None:
            dist_adjacent = next_ex[0] - prev_ex[1]   # intron length
        elif prev_ex is not None:
            dist_adjacent = 0
        elif next_ex is not None:
            dist_adjacent = 0

        return 'intronic', dist_closest, dist_adjacent

    # 3. Intergenic
    up_gene, down_gene = flanking_genes(pos, chrom, genes, gene_starts, gene_ends)

    dist_up   = (pos - up_gene[1])   if up_gene   else None
    dist_down = (down_gene[0] - pos) if down_gene else None

    # Determine upstream / downstream by distance
    if dist_up is None and dist_down is None:
        return 'intergenic_upstream', '.', '.'

    if dist_up is None:
        cls = 'intergenic_downstream'
        dist_closest  = dist_down
        dist_adjacent = 0
    elif dist_down is None:
        cls = 'intergenic_upstream'
        dist_closest  = dist_up
        dist_adjacent = 0
    else:
        # Closer to upstream gene → downstream of it; closer to downstream → upstream of it
        if dist_up <= dist_down:
            cls = 'intergenic_downstream'
            dist_closest = dist_up
        else:
            cls = 'intergenic_upstream'
            dist_closest = dist_down
        dist_adjacent = down_gene[0] - up_gene[1]

    return cls, dist_closest, dist_adjacent


# ── VCF parser ────────────────────────────────────────────────────────

def parse_vcf(vcf_path):
    """
    Yield dicts with keys: chrom, pos(0-based int), end, id, ref, alt, filter
    Handles multi-allelic ALT by joining with comma.
    """
    with open_file(vcf_path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 7:
                continue
            chrom  = parts[0]
            pos1   = int(parts[1])          # 1-based VCF
            pos0   = pos1 - 1               # convert to 0-based
            vid    = parts[2] if parts[2] != '.' else f"{chrom}:{pos1}"
            ref    = parts[3]
            alt    = parts[4]
            filt   = parts[6] if len(parts) > 6 else '.'
            end    = pos0 + len(ref)        # half-open end
            yield {
                'chrom': chrom, 'pos0': pos0, 'end': end,
                'id': vid, 'ref': ref, 'alt': alt, 'filter': filt
            }


# ── main ──────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description='Classify VCF variants relative to GTF gene features.')
    ap.add_argument('--vcf',  required=True,  help='Input VCF or VCF.gz')
    ap.add_argument('--gtf',  required=True,  help='GTF or GTF.gz annotation')
    ap.add_argument('--out',  required=True,  help='Output BED-like TSV file')
    ap.add_argument('--upstream-bp',   type=int, default=2000,
                    help='Max bp upstream to label as intergenic_upstream (default: 2000)')
    ap.add_argument('--downstream-bp', type=int, default=2000,
                    help='Max bp downstream to label as intergenic_downstream (default: 2000)')
    args = ap.parse_args()

    # Load annotation
    genes, exons, gene_starts, gene_ends, exon_starts, exon_ends = load_gtf(args.gtf)

    # Process variants
    header = [
        '#chrom', 'start', 'end', 'ID', 'REF', 'ALT', 'FILTER',
        'class', 'dist_closest', 'dist_adjacent'
    ]

    n_written = 0
    print(f"Classifying variants from: {args.vcf}", file=sys.stderr)

    with open(args.out, 'w', newline='') as out_fh:
        writer = csv.writer(out_fh, delimiter='\t')
        writer.writerow(header)

        for var in parse_vcf(args.vcf):
            cls, dist_c, dist_a = classify_variant(
                var['pos0'], var['chrom'],
                genes, exons,
                gene_starts, gene_ends,
                exon_starts, exon_ends,
                upstream_bp=args.upstream_bp,
                downstream_bp=args.downstream_bp
            )
            writer.writerow([
                var['chrom'],
                var['pos0'],
                var['end'],
                var['id'],
                var['ref'],
                var['alt'],
                var['filter'],
                cls,
                dist_c,
                dist_a
            ])
            n_written += 1
            if n_written % 50_000 == 0:
                print(f"  Processed {n_written:,} variants...", file=sys.stderr)

    print(f"Done. Wrote {n_written:,} variants to {args.out}", file=sys.stderr)


if __name__ == '__main__':
    main()
