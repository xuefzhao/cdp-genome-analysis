#!/usr/bin/env python3
"""
rename_bam_chroms.py
────────────────────
Convert BAM chromosome names from Ensembl style (1, 2, KI270419.1 …)
to UCSC hg38 style (chr1, chr2, chrUn_KI270419v1 …) using the
official UCSC chromAlias table.

Steps
-----
1. Download (or use a local copy of) the UCSC chromAlias.txt.gz for hg38
2. Build an Ensembl → UCSC name mapping
3. Rewrite the BAM @SQ header lines with the new names
4. Call `samtools reheader` to produce the output BAM

Requirements
------------
  - samtools ≥ 1.10  (in PATH)
  - Python ≥ 3.8 (stdlib only)

Usage
-----
  python rename_bam_chroms.py \\
      --bam       input.bam \\
      --out       output.bam \\
      [--alias    chromAlias.txt.gz]   # downloaded automatically if omitted
      [--unmapped keep|drop]           # what to do with contigs not in alias table
                                       # keep = leave name unchanged (default)
                                       # drop = remove @SQ line (reads on that contig lost)
"""

import argparse
import gzip
import os
import subprocess
import sys
import tempfile
import urllib.request

ALIAS_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/chromAlias.txt.gz"


# ── alias table ──────────────────────────────────────────────────────

def download_alias(dest_path):
    print(f"Downloading chromAlias table from UCSC → {dest_path}", file=sys.stderr)
    urllib.request.urlretrieve(ALIAS_URL, dest_path)
    print("  Download complete.", file=sys.stderr)


def load_alias(alias_path):
    """
    Parse chromAlias.txt(.gz).

    UCSC chromAlias files come in two formats depending on release:

    Format A – older (3 columns):
        alias           source      ucscName
        1               ensembl     chr1
        KI270419.1      ensembl     chrUn_KI270419v1

    Format B – newer (columns may be reordered, ucsc name is first):
        #ucscName       alias       source
        chr1            1           ensembl
        chrUn_KI270419v1 KI270419.1 ensembl

    Strategy: for each data row, identify which field looks like a UCSC name
    (starts with 'chr') and which looks like an Ensembl name (does not start
    with 'chr'). Build mapping: ensembl_name → ucsc_name regardless of column
    order.  Rows where both or neither field starts with 'chr' are skipped.
    """
    mapping = {}
    opener = gzip.open if alias_path.endswith('.gz') else open
    skipped = 0

    with opener(alias_path, 'rt') as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                continue

            # Find the UCSC-style name (starts with 'chr') among all fields
            ucsc_candidates = [p for p in parts if p.startswith('chr')]
            # Find the non-UCSC names (Ensembl / RefSeq / etc.)
            other_candidates = [p for p in parts if not p.startswith('chr')
                                 and p not in ('ensembl', 'refseq', 'assembly',
                                               'ucsc', 'genbank')]

            if len(ucsc_candidates) != 1 or not other_candidates:
                skipped += 1
                continue

            ucsc_name = ucsc_candidates[0]
            for alias in other_candidates:
                if alias != ucsc_name:
                    mapping[alias] = ucsc_name

    if skipped:
        print(f"  Skipped {skipped} ambiguous alias rows.", file=sys.stderr)
    print(f"  Loaded {len(mapping):,} alias→UCSC name mappings.", file=sys.stderr)
    if len(mapping) == 0:
        print("  WARNING: 0 mappings loaded — check alias file format.",
              file=sys.stderr)
        print("  First few raw lines:", file=sys.stderr)
        with opener(alias_path, 'rt') as fh:
            for i, line in enumerate(fh):
                print(f"    {repr(line.rstrip())}", file=sys.stderr)
                if i >= 5:
                    break
    return mapping


# ── header rewriting ──────────────────────────────────────────────────

def rewrite_header(old_header_lines, mapping, unmapped='keep'):
    """
    Rewrite @SQ SN: fields using mapping.

    Parameters
    ----------
    old_header_lines : list[str]  (each line already stripped of newline)
    mapping          : dict  ensembl → ucsc
    unmapped         : 'keep' | 'drop'

    Returns list[str] of new header lines.
    Raises ValueError if a contig is encountered and unmapped='drop'
    (just skips the @SQ line; reads mapped there will be unmapped in output).
    """
    new_lines = []
    renamed, kept, dropped = 0, 0, 0

    for line in old_header_lines:
        if not line.startswith('@SQ'):
            new_lines.append(line)
            continue

        # Parse SN field
        fields = line.split('\t')
        sn_idx = next((i for i, f in enumerate(fields) if f.startswith('SN:')), None)
        if sn_idx is None:
            new_lines.append(line)
            continue

        old_name = fields[sn_idx][3:]   # strip 'SN:'

        if old_name in mapping:
            fields[sn_idx] = 'SN:' + mapping[old_name]
            new_lines.append('\t'.join(fields))
            renamed += 1
        elif unmapped == 'drop':
            # silently drop this @SQ line
            dropped += 1
        else:
            # keep as-is
            new_lines.append(line)
            kept += 1

    print(f"  Header: {renamed} contigs renamed, {kept} kept unchanged, "
          f"{dropped} dropped.", file=sys.stderr)
    return new_lines


# ── main ──────────────────────────────────────────────────────────────

def run(cmd, check=True):
    print(f"  $ {' '.join(cmd)}", file=sys.stderr)
    return subprocess.run(cmd, check=check)


def main():
    ap = argparse.ArgumentParser(
        description='Rename BAM chromosomes from Ensembl to UCSC hg38 style.')
    ap.add_argument('--bam',      required=True,
                    help='Input BAM (Ensembl contig names)')
    ap.add_argument('--out',      required=True,
                    help='Output BAM (UCSC contig names)')
    ap.add_argument('--alias',    default=None,
                    help='Path to chromAlias.txt or chromAlias.txt.gz '
                         '(downloaded automatically if not provided)')
    ap.add_argument('--unmapped', choices=['keep', 'drop'], default='keep',
                    help='What to do with contigs absent from alias table '
                         '(default: keep — leave name unchanged)')
    ap.add_argument('--ref',      default=None,
                    help='Reference FASTA — required when input is a CRAM file '
                         'so samtools can decode reads during reheader')
    args = ap.parse_args()

    # ── Step 1: get alias table ────────────────────────────────────────
    alias_path = args.alias
    tmp_alias  = None
    if alias_path is None:
        tmp_alias  = tempfile.NamedTemporaryFile(suffix='.chromAlias.txt.gz',
                                                 delete=False)
        alias_path = tmp_alias.name
        tmp_alias.close()
        download_alias(alias_path)

    print("Loading alias table…", file=sys.stderr)
    mapping = load_alias(alias_path)

    if tmp_alias:
        os.unlink(alias_path)

    # ── Step 2: extract current header ────────────────────────────────
    print("Extracting BAM header…", file=sys.stderr)
    result = subprocess.run(
        ['samtools', 'view', '-H', args.bam],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        universal_newlines=True, check=True
    )
    old_lines = result.stdout.rstrip('\n').split('\n')

    # ── Step 3: rewrite header ─────────────────────────────────────────
    print("Rewriting contig names in header…", file=sys.stderr)
    # Debug: show first few contig names before renaming
    sq_lines = [l for l in old_lines if l.startswith('@SQ')]
    sample_names = []
    for l in sq_lines[:5]:
        for f in l.split('\t'):
            if f.startswith('SN:'):
                sample_names.append(f[3:])
    print(f"  First contig names in header: {sample_names}", file=sys.stderr)
    new_lines = rewrite_header(old_lines, mapping, unmapped=args.unmapped)

    # ── Step 4: write new header to temp file ─────────────────────────
    with tempfile.NamedTemporaryFile(mode='w', suffix='.sam',
                                     delete=False) as tmp_hdr:
        tmp_hdr.write('\n'.join(new_lines) + '\n')
        tmp_hdr_path = tmp_hdr.name

    # ── Step 5: samtools reheader ──────────────────────────────────────
    # samtools reheader always writes to stdout (no --output flag).
    # Redirect stdout directly to the output file.
    print("Running samtools reheader…", file=sys.stderr)
    cmd = ['samtools', 'reheader', tmp_hdr_path, args.bam]
    # For CRAM input, samtools needs the reference to decode reads
    if args.ref:
        cmd = ['samtools', 'reheader', '-T', args.ref, tmp_hdr_path, args.bam]
    print("  $ " + ' '.join(cmd) + " > " + args.out, file=sys.stderr)
    try:
        with open(args.out, 'wb') as out_fh:
            subprocess.run(cmd, stdout=out_fh, check=True)
    finally:
        os.unlink(tmp_hdr_path)

    print(f"\nDone. Output: {args.out}", file=sys.stderr)


if __name__ == '__main__':
    main()
