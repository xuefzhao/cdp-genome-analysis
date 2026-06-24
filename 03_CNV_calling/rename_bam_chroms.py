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

    UCSC format (3 columns, tab-separated):
        alias           source      ucscName
        1               ensembl     chr1
        KI270419.1      ensembl     chrUn_KI270419v1
        chrUn_KI270419v1 ucsc       chrUn_KI270419v1   ← identity rows exist too

    Returns dict: ensembl_name → ucsc_name
    Only rows whose source contains 'ensembl' are used; identity rows are skipped.
    """
    mapping = {}
    opener = gzip.open if alias_path.endswith('.gz') else open

    with opener(alias_path, 'rt') as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 3:
                continue
            alias, source, ucsc_name = parts[0], parts[1], parts[2]
            if 'ensembl' in source.lower() and alias != ucsc_name:
                mapping[alias] = ucsc_name

    print(f"  Loaded {len(mapping):,} Ensembl→UCSC name mappings.", file=sys.stderr)
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
        capture_output=True, text=True, check=True
    )
    old_lines = result.stdout.rstrip('\n').split('\n')

    # ── Step 3: rewrite header ─────────────────────────────────────────
    print("Rewriting contig names in header…", file=sys.stderr)
    new_lines = rewrite_header(old_lines, mapping, unmapped=args.unmapped)

    # ── Step 4: write new header to temp file ─────────────────────────
    with tempfile.NamedTemporaryFile(mode='w', suffix='.sam',
                                     delete=False) as tmp_hdr:
        tmp_hdr.write('\n'.join(new_lines) + '\n')
        tmp_hdr_path = tmp_hdr.name

    # ── Step 5: samtools reheader ──────────────────────────────────────
    print("Running samtools reheader…", file=sys.stderr)
    try:
        run(['samtools', 'reheader', tmp_hdr_path, args.bam,
             '--output', args.out])
    except TypeError:
        # older samtools: reheader writes to stdout, no --output flag
        print("  (falling back to stdout redirect for older samtools)",
              file=sys.stderr)
        with open(args.out, 'wb') as out_fh:
            subprocess.run(
                ['samtools', 'reheader', tmp_hdr_path, args.bam],
                stdout=out_fh, check=True
            )
    finally:
        os.unlink(tmp_hdr_path)

    # ── Step 6: index ──────────────────────────────────────────────────
    print("Indexing output BAM…", file=sys.stderr)
    run(['samtools', 'index', args.out])

    print(f"\nDone. Output: {args.out}", file=sys.stderr)


if __name__ == '__main__':
    main()
