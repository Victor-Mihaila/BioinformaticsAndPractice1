#!/usr/bin/env python3
# backprob_by_namespace_v3.py
#
# Compute background GO-term probabilities per namespace.
# Works on:
#   • GAF      (--fmt gaf)  – standard 17-column GO annotation files
#   • 2-col TSV (--fmt tsv) – <identifier>\t<GO_id> rows

import argparse
from collections import Counter, defaultdict
from pathlib import Path
from goatools.obo_parser import GODag


def parse_args():
    p = argparse.ArgumentParser(
        description="Background GO-term probabilities (roots = 1.0)"
    )
    p.add_argument("--infile", required=True, type=Path,
                   help="Input annotation file (GAF or 2-column TSV)")
    p.add_argument("--fmt", choices=["gaf", "tsv"], default="gaf",
                   help="Input file format: gaf (default) or tsv")
    p.add_argument("--obo", required=True, type=Path,
                   help="GO ontology in OBO format")
    p.add_argument("--out", required=True, type=Path,
                   help="Output file prefix (writes *_BP.tsv, *_MF.tsv, *_CC.tsv)")
    return p.parse_args()


def ancestors_plus_self(go_id, godag):
    """Return {go_id} ∪ ancestors ; empty set if unknown/obsolete."""
    if go_id not in godag:
        return set()
    return {go_id} | godag[go_id].get_all_parents()


def iterate_go_ids(path: Path, file_fmt: str):
    """Yield GO IDs, one per input line, according to file format."""
    with path.open() as fh:
        for ln in fh:
            if not ln.strip():
                continue
            if file_fmt == "gaf":
                if ln.startswith("!"):
                    continue
                cols = ln.rstrip("\n").split("\t")
                if len(cols) < 5:
                    continue
                yield cols[4]                         # column 5 = GO ID
            else:  # 2-col TSV
                go_id = ln.rstrip("\n").split("\t")[1]
                yield go_id


def main():
    args = parse_args()
    godag = GODag(str(args.obo), optional_attrs={'relationship'})

    tally = defaultdict(Counter)    # tally[ns][GO] = hits incl. ancestors
    lines_per_ns = Counter()        # denominator per namespace

    for go_id in iterate_go_ids(args.infile, args.fmt):
        if go_id not in godag:
            continue                       # skip obsolete / unknown IDs
        ns = godag[go_id].namespace        # "biological_process", …
        lines_per_ns[ns] += 1

        for anc in ancestors_plus_self(go_id, godag):
            tally[ns][anc] += 1

    ns2tag = {"biological_process": "BP",
              "molecular_function": "MF",
              "cellular_component": "CC"}

    for ns, counter in tally.items():
        denom = lines_per_ns[ns]
        if denom == 0:
            continue

        outfile = Path(f'{args.out}_{ns2tag[ns]}.tsv')
        print(f"Writing {outfile} (normalised by {denom} lines)")

        with outfile.open("w") as out_fh:
            for go_id, cnt in sorted(counter.items(),
                                     key=lambda kv: kv[1],
                                     reverse=True):
                out_fh.write(f"{go_id}\t{cnt/denom:.8f}\n")


if __name__ == "__main__":
    main()
