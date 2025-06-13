#!/usr/bin/env python3
# backprob_by_namespace_v2.py
#
# Background probabilities for GO terms, *per namespace*, normalised so
# each root (GO:0008150/3674/5575) gets probability 1.

import argparse
from collections import Counter, defaultdict
from pathlib import Path
from goatools.obo_parser import GODag


def parse_args():
    ap = argparse.ArgumentParser(
        description="Background GO-term probabilities (roots = 1.0)"
    )
    ap.add_argument("--gaf", required=True, type=Path, help="Input GAF file")
    ap.add_argument("--obo", required=True, type=Path, help="GO ontology in OBO format")
    ap.add_argument("--out", required=True, type=Path,
                    help="Output file prefix (three TSVs will be written)")
    return ap.parse_args()


def ancestors_including_self(go_id, godag):
    """Return {term itself} U {all its ancestors}.  Empty if term unknown."""
    if go_id not in godag:
        return set()
    term = godag[go_id]
    return {go_id} | term.get_all_parents()


def main():
    args = parse_args()

    godag = GODag(str(args.obo), optional_attrs={'relationship'})

    # tally[namespace][go_id]  ,  lines_per_ns[namespace] = number of *original* lines
    tally = defaultdict(Counter)
    lines_per_ns = Counter()

    with args.gaf.open() as gaf_fh:
        for line in gaf_fh:
            if not line.strip() or line.startswith("!"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 5:
                continue
            go_id = cols[4]

            if go_id not in godag:          # skip obsolete / unknown
                continue
            ns = godag[go_id].namespace     # 'biological_process', ...
            lines_per_ns[ns] += 1           # 1 unit for the denominator

            for anc in ancestors_including_self(go_id, godag):
                tally[ns][anc] += 1         # ancestor counts

    ns2tag = {"biological_process": "BP",
              "molecular_function": "MF",
              "cellular_component": "CC"}

    for ns, counter in tally.items():
        tot_lines = lines_per_ns[ns]
        if tot_lines == 0:
            continue

        outfile = Path(f'{args.out}_{ns2tag[ns]}.tsv')
        print(f"Writing {outfile}  (normalised by {tot_lines} lines)")

        with outfile.open("w") as out_fh:
            for go_id, cnt in sorted(counter.items(),
                                     key=lambda kv: kv[1],
                                     reverse=True):
                prob = cnt / tot_lines
                out_fh.write(f"{go_id}\t{prob:.8f}\n")


if __name__ == "__main__":
    main()
