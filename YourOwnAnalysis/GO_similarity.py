#!/usr/bin/env python

import sys
import os
import math
from goatools.obo_parser import GODag

def parse_background_probs(prob_file):
    """
    Reads a tab- or space-separated file where each line is:
        <GO_term>  <probability>
    Returns a dictionary { go_id: float_probability }.
    Skips lines that don't parse correctly.
    """
    bg_probs = {}
    with open(prob_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            go_id, prob_str = parts[0], parts[1]
            try:
                bg_probs[go_id] = float(prob_str)
            except ValueError:
                pass  # skip lines with invalid float
    return bg_probs

def sim_go_terms(a, b, go_dag, background_probs):
    """
    Computes the similarity score (based on your chosen formula):
        sim(a,b) = 2 * [ ln( P(LCA(a,b)) ) / ln( P(a)*P(b) ) ]

    or any other variant you prefer. Returns the float score or None if missing data.
    """
    if a not in background_probs or b not in background_probs:
        return None  # missing probabilities
    if a not in go_dag or b not in go_dag:
        return None  # missing in DAG

    ancestors_a = go_dag[a].get_all_parents()
    ancestors_a.add(a)
    ancestors_b = go_dag[b].get_all_parents()
    ancestors_b.add(b)

    common_ancestors = ancestors_a.intersection(ancestors_b)
    if not common_ancestors:
        return None

    # Pick the LCA with greatest depth
    lca = max(common_ancestors, key=lambda go_id: go_dag[go_id].depth)
    if lca not in background_probs:
        return None

    p_a = background_probs[a]
    p_b = background_probs[b]
    p_lca = background_probs[lca]

    if p_a <= 0 or p_b <= 0 or p_lca <= 0:
        return None

    num = math.log(p_lca)
    denom = math.log(p_a * p_b)
    if denom == 0:
        return None

    ratio = num / denom
    return 2.0 * ratio

def sim_go_term_family(a, family_go_terms, go_dag, background_probs):
    """
    Returns sim(a, G) = max( sim(a,g) for g in G ), or None if none are computable.
    """
    best_score = None
    for g in family_go_terms:
        score = sim_go_terms(a, g, go_dag, background_probs)
        if score is not None:
            if best_score is None or score > best_score:
                best_score = score
    return best_score

def sim_go_families(A, B, go_dag, background_probs):
    """
    Computes the similarity between two families A, B:
      sim(A,B) = ( sum_{a in A}[sim(a,B)] + sum_{b in B}[sim(b,A)] ) / (|A| + |B|)

    where sim(a,B) = max_{g in B}[ sim(a,g) ] and sim(b,A) = max_{g in A}[ sim(b,g) ].
    Returns the average, or None if A and B are both empty.
    """
    # If both sets/families are empty, undefined
    if not A and not B:
        return None

    sumA = 0.0
    for a in A:
        s = sim_go_term_family(a, B, go_dag, background_probs)
        # If sim_go_term_family returns None, skip or treat as 0
        if s is not None:
            sumA += s

    sumB = 0.0
    for b in B:
        s = sim_go_term_family(b, A, go_dag, background_probs)
        if s is not None:
            sumB += s

    denom = len(A) + len(B)
    if denom == 0:
        return None  # avoid division by zero if one is empty but the other is not

    return (sumA + sumB) / denom

def main():
    # Example usage:
    #   python GO_similarity.py go-basic.obo mf_backprobs.tsv "GO:0003674,GO:0005488" "GO:0016491"
    #
    # Then interpret each quoted argument as a comma-separated list of GO terms.
    # We'll parse them into two families, A and B, then compute sim(A,B).

    if len(sys.argv) != 4:
        script_name = os.path.basename(__file__)
        print(f"Usage: python {script_name} <obo_file> <prob_file> <A_list> <B_list>")
        print("Where <A_list> and <B_list> are comma-separated GO terms.")
        sys.exit(1)

    obo_file = sys.argv[1]
    prob_file = sys.argv[2]
    familyA_str = sys.argv[3]
    familyB_str = sys.argv[4]

    # Parse families
    familyA = [t.strip() for t in familyA_str.split(",") if t.strip()]
    familyB = [t.strip() for t in familyB_str.split(",") if t.strip()]

    # Load DAG
    print(f"Loading GO DAG from: {obo_file}")
    go_dag = GODag(obo_file)

    # Load background probabilities
    print(f"Parsing background probabilities from: {prob_file}")
    background_probs = parse_background_probs(prob_file)
    if not background_probs:
        print("No valid probabilities found.")
        sys.exit(1)

    # Compute family-familiy similarity
    score = sim_go_families(familyA, familyB, go_dag, background_probs)
    if score is None:
        print(f"Could not compute sim(A,B). Possibly no overlap or missing data.")
    else:
        print(f"sim(A,B) = {score:.4f}")

if __name__ == "__main__":
    main()
