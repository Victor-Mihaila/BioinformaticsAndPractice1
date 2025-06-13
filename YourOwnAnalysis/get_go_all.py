
import pandas as pd
import pickle
import numpy as np
import plotly.express as px
from GO_similarity import parse_background_probs, sim_go_families, sim_go_term_family
from goatools.obo_parser import GODag
import matplotlib.pyplot as plt
import re
import gseapy as gp
cnts = pd.read_csv('/home/victor/bioinfo1/binfo1-work/read-counts.txt', sep='\t', comment='#')
cnts['Geneid'] = cnts['Geneid'].str.split('.').str[0]
cnts['clip_enrichment'] = (cnts['CLIP-35L33G.bam'] / cnts['RNA-control.bam'])
cnts['ribosome_density_change'] = (cnts['RPF-siLin28a.bam'] / cnts['RNA-siLin28a.bam'] / (cnts['RPF-siLuc.bam'] / cnts['RNA-siLuc.bam']))
cnts.dropna(axis=0, inplace=True)
cnts.set_index("Geneid", inplace=True, drop=True)
cnts = cnts.replace([np.inf, -np.inf], np.nan).dropna(subset=["clip_enrichment"])
cnt_ids = cnts.index.to_series().tolist() 
print(len(cnt_ids))
# cnts.head()
# cnts['clip_enrichment'].hist(range=[0,15])
# cnts['ribosome_density_change'].hist(range=[0,15])

def batched(iterable, n): 
    items = []
    for item in iterable:
        items.append(item)
        if len(items) == n:
            yield items
            items = []
    if items:
        yield items
        
bm = gp.biomart.Biomart()
results = []

print("gettig go")
for batch in batched(cnt_ids, 100):
    res = bm.query(dataset='mmusculus_gene_ensembl', #7
                   attributes=['ensembl_gene_id', 'external_gene_name', 'go_id'], #8
                   filters={'ensembl_gene_id': batch})
    results.append(res)

print("finished getting go")
results = pd.concat(results)
results.to_csv("go.csv")
results = results.replace("", np.nan).dropna(subset=["go_id"])
gene2go = (
    results.groupby("ensembl_gene_id")
           .agg({
               "external_gene_name": "first",
               "go_id": lambda g: list(set(g))     # unique GO IDs
           })
           .rename(columns={
               "external_gene_name": "gene_symbol",
               "go_id": "GO_list"
           })
           .reset_index()
)

print(gene2go.head())

cnts = cnts.reset_index()                     
cnts = cnts.merge(
    gene2go,
    left_on = "Geneid",
    right_on = "ensembl_gene_id",
    how = "left"
)

# if you want GO as a single semicolon-separated string instead of list:
cnts["GO"] = cnts["GO_list"].apply(
    lambda lst: ";".join(lst) if isinstance(lst, list) else ""
)

cnts.drop(columns=["ensembl_gene_id"], inplace=True)
cnts.set_index("Geneid", inplace=True)

cnts.to_csv("counts_and_GO.csv")