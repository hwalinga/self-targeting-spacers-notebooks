#!/usr/bin/env python3
"""
Remove all genomes that have the same hit patterns
(thus indicating a double genome)
Also add information about spacer position.
"""


import sys
from itertools import islice

import pandas as pd

df = pd.read_csv(
    sys.stdin,
    sep="\t",
    names=[
        "genome_key",
        "genome_id",
        "contig_hit",
        "c_1",
        "c_2",
        "spacer_id",
        "hit_ident",
        "spacer_size",
        "contig_size",
        "3_prime",
        "5_prime",
        "array_id",
        "array_size",
        "array_confid",
        "repeat_size",
        "array_type",
        "genome_type",
        "PAM_confid",
        "poss_PAM",
        "PAM_side",
        "phage_class",
        "with_phage",
        "gene_hit",
        "gene_id",
        "gene_class",
    ],
    index_col=False,
    keep_default_na=False,
    dtype={"genome_id": str},
).set_index("genome_id")

double_genome = df.groupby("spacer_id").apply(
    lambda x: "-".join((islice(set(x.index), 1, None)))
)
remove_genome_id = set(double_genome[double_genome != ""].unique())
df = df[~df.index.isin(remove_genome_id)]

# By the definition of the spacer ID, the last number indicates the position in the array
# The orientation of the array is defined by CRISPRDetect

df = df.assign(spacer_pos=lambda x: x["spacer_id"].str.split("_").str[-1].astype(int))

# relative spacer position
# to change a number n range a-b with fixed a and variable b,
# to a relative position 0-100 use
# (n-a)/(b-a)
# 0: leader
# 1: tail
df["rel_spacer_pos"] = df["spacer_pos"].subtract(1).divide(df["array_size"].subtract(1))

df["name_spacer_pos"] = df["rel_spacer_pos"].apply(
    lambda x: {0: "leader", 1: "tail"}.get(x) or "middle"
)

df.to_csv(sys.stdout, sep="\t")
