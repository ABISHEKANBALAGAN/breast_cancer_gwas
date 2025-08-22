#!/usr/bin/env python3
# scripts/1_preprocess.py
import pandas as pd

CHUNKSIZE = 500000  # Process 500K lines at a time
COLS = {
    'chromosome': 'CHR',
    'base_pair_location': 'POS',
    'effect_allele': 'A1',
    'other_allele': 'A2',
    'beta': 'BETA',
    'standard_error': 'SE',
    'p_value': 'P',
    'variant_id': 'SNP',
    'effect_allele_frequency': 'MAF'
}

# Process in chunks
for i, chunk in enumerate(pd.read_csv("data/GCST90435600.tsv", sep="\t", chunksize=CHUNKSIZE)):
    chunk = chunk[COLS.keys()].rename(columns=COLS)
    chunk.to_csv(
        "results/cleaned_gwas.tsv",
        sep="\t",
        mode="w" if i == 0 else "a",
        header=(i == 0),
        index=False
    )
    print(f"Processed chunk {i+1}: {len(chunk)} variants")