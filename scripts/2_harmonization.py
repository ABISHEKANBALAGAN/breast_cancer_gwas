#!/usr/bin/env python3
import pandas as pd

# Load eQTL data
eqtl = pd.read_csv(
    "data/Breast_Mammary_Tissue.v8.egenes.txt",
    sep="\t",
    usecols=["variant_id", "chr", "variant_pos", "ref", "alt", "maf", "rs_id_dbSNP151_GRCh38p7"]
)

# Clean chromosome column (remove 'chr' prefix and spaces)
eqtl["chr"] = eqtl["chr"].str.replace("chr", "").str.strip()  # Remove 'chr' and any whitespace

# Filter to keep only chromosomes 1-22 (exclude X, Y, MT etc.)
eqtl = eqtl[eqtl["chr"].isin([str(x) for x in range(1, 23)])]
eqtl["chr"] = eqtl["chr"].astype(int)  # Now safe to convert to integer
eqtl["variant_pos"] = eqtl["variant_pos"].astype(int)

# Process GWAS in chunks
for i, gwas in enumerate(pd.read_csv("results/cleaned_gwas.tsv", sep="\t", chunksize=500000)):
    # Ensure numeric types (GWAS should already be 1-22 per your note)
    gwas["CHR"] = gwas["CHR"].astype(int)
    gwas["POS"] = gwas["POS"].astype(int)
    
    # Merge on chromosome and position
    merged = pd.merge(
        gwas,
        eqtl,
        left_on=["CHR", "POS"],
        right_on=["chr", "variant_pos"],
        how="inner"
    )
    
    # Relaxed allele matching
    merged = merged[
        (merged["A1"].str.upper() == merged["alt"].str.upper()) |
        (merged["A1"].str.upper() == merged["ref"].str.upper())
    ]
    
    # Flip effect sizes where needed
    merged.loc[merged["A1"].str.upper() == merged["alt"].str.upper(), "BETA"] *= -1
    
    # Rename columns for clarity
    merged = merged.rename(columns={
        "maf": "MAF_tissue",
        "rs_id_dbSNP151_GRCh38p7": "rsID"
    })
    
    # Save results
    merged.to_csv(
    "results/harmonized.tsv",
    sep="\t",
    mode="w" if i == 0 else "a",
    header=(i == 0),
    index=False
)
    print(f"Processed chunk {i+1}: {len(merged)} variants matched")

print("Harmonization complete!")