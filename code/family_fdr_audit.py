import numpy as np
import pandas as pd


# ============================================================
# BENJAMINI-HOCHBERG FDR
# ============================================================

def bh_fdr(pvalues):

    pvalues = np.asarray(pvalues, dtype=float)

    n = len(pvalues)

    order = np.argsort(pvalues)
    ranked = pvalues[order]

    q = ranked * n / np.arange(1, n + 1)

    # enforce monotonicity
    q = np.minimum.accumulate(q[::-1])[::-1]

    result = np.empty(n)
    result[order] = q

    return result


# ============================================================
# TEST FAMILIES
# ============================================================

tests = [

    # --------------------------------------------------------
    # PROXIMITY / PRIMARY
    # --------------------------------------------------------

    ("Proximity",
     "Whole-UTR permutation",
     0.0001999600079984003),

    ("Proximity",
     "CDS-overlap excluded",
     0.0001999600079984003),

    ("Proximity",
     "ClinVar observed-space",
     0.41711657668466307),

    ("Proximity",
     "High-confidence P/LP whole-UTR",
     0.0001999600079984003),

    ("Proximity",
     "High-confidence P/LP observed-space",
     0.2403519296140772),

    # --------------------------------------------------------
    # REVIEW STATUS / CURATION
    # --------------------------------------------------------

    ("Review-status",
     "Within-group review-status permutation",
     0.0001999600079984003),

    ("Review-status",
     "Coordinate-level review-status permutation",
     0.0001999600079984003),


    # --------------------------------------------------------
    # REGRESSION
    # --------------------------------------------------------

    ("Regression",
     "Original regression distance",
     7.1787148842081996e-33),

    ("Regression",
     "Original regression UTR length",
     0.012184198744742206),

    ("Regression",
     "Original regression UTR type",
     0.018973053720891114),

    # --------------------------------------------------------
    # GENE LEVEL
    # --------------------------------------------------------

    ("Gene-level",
     "Gene-burden distance",
     1.35533e-15),

    ("Gene-level",
     "Gene-burden UTR length",
     0.00102394),

    ("Gene-level",
     "Gene-burden UTR type",
     1.70288e-07),

    ("Gene-level",
     "Gene-burden HC P/LP burden",
     3.98485e-08),

    # --------------------------------------------------------
    # POPULATION / AF
    # --------------------------------------------------------

    ("Population",
     "WES-AF distance",
     1.93983e-06),

    ("Population",
     "WES-AF UTR type",
     0.00345128),

    ("Population",
     "WES-AF allele frequency",
     0.00090023),

    # --------------------------------------------------------
    # CONSERVATION
    # --------------------------------------------------------

    ("Conservation",
     "phyloP Single vs Multiple",
     0.1967606478704259),
]


# ============================================================
# DATAFRAME
# ============================================================

df = pd.DataFrame(
    tests,
    columns=[
        "Family",
        "Test",
        "P_value"
    ]
)


# ============================================================
# FAMILY-SPECIFIC BH
# ============================================================

df["FDR_q"] = np.nan

for family in df["Family"].unique():

    mask = df["Family"] == family

    df.loc[mask, "FDR_q"] = bh_fdr(
        df.loc[mask, "P_value"].values
    )


df["Significant_FDR_0.05"] = (
    df["FDR_q"] < 0.05
)


# ============================================================
# SORT
# ============================================================

df = df.sort_values(
    ["Family", "P_value"]
).reset_index(drop=True)


# ============================================================
# PRINT
# ============================================================

print("=" * 100)
print("FAMILY-SPECIFIC FDR AUDIT")
print("=" * 100)

for family in df["Family"].unique():

    print("\n" + "-" * 100)
    print(f"FAMILY: {family}")
    print("-" * 100)

    sub = df[df["Family"] == family]

    print(
        sub[
            [
                "Test",
                "P_value",
                "FDR_q",
                "Significant_FDR_0.05"
            ]
        ].to_string(index=False)
    )


# ============================================================
# SUMMARY
# ============================================================

print("\n" + "=" * 100)
print("SUMMARY")
print("=" * 100)

print(
    "Total tests:",
    len(df)
)

print(
    "FDR significant:",
    int(df["Significant_FDR_0.05"].sum())
)

print(
    "FDR non-significant:",
    int((~df["Significant_FDR_0.05"]).sum())
)


# ============================================================
# SAVE
# ============================================================

out = r"..\data\family_specific_fdr_audit.tsv"

df.to_csv(
    out,
    sep="\t",
    index=False
)

print("\nSaved:")
print(out)
