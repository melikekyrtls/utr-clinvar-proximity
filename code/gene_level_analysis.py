import pandas as pd
import numpy as np
from scipy.stats import spearmanr, pearsonr

BASE = r"..\data"

CANONICAL = BASE + r"\canonical_utr_vus_dataset.tsv"
CLINVAR = BASE + r"\ClinVar_UTR_GRCh38.tsv"

PLP_VALUES = {
    "Pathogenic",
    "Likely pathogenic",
    "Pathogenic/Likely pathogenic",
    "Likely pathogenic/Pathogenic"
}

GROUP_COLS = ["Gene", "UTR_Type", "Transcript"]


# ============================================================
# LOAD
# ============================================================

vus = pd.read_csv(CANONICAL, sep="\t")
cv = pd.read_csv(CLINVAR, sep="\t", low_memory=False)


# ============================================================
# STANDARDIZE COORDINATES
# IMPORTANT: USE START, NOT PositionVCF
# ============================================================

cv["Start"] = pd.to_numeric(
    cv["Start"],
    errors="coerce"
)

cv = cv[
    cv["Start"].notna()
].copy()

cv["Position"] = cv["Start"].astype(int)


# ============================================================
# P/LP DEFINITION
# ============================================================

cv["ClinicalSignificance"] = (
    cv["ClinicalSignificance"]
    .astype(str)
    .str.strip()
)

cv["Is_PLP"] = cv[
    "ClinicalSignificance"
].isin(PLP_VALUES)


# ============================================================
# BASIC MAPPING
# ============================================================

cv = cv[
    cv["UTR_Transcript"].notna()
    &
    (cv["UTR_Transcript"].astype(str).str.lower() != "na")
    &
    cv["UTR_Type"].notna()
    &
    cv["GeneSymbol"].notna()
].copy()

cv["Gene"] = cv["GeneSymbol"].astype(str)
cv["Transcript"] = cv["UTR_Transcript"].astype(str)


# ============================================================
# CANONICAL VUS
# ============================================================

vus["Position"] = pd.to_numeric(
    vus["Position"],
    errors="coerce"
)

vus = vus.dropna(
    subset=["Position"]
).copy()

vus["Position"] = vus["Position"].astype(int)


# ============================================================
# UNIQUE CLINVAR POSITIONS
# Position-only deduplication
# ============================================================

cv_positions = (
    cv[
        GROUP_COLS + ["Position"]
    ]
    .drop_duplicates()
)

plp_positions = (
    cv[
        cv["Is_PLP"]
    ][
        GROUP_COLS + ["Position"]
    ]
    .drop_duplicates()
)


# ============================================================
# P/LP MAP BY GROUP
# ============================================================

plp_map = {}

for key, group in plp_positions.groupby(
    GROUP_COLS
):

    plp_map[key] = np.array(
        sorted(
            group["Position"]
            .astype(int)
            .unique()
        ),
        dtype=np.int64
    )


# ============================================================
# GROUP-LEVEL PROXIMITY
# ============================================================

group_results = []


for key, group in vus.groupby(GROUP_COLS):

    gene, utr_type, transcript = key

    positions = (
        group["Position"]
        .astype(int)
        .to_numpy()
    )

    unique_vus = np.unique(
        positions
    )

    plp = plp_map.get(key)

    if plp is None or len(plp) == 0:
        continue

    distances = []

    for pos in positions:

        diffs = np.abs(
            plp - pos
        )

        # Same-position variants are excluded,
        # matching the primary proximity analysis.
        diffs = diffs[
            diffs > 0
        ]

        if len(diffs) > 0:
            distances.append(
                int(diffs.min())
            )

    if not distances:
        continue

    # ClinVar burden within this exact group
    group_cv = cv_positions[
        (cv_positions["Gene"] == gene)
        &
        (cv_positions["UTR_Type"] == utr_type)
        &
        (cv_positions["Transcript"] == transcript)
    ]

    group_plp = plp_positions[
        (plp_positions["Gene"] == gene)
        &
        (plp_positions["UTR_Type"] == utr_type)
        &
        (plp_positions["Transcript"] == transcript)
    ]

    group_results.append({

        "Gene": gene,
        "UTR_Type": utr_type,
        "Transcript": transcript,

        "VUS_rows":
            len(positions),

        "VUS_unique_positions":
            len(unique_vus),

        "PLP_unique_positions":
            len(group_plp),

        "ClinVar_unique_positions":
            len(group_cv),

        "ClinVar_nonPLP_positions":
            len(group_cv) - len(group_plp),

        "Mean_nearest_PLP":
            np.mean(distances),

        "Median_nearest_PLP":
            np.median(distances),

        "Min_nearest_PLP":
            np.min(distances),

        "Max_nearest_PLP":
            np.max(distances)
    })


groups = pd.DataFrame(
    group_results
)


# ============================================================
# TRUE GENE-LEVEL AGGREGATION
# ============================================================

gene_results = []


for gene, g in groups.groupby("Gene"):

    total_vus_rows = g["VUS_rows"].sum()
    total_vus_unique = g["VUS_unique_positions"].sum()
    total_plp = g["PLP_unique_positions"].sum()
    total_clinvar = g["ClinVar_unique_positions"].sum()

    # Weighted mean proximity across VUS rows
    weighted_distance = (
        (
            g["Mean_nearest_PLP"]
            *
            g["VUS_rows"]
        ).sum()
        /
        total_vus_rows
    )

    gene_results.append({

        "Gene": gene,

        "Groups":
            len(g),

        "VUS_rows":
            int(total_vus_rows),

        "VUS_unique_positions":
            int(total_vus_unique),

        "PLP_unique_positions":
            int(total_plp),

        "ClinVar_unique_positions":
            int(total_clinvar),

        "ClinVar_nonPLP_positions":
            int(total_clinvar - total_plp),

        "VUS_fraction_of_ClinVar":
            total_vus_unique / total_clinvar
            if total_clinvar > 0
            else np.nan,

        "PLP_fraction_of_ClinVar":
            total_plp / total_clinvar
            if total_clinvar > 0
            else np.nan,

        "VUS_to_PLP_ratio":
            total_vus_unique / total_plp
            if total_plp > 0
            else np.nan,

        "Mean_nearest_PLP":
            weighted_distance,

        "Median_group_mean_distance":
            g["Mean_nearest_PLP"].median(),

        "Min_nearest_PLP":
            g["Min_nearest_PLP"].min(),

        "Max_nearest_PLP":
            g["Max_nearest_PLP"].max()
    })


df = pd.DataFrame(
    gene_results
)


# ============================================================
# SUMMARY
# ============================================================

print()
print("=" * 75)
print("TRUE GENE-LEVEL CLINVAR BIAS ANALYSIS")
print("=" * 75)

print()
print("Analyzed gene groups:", len(groups))
print("Analyzed genes:", len(df))

print(
    "Total VUS rows:",
    int(df["VUS_rows"].sum())
)

print(
    "Total VUS unique group-positions:",
    int(df["VUS_unique_positions"].sum())
)

print(
    "Total P/LP unique group-positions:",
    int(df["PLP_unique_positions"].sum())
)

print(
    "Total ClinVar unique group-positions:",
    int(df["ClinVar_unique_positions"].sum())
)


# ============================================================
# TOP GENES BY VUS BURDEN
# ============================================================

print()
print("=" * 75)
print("TOP 20 GENES BY VUS BURDEN")
print("=" * 75)

print(
    df[
        [
            "Gene",
            "Groups",
            "VUS_rows",
            "VUS_unique_positions",
            "PLP_unique_positions",
            "ClinVar_unique_positions",
            "Mean_nearest_PLP"
        ]
    ]
    .sort_values(
        "VUS_unique_positions",
        ascending=False
    )
    .head(20)
    .to_string(index=False)
)


# ============================================================
# TOP GENES BY P/LP BURDEN
# ============================================================

print()
print("=" * 75)
print("TOP 20 GENES BY P/LP BURDEN")
print("=" * 75)

print(
    df[
        [
            "Gene",
            "Groups",
            "VUS_unique_positions",
            "PLP_unique_positions",
            "ClinVar_unique_positions",
            "Mean_nearest_PLP"
        ]
    ]
    .sort_values(
        "PLP_unique_positions",
        ascending=False
    )
    .head(20)
    .to_string(index=False)
)


# ============================================================
# CORRELATIONS
# ============================================================

def report_spearman(x, y, label):

    mask = (
        np.isfinite(x)
        &
        np.isfinite(y)
    )

    rho, p = spearmanr(
        x[mask],
        y[mask]
    )

    print(
        f"{label}: rho={rho:.4f}, p={p:.6g}"
    )


print()
print("=" * 75)
print("GENE-LEVEL CORRELATIONS")
print("=" * 75)

report_spearman(
    df["VUS_unique_positions"].to_numpy(),
    df["Mean_nearest_PLP"].to_numpy(),
    "VUS burden vs mean proximity"
)

report_spearman(
    df["PLP_unique_positions"].to_numpy(),
    df["Mean_nearest_PLP"].to_numpy(),
    "P/LP burden vs mean proximity"
)

report_spearman(
    df["ClinVar_unique_positions"].to_numpy(),
    df["Mean_nearest_PLP"].to_numpy(),
    "Total ClinVar burden vs mean proximity"
)

report_spearman(
    df["VUS_fraction_of_ClinVar"].to_numpy(),
    df["Mean_nearest_PLP"].to_numpy(),
    "VUS fraction vs mean proximity"
)


# ============================================================
# DISTRIBUTION
# ============================================================

print()
print("=" * 75)
print("GENE-LEVEL VARIANT DISTRIBUTION")
print("=" * 75)

print(
    "Median VUS unique positions:",
    df["VUS_unique_positions"].median()
)

print(
    "Median P/LP unique positions:",
    df["PLP_unique_positions"].median()
)

print(
    "Median ClinVar unique positions:",
    df["ClinVar_unique_positions"].median()
)

print(
    "Median VUS fraction:",
    df["VUS_fraction_of_ClinVar"].median()
)

print(
    "Median P/LP fraction:",
    df["PLP_fraction_of_ClinVar"].median()
)


# ============================================================
# SAVE
# ============================================================

out = BASE + r"\gene_level_analysis.tsv"

df.to_csv(
    out,
    sep="\t",
    index=False
)

group_out = BASE + r"\gene_level_group_results.tsv"

groups.to_csv(
    group_out,
    sep="\t",
    index=False
)

print()
print("=" * 75)
print("FILES SAVED")
print("=" * 75)

print(out)
print(group_out)

print()
print("=" * 75)
print("DONE")
print("=" * 75)
