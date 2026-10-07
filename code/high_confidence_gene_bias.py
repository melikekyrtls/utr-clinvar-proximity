import pandas as pd
import numpy as np

BASE = r"..\data"

CANONICAL = BASE + r"\canonical_utr_vus_dataset.tsv"
CLINVAR = BASE + r"\ClinVar_UTR_GRCh38.tsv"
MANE = BASE + r"\MANE_UTR_coordinates_v2.tsv"

N_PERM = 5000
SEED = 42

# ------------------------------------------------------------
# LOAD
# ------------------------------------------------------------

vus = pd.read_csv(
    CANONICAL,
    sep="\t"
)

cv = pd.read_csv(
    CLINVAR,
    sep="\t"
)

mane = pd.read_csv(
    MANE,
    sep="\t"
)

group_cols = [
    "Gene",
    "UTR_Type",
    "Transcript"
]

# ------------------------------------------------------------
# CLINVAR POSITION
# ------------------------------------------------------------

cv["PositionVCF"] = pd.to_numeric(
    cv["PositionVCF"],
    errors="coerce"
)

cv["Start"] = pd.to_numeric(
    cv["Start"],
    errors="coerce"
)

cv["Position"] = cv["PositionVCF"]

fallback = (
    cv["Position"].isna()
    |
    (cv["Position"] == -1)
)

cv.loc[fallback, "Position"] = cv.loc[
    fallback,
    "Start"
]

# ------------------------------------------------------------
# P/LP
# ------------------------------------------------------------

plp_values = {
    "Pathogenic",
    "Likely pathogenic",
    "Pathogenic/Likely pathogenic",
    "Likely pathogenic/Pathogenic"
}

cv["Is_PLP"] = (
    cv["ClinicalSignificance"]
    .isin(plp_values)
)

# ------------------------------------------------------------
# REVIEW STATUS
# ------------------------------------------------------------

cv["ReviewStatus_clean"] = (
    cv["ReviewStatus"]
    .astype(str)
    .str.strip()
    .str.lower()
)

high_conf_status = {
    "criteria provided, multiple submitters, no conflicts",
    "reviewed by expert panel",
    "practice guideline"
}

cv["HighConfidence"] = (
    cv["ReviewStatus_clean"]
    .isin(high_conf_status)
)

# ------------------------------------------------------------
# UTR MAPPING FILTER
# ------------------------------------------------------------

cv = cv[
    cv["UTR_Transcript"].notna()
    &
    (
        cv["UTR_Transcript"]
        .astype(str)
        .str.lower()
        != "na"
    )
    &
    cv["UTR_Type"].notna()
    &
    cv["GeneSymbol"].notna()
    &
    cv["Position"].notna()
].copy()

cv["Gene"] = cv["GeneSymbol"]
cv["Transcript"] = cv["UTR_Transcript"]

# ------------------------------------------------------------
# BASIC COUNTS
# ------------------------------------------------------------

print("=" * 75)
print("HIGH-CONFIDENCE P/LP + GENE-LEVEL BIAS ANALYSIS")
print("=" * 75)

print("Canonical VUS rows:", len(vus))

# ------------------------------------------------------------
# GENE-LEVEL VUS DISTRIBUTION
#
# One genomic position counted once per gene.
# ------------------------------------------------------------

vus_gene = (
    vus[
        ["Gene", "Position"]
    ]
    .drop_duplicates()
)

vus_gene_counts = (
    vus_gene
    .groupby("Gene")
    .size()
    .reset_index(
        name="VUS_unique_positions"
    )
)

# ------------------------------------------------------------
# GENE-LEVEL HIGH-CONFIDENCE P/LP DISTRIBUTION
# ------------------------------------------------------------

hc_plp = cv[
    cv["Is_PLP"]
    &
    cv["HighConfidence"]
].copy()

hc_plp_gene = (
    hc_plp[
        ["Gene", "Position"]
    ]
    .drop_duplicates()
)

hc_gene_counts = (
    hc_plp_gene
    .groupby("Gene")
    .size()
    .reset_index(
        name="HC_PLP_unique_positions"
    )
)

# ------------------------------------------------------------
# MERGE
# ------------------------------------------------------------

gene_table = vus_gene_counts.merge(
    hc_gene_counts,
    on="Gene",
    how="outer"
).fillna(0)

gene_table["VUS_unique_positions"] = (
    gene_table["VUS_unique_positions"]
    .astype(int)
)

gene_table["HC_PLP_unique_positions"] = (
    gene_table["HC_PLP_unique_positions"]
    .astype(int)
)

gene_table = gene_table.sort_values(
    "HC_PLP_unique_positions",
    ascending=False
)

# ------------------------------------------------------------
# GINI FUNCTION
# ------------------------------------------------------------

def gini(values):

    x = np.asarray(
        values,
        dtype=float
    )

    x = np.sort(x)

    if len(x) == 0:
        return np.nan

    total = x.sum()

    if total == 0:
        return 0.0

    n = len(x)

    return (
        (2 * np.sum(
            (np.arange(1, n + 1)) * x
        ))
        / (n * total)
        - (n + 1) / n
    )

# ------------------------------------------------------------
# TOP-10% CONCENTRATION
# ------------------------------------------------------------

n_genes = len(gene_table)

top_n = max(
    1,
    int(np.ceil(0.10 * n_genes))
)

hc_sorted = gene_table.sort_values(
    "HC_PLP_unique_positions",
    ascending=False
)

vus_sorted = gene_table.sort_values(
    "VUS_unique_positions",
    ascending=False
)

top_hc_share = (
    hc_sorted
    .head(top_n)["HC_PLP_unique_positions"]
    .sum()
    /
    max(1, hc_sorted["HC_PLP_unique_positions"].sum())
)

top_vus_share = (
    vus_sorted
    .head(top_n)["VUS_unique_positions"]
    .sum()
    /
    max(1, vus_sorted["VUS_unique_positions"].sum())
)

print("\n" + "=" * 75)
print("GENE-LEVEL CONCENTRATION")
print("=" * 75)

print("Genes represented by VUS:", len(vus_gene_counts))
print(
    "Genes represented by high-confidence P/LP:",
    len(hc_gene_counts)
)

print(
    "Total unique VUS gene-position combinations:",
    len(vus_gene)
)

print(
    "Total unique HC P/LP gene-position combinations:",
    len(hc_plp_gene)
)

print(
    "VUS Gini:",
    f"{gini(gene_table['VUS_unique_positions']):.4f}"
)

print(
    "HC P/LP Gini:",
    f"{gini(gene_table['HC_PLP_unique_positions']):.4f}"
)

print(
    "Top 10% genes share of VUS positions:",
    f"{top_vus_share:.4%}"
)

print(
    "Top 10% genes share of HC P/LP positions:",
    f"{top_hc_share:.4%}"
)

print("\nTop genes by HC P/LP count:")

print(
    gene_table[
        [
            "Gene",
            "VUS_unique_positions",
            "HC_PLP_unique_positions"
        ]
    ]
    .head(20)
    .to_string(index=False)
)

# ------------------------------------------------------------
# GROUP-LEVEL HIGH-CONFIDENCE P/LP MAP
# ------------------------------------------------------------

hc_plp_group = (
    hc_plp[
        group_cols + ["Position"]
    ]
    .drop_duplicates()
)

hc_plp_positions = {}

for key, group in hc_plp_group.groupby(
    group_cols
):

    hc_plp_positions[key] = np.array(
        sorted(
            set(
                group["Position"]
                .astype(int)
            )
        ),
        dtype=np.int64
    )

# ------------------------------------------------------------
# ALL P/LP POSITIONS
#
# Used only to prevent randomization into any clinically
# classified P/LP coordinate.
# ------------------------------------------------------------

all_plp_group = (
    cv[
        cv["Is_PLP"]
    ][
        group_cols + ["Position"]
    ]
    .drop_duplicates()
)

all_plp_positions = {}

for key, group in all_plp_group.groupby(
    group_cols
):

    all_plp_positions[key] = set(
        group["Position"]
        .astype(int)
        .tolist()
    )

# ------------------------------------------------------------
# MANE UTR POSITION POOLS
# ------------------------------------------------------------

mane["Start"] = pd.to_numeric(
    mane["Start"],
    errors="coerce"
)

mane["Stop"] = pd.to_numeric(
    mane["Stop"],
    errors="coerce"
)

mane = mane.dropna(
    subset=[
        "Start",
        "Stop",
        "Transcript",
        "UTR_Type"
    ]
).copy()

mane["Start"] = mane["Start"].astype(int)
mane["Stop"] = mane["Stop"].astype(int)

# Build position arrays for each transcript/UTR.
mane_pools = {}

for key, group in mane.groupby(
    ["Transcript", "UTR_Type"]
):

    intervals = (
        group[
            ["Start", "Stop"]
        ]
        .sort_values("Start")
        .to_numpy()
    )

    merged = []

    for start, stop in intervals:

        start = int(start)
        stop = int(stop)

        if not merged:
            merged.append(
                [start, stop]
            )
        else:
            if start <= merged[-1][1] + 1:
                merged[-1][1] = max(
                    merged[-1][1],
                    stop
                )
            else:
                merged.append(
                    [start, stop]
                )

    pieces = []

    for start, stop in merged:
        pieces.append(
            np.arange(
                start,
                stop + 1,
                dtype=np.int64
            )
        )

    if pieces:
        mane_pools[key] = np.concatenate(
            pieces
        )

# ------------------------------------------------------------
# BUILD OBSERVED DATASET
# ------------------------------------------------------------

vus["Position"] = pd.to_numeric(
    vus["Position"],
    errors="coerce"
)

vus = vus.dropna(
    subset=["Position"]
).copy()

vus["Position"] = vus["Position"].astype(int)

observed_distances = []
group_data = []

for key, group in vus.groupby(
    group_cols
):

    # Only groups containing at least one
    # high-confidence P/LP anchor.
    if key not in hc_plp_positions:
        continue

    hc_positions = hc_plp_positions[key]

    positions = group[
        "Position"
    ].to_numpy(
        dtype=np.int64
    )

    unique_positions, multiplicities = np.unique(
        positions,
        return_counts=True
    )

    group_distances = []

    for pos in positions:

        diffs = np.abs(
            hc_positions
            -
            np.int64(pos)
        )

        diffs = diffs[
            diffs > 0
        ]

        if len(diffs) == 0:
            continue

        group_distances.append(
            int(diffs.min())
        )

    if len(group_distances) == 0:
        continue

    observed_distances.extend(
        group_distances
    )

    group_data.append(
        (
            key,
            unique_positions,
            multiplicities
        )
    )

print("\n" + "=" * 75)
print("HIGH-CONFIDENCE P/LP PROXIMITY DATASET")
print("=" * 75)

print(
    "High-confidence P/LP rows:",
    len(hc_plp)
)

print(
    "High-confidence P/LP unique group-position coordinates:",
    len(hc_plp_group)
)

print(
    "Groups with HC P/LP anchors:",
    len(group_data)
)

print(
    "Observed VUS distances:",
    len(observed_distances)
)

if len(observed_distances) == 0:
    raise RuntimeError(
        "No VUS observations could be matched to HC P/LP anchors."
    )

observed_mean = np.mean(
    observed_distances
)

print(
    "Observed mean distance:",
    f"{observed_mean:.4f}"
)

# ------------------------------------------------------------
# WHOLE-UTR PERMUTATION
# ------------------------------------------------------------

rng = np.random.default_rng(
    SEED
)

whole_null = []

skipped_groups = 0
replacement_groups = 0

for iteration in range(N_PERM):

    perm_distances = []

    for key, unique_positions, multiplicities in group_data:

        transcript = key[2]
        utr_type = key[1]

        mane_key = (
            transcript,
            utr_type
        )

        if mane_key not in mane_pools:
            skipped_groups += 1
            continue

        pool = mane_pools[mane_key]

        plp_set = all_plp_positions.get(
            key,
            set()
        )

        if len(plp_set) > 0:

            pool = pool[
                ~np.isin(
                    pool,
                    np.array(
                        list(plp_set),
                        dtype=np.int64
                    )
                )
            ]

        n_unique = len(
            unique_positions
        )

        if len(pool) == 0:
            continue

        replace = (
            n_unique > len(pool)
        )

        if replace:
            replacement_groups += 1

        sampled = rng.choice(
            pool,
            size=n_unique,
            replace=replace
        )

        shuffled_mult = rng.permutation(
            multiplicities
        )

        hc_positions = hc_plp_positions[key]

        for pos, mult in zip(
            sampled,
            shuffled_mult
        ):

            diffs = np.abs(
                hc_positions
                -
                np.int64(pos)
            )

            diffs = diffs[
                diffs > 0
            ]

            if len(diffs) == 0:
                continue

            nearest = int(
                diffs.min()
            )

            perm_distances.extend(
                [nearest] * int(mult)
            )

    if len(perm_distances) > 0:

        whole_null.append(
            np.mean(perm_distances)
        )

whole_null = np.array(
    whole_null,
    dtype=float
)

# ------------------------------------------------------------
# EMPIRICAL CLINVAR SAMPLE-SPACE PERMUTATION
# ------------------------------------------------------------

# Build empirical non-P/LP ClinVar pools
all_cv_positions = (
    cv[
        group_cols + ["Position"]
    ]
    .drop_duplicates()
)

all_plp_df = (
    cv[
        cv["Is_PLP"]
    ][
        group_cols + ["Position"]
    ]
    .drop_duplicates()
)

empirical_pool_df = all_cv_positions.merge(
    all_plp_df.assign(
        IS_PLP_POSITION=True
    ),
    on=group_cols + ["Position"],
    how="left"
)

empirical_pool_df["IS_PLP_POSITION"] = (
    empirical_pool_df["IS_PLP_POSITION"]
    .fillna(False)
    .astype(bool)
)

empirical_pool_df = empirical_pool_df[
    empirical_pool_df["IS_PLP_POSITION"] == False
]

empirical_pools = {}

for key, group in empirical_pool_df.groupby(
    group_cols
):

    empirical_pools[key] = np.array(
        sorted(
            set(
                group["Position"]
                .astype(int)
            )
        ),
        dtype=np.int64
    )

empirical_null = []

for iteration in range(N_PERM):

    perm_distances = []

    for key, unique_positions, multiplicities in group_data:

        if key not in empirical_pools:
            continue

        pool = empirical_pools[key]

        if len(pool) == 0:
            continue

        n_unique = len(
            unique_positions
        )

        sampled = rng.choice(
            pool,
            size=n_unique,
            replace=True
        )

        shuffled_mult = rng.permutation(
            multiplicities
        )

        hc_positions = hc_plp_positions[key]

        for pos, mult in zip(
            sampled,
            shuffled_mult
        ):

            diffs = np.abs(
                hc_positions
                -
                np.int64(pos)
            )

            diffs = diffs[
                diffs > 0
            ]

            if len(diffs) == 0:
                continue

            nearest = int(
                diffs.min()
            )

            perm_distances.extend(
                [nearest] * int(mult)
            )

    if len(perm_distances) > 0:

        empirical_null.append(
            np.mean(perm_distances)
        )

empirical_null = np.array(
    empirical_null,
    dtype=float
)

# ------------------------------------------------------------
# REPORT FUNCTION
# ------------------------------------------------------------

def report_null(name, values):

    null_mean = values.mean()

    low = np.percentile(
        values,
        2.5
    )

    high = np.percentile(
        values,
        97.5
    )

    extreme = np.sum(
        values <= observed_mean
    )

    p = (
        extreme + 1
    ) / (
        len(values) + 1
    )

    ratio = (
        null_mean
        /
        observed_mean
    )

    print(
        "\n" +
        "=" * 75
    )

    print(name)

    print("=" * 75)

    print(
        "Observed mean:",
        f"{observed_mean:.4f}"
    )

    print(
        "Null mean:",
        f"{null_mean:.4f}"
    )

    print(
        "Null 95% interval:",
        f"{low:.4f} to {high:.4f}"
    )

    print(
        "Random / observed ratio:",
        f"{ratio:.4f}"
    )

    print(
        "Valid permutations:",
        len(values)
    )

    print(
        "Extreme permutations:",
        int(extreme)
    )

    print(
        "One-sided empirical p:",
        p
    )


report_null(
    "WHOLE-UTR HIGH-CONFIDENCE P/LP NULL",
    whole_null
)

report_null(
    "EMPIRICAL CLINVAR SAMPLE-SPACE HIGH-CONFIDENCE P/LP NULL",
    empirical_null
)

# ------------------------------------------------------------
# SAVE GENE TABLE
# ------------------------------------------------------------

gene_out = BASE + r"\gene_level_hc_plp_bias.tsv"

gene_table.to_csv(
    gene_out,
    sep="\t",
    index=False
)

print(
    "\nGene-level table saved:",
    gene_out
)

print(
    "\n" +
    "=" * 75
)

print("DONE")

print("=" * 75)
