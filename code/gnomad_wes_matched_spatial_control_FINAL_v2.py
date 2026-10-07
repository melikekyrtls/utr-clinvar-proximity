import pandas as pd
import numpy as np

# ============================================================
# FINAL gNOMAD WES MATCHED SPATIAL CONTROL
# ============================================================

BASE = r"..\data"

CLINVAR = BASE + r"\canonical_utr_vus_gnomad_annotations_clean.tsv"
WES = BASE + r"\gnomad_wes_true_utr_variants.tsv"
MANE = BASE + r"\MANE_UTR_coordinates_v2.tsv"

OUT = BASE + r"\gnomad_wes_matched_spatial_control_FINAL_v2.tsv"

N_PERM = 5000
SEED = 42

rng = np.random.default_rng(SEED)

print("=" * 80)
print("FINAL CLINVAR vs gNOMAD WES MATCHED SPATIAL CONTROL")
print("=" * 80)
print("Permutations:", N_PERM)
print("Seed:", SEED)

# ============================================================
# LOAD
# ============================================================

clin = pd.read_csv(
    CLINVAR,
    sep="\t",
    low_memory=False
)

wes = pd.read_csv(
    WES,
    sep="\t",
    low_memory=False
)

mane = pd.read_csv(
    MANE,
    sep="\t",
    low_memory=False
)

print()
print("ClinVar rows:", len(clin))
print("WES rows:", len(wes))
print("MANE rows:", len(mane))

# ============================================================
# NORMALIZE
# ============================================================

for df in [clin, wes, mane]:

    df["Chromosome"] = (
        df["Chromosome"]
        .astype(str)
        .str.replace("chr", "", regex=False)
        .str.strip()
    )

clin["Position"] = pd.to_numeric(
    clin["Position"],
    errors="coerce"
)

wes["Position"] = pd.to_numeric(
    wes["Position"],
    errors="coerce"
)

mane["Start"] = pd.to_numeric(
    mane["Start"],
    errors="coerce"
)

mane["Stop"] = pd.to_numeric(
    mane["Stop"],
    errors="coerce"
)

clin = clin.dropna(
    subset=["Chromosome", "Position"]
).copy()

wes = wes.dropna(
    subset=["Chromosome", "Position"]
).copy()

mane = mane.dropna(
    subset=["Chromosome", "Start", "Stop"]
).copy()

clin["Position"] = clin["Position"].astype(int)
wes["Position"] = wes["Position"].astype(int)
mane["Start"] = mane["Start"].astype(int)
mane["Stop"] = mane["Stop"].astype(int)

# ============================================================
# BUILD TRUE TRANSCRIPT-ORIENTED MANE MAP
# ============================================================

print()
print("Building transcript-oriented MANE UTR map...")

segments = []

for (tx, utr), g in mane.groupby(
    ["Transcript", "UTR_Type"],
    sort=False
):

    strand = str(g["Strand"].iloc[0])
    chrom = str(g["Chromosome"].iloc[0])

    intervals = sorted(
        [
            (int(r.Start), int(r.Stop))
            for _, r in g.iterrows()
        ]
    )

    # Merge overlapping or adjacent UTR segments
    merged = []

    for start, stop in intervals:

        if (
            not merged
            or start > merged[-1][1] + 1
        ):
            merged.append([start, stop])

        else:
            merged[-1][1] = max(
                merged[-1][1],
                stop
            )

    # Transcript orientation
    if strand == "+":
        ordered = merged
    else:
        ordered = merged[::-1]

    total_length = sum(
        stop - start + 1
        for start, stop in merged
    )

    cumulative = 0

    for start, stop in ordered:

        length = stop - start + 1

        segments.append({
            "Chromosome": chrom,
            "Transcript": tx,
            "UTR_Type": utr,
            "Strand": strand,
            "Start": start,
            "Stop": stop,
            "CumStart": cumulative + 1,
            "CumStop": cumulative + length,
            "UTR_Length": total_length
        })

        cumulative += length

segments = pd.DataFrame(segments)

print(
    "Transcript/UTR segments:",
    len(segments)
)

# ============================================================
# CREATE FAST SEGMENT LOOKUP
# ============================================================

print()
print("Creating UTR lookup...")

segment_lookup = {}

for key, g in segments.groupby(
    ["Chromosome", "Transcript", "UTR_Type"],
    sort=False
):

    segment_lookup[key] = g[
        [
            "Start",
            "Stop",
            "CumStart",
            "CumStop",
            "UTR_Length",
            "Strand"
        ]
    ].to_numpy()

print(
    "UTR lookup groups:",
    len(segment_lookup)
)

# ============================================================
# UNIQUE CLINVAR COORDINATES
# ============================================================

clin_u = clin[
    [
        "Gene",
        "UTR_Type",
        "Transcript",
        "Position",
        "Chromosome"
    ]
].drop_duplicates(
    subset=[
        "Transcript",
        "UTR_Type",
        "Position"
    ]
).copy()

clin_u["Group"] = (
    clin_u["UTR_Type"].astype(str)
    + "|"
    + clin_u["Transcript"].astype(str)
)

print()
print(
    "Unique ClinVar spatial coordinates:",
    len(clin_u)
)

# ============================================================
# MAP ALL CLINVAR COORDINATES TO MANE
# ============================================================

print()
print("Mapping all ClinVar coordinates to MANE...")

mapped_rows = []

for row in clin_u.itertuples(index=False):

    key = (
        str(row.Chromosome),
        str(row.Transcript),
        str(row.UTR_Type)
    )

    if key not in segment_lookup:
        continue

    pos = int(row.Position)

    arr = segment_lookup[key]

    for (
        start,
        stop,
        cum_start,
        cum_stop,
        utr_length,
        strand
    ) in arr:

        if start <= pos <= stop:

            if strand == "+":

                cumulative_pos = (
                    cum_start
                    + pos
                    - start
                )

            else:

                cumulative_pos = (
                    cum_start
                    + stop
                    - pos
                )

            d5 = cumulative_pos - 1

            d3 = (
                utr_length
                - cumulative_pos
            )

            nearest = min(d5, d3)

            frac = (
                cumulative_pos
                / utr_length
            )

            edge10 = (
                frac <= 0.10
                or frac >= 0.90
            )

            mapped_rows.append({
                "Gene": row.Gene,
                "UTR_Type": row.UTR_Type,
                "Transcript": row.Transcript,
                "Chromosome": row.Chromosome,
                "Position": pos,
                "Group": row.Group,
                "UTR_Length": utr_length,
                "Distance_to_5prime_edge": d5,
                "Distance_to_3prime_edge": d3,
                "Distance_to_Nearest_Edge": nearest,
                "Relative_Position": frac,
                "Edge10": edge10
            })

            break

clin_spatial = pd.DataFrame(mapped_rows)

print(
    "ClinVar coordinates mapped to MANE:",
    len(clin_spatial)
)

print(
    "ClinVar coordinates NOT mapped:",
    len(clin_u) - len(clin_spatial)
)

if len(clin_spatial) == 0:
    raise RuntimeError(
        "No ClinVar coordinates mapped to MANE."
    )

# ============================================================
# MAP gNOMAD WES TO SAME UTR SPACE
# ============================================================

print()
print("Mapping gNOMAD WES variants to same UTR space...")

wes_groups = {
    key: g
    for key, g in wes.groupby(
        ["Transcript", "UTR_Type"],
        sort=False
    )
}

wes_parts = []

for _, seg in segments.iterrows():

    key = (
        seg["Transcript"],
        seg["UTR_Type"]
    )

    if key not in wes_groups:
        continue

    w = wes_groups[key]

    mask = (
        (w["Chromosome"] == seg["Chromosome"])
        &
        (w["Position"] >= seg["Start"])
        &
        (w["Position"] <= seg["Stop"])
    )

    if not mask.any():
        continue

    x = w.loc[mask].copy()

    if seg["Strand"] == "+":

        x["UTR_Cumulative_Position"] = (
            seg["CumStart"]
            + x["Position"]
            - seg["Start"]
        )

    else:

        x["UTR_Cumulative_Position"] = (
            seg["CumStart"]
            + seg["Stop"]
            - x["Position"]
        )

    x["UTR_Length_corrected"] = (
        seg["UTR_Length"]
    )

    wes_parts.append(x)

if not wes_parts:
    raise RuntimeError(
        "No WES variants mapped to MANE UTRs."
    )

wes2 = pd.concat(
    wes_parts,
    ignore_index=True
)

# Remove exact spatial duplicates
wes2 = wes2.drop_duplicates(
    subset=[
        "Chromosome",
        "Position",
        "Ref",
        "Alt",
        "Transcript",
        "UTR_Type"
    ]
)

wes2["WES_AF"] = pd.to_numeric(
    wes2["WES_AF"],
    errors="coerce"
)

# Population-observed variants only
wes2 = wes2[
    wes2["WES_AF"] > 0
].copy()

wes2["Distance_to_5prime_edge"] = (
    wes2["UTR_Cumulative_Position"] - 1
)

wes2["Distance_to_3prime_edge"] = (
    wes2["UTR_Length_corrected"]
    - wes2["UTR_Cumulative_Position"]
)

wes2["Distance_to_Nearest_Edge"] = np.minimum(
    wes2["Distance_to_5prime_edge"],
    wes2["Distance_to_3prime_edge"]
)

wes2["Relative_Position_Fraction"] = (
    wes2["UTR_Cumulative_Position"]
    / wes2["UTR_Length_corrected"]
)

wes2["Edge10"] = (
    (wes2["Relative_Position_Fraction"] <= 0.10)
    |
    (wes2["Relative_Position_Fraction"] >= 0.90)
)

wes2["Group"] = (
    wes2["UTR_Type"].astype(str)
    + "|"
    + wes2["Transcript"].astype(str)
)

print(
    "gNOMAD WES AF>0 UTR rows:",
    len(wes2)
)

# ============================================================
# CREATE WES SAMPLING CACHE
# ============================================================

print()
print("Creating WES sampling cache...")

wes_cache = {}

for group, g in wes2.groupby(
    "Group",
    sort=False
):

    # One spatial observation per genomic coordinate
    g = g.drop_duplicates(
        subset=[
            "Chromosome",
            "Position"
        ]
    )

    wes_cache[group] = {
        "edge": g[
            "Edge10"
        ].to_numpy(
            dtype=np.int8
        ),
        "distance": g[
            "Distance_to_Nearest_Edge"
        ].to_numpy(
            dtype=float
        )
    }

print(
    "WES groups available:",
    len(wes_cache)
)

# ============================================================
# MATCH CLINVAR GROUPS TO WES GROUPS
# ============================================================

print()
print("Matching ClinVar groups to WES sampling space...")

group_plan = []

for group, g in clin_spatial.groupby(
    "Group",
    sort=False
):

    if group not in wes_cache:
        continue

    n_clin = len(g)
    n_wes = len(
        wes_cache[group]["edge"]
    )

    if n_wes == 0:
        continue

    group_plan.append(
        (
            group,
            n_clin,
            n_wes
        )
    )

matched_n = sum(
    x[1]
    for x in group_plan
)

print(
    "Matched groups:",
    len(group_plan)
)

print(
    "ClinVar coordinates represented:",
    matched_n
)

print(
    "ClinVar coordinates not represented:",
    len(clin_spatial) - matched_n
)

if len(group_plan) == 0:
    raise RuntimeError(
        "No ClinVar groups could be matched to WES."
    )

# ============================================================
# REPLACEMENT GROUP COUNT
# ============================================================

replacement_groups = sum(
    n_wes < n_clin
    for _, n_clin, n_wes in group_plan
)

print(
    "Groups requiring replacement sampling:",
    replacement_groups
)

# ============================================================
# OBSERVED CLINVAR DISTRIBUTION
# ============================================================

matched_group_names = {
    x[0]
    for x in group_plan
}

clin_obs = clin_spatial[
    clin_spatial["Group"].isin(
        matched_group_names
    )
].copy()

obs_edge = (
    clin_obs["Edge10"]
    .astype(int)
    .mean()
)

obs_dist = (
    clin_obs[
        "Distance_to_Nearest_Edge"
    ]
    .astype(float)
    .mean()
)

print()
print(
    "Observed matched ClinVar Edge10:",
    obs_edge
)

print(
    "Observed matched ClinVar mean edge distance:",
    obs_dist
)

# ============================================================
# PERMUTATION
# ============================================================

print()
print(
    "Running",
    N_PERM,
    "WES matched permutations..."
)

perm_edge = np.empty(
    N_PERM,
    dtype=float
)

perm_dist = np.empty(
    N_PERM,
    dtype=float
)

replacement_operations = 0

for p in range(N_PERM):

    total_edge = 0
    total_dist = 0.0
    total_n = 0

    for group, n_clin, n_wes in group_plan:

        cache = wes_cache[group]

        replace = (
            n_wes < n_clin
        )

        if replace:
            replacement_operations += 1

        idx = rng.choice(
            n_wes,
            size=n_clin,
            replace=replace
        )

        total_edge += int(
            cache["edge"][idx].sum()
        )

        total_dist += float(
            cache["distance"][idx].sum()
        )

        total_n += n_clin

    perm_edge[p] = (
        total_edge / total_n
    )

    perm_dist[p] = (
        total_dist / total_n
    )

    if (p + 1) % 500 == 0:
        print(
            f"{p + 1}/{N_PERM}"
        )

# ============================================================
# EMPIRICAL P VALUES
# ============================================================

# Hypothesis:
# ClinVar VUS are MORE edge-proximal than WES population variation.

# Edge10:
# Larger = more edge-proximal

p_edge = (
    np.sum(
        perm_edge >= obs_edge
    ) + 1
) / (
    N_PERM + 1
)

# Nearest-edge distance:
# Smaller = more edge-proximal

p_dist = (
    np.sum(
        perm_dist <= obs_dist
    ) + 1
) / (
    N_PERM + 1
)

# ============================================================
# NULL INTERVALS
# ============================================================

edge_low, edge_high = np.percentile(
    perm_edge,
    [2.5, 97.5]
)

dist_low, dist_high = np.percentile(
    perm_dist,
    [2.5, 97.5]
)

# ============================================================
# FINAL RESULTS
# ============================================================

print()
print("=" * 80)
print("FINAL WES MATCHED CONTROL RESULTS")
print("=" * 80)

print(
    "Unique ClinVar coordinates:",
    len(clin_u)
)

print(
    "ClinVar coordinates mapped to MANE:",
    len(clin_spatial)
)

print(
    "Matched ClinVar coordinates:",
    matched_n
)

print(
    "Matched groups:",
    len(group_plan)
)

print(
    "Groups requiring replacement:",
    replacement_groups
)

print()

print(
    "Observed ClinVar Edge10:",
    obs_edge
)

print(
    "WES matched null mean Edge10:",
    perm_edge.mean()
)

print(
    "WES null 95% interval Edge10:",
    edge_low,
    edge_high
)

print(
    "Edge10 enrichment p:",
    p_edge
)

print()

print(
    "Observed ClinVar mean edge distance:",
    obs_dist
)

print(
    "WES matched null mean edge distance:",
    perm_dist.mean()
)

print(
    "WES null 95% interval distance:",
    dist_low,
    dist_high
)

print(
    "Distance-to-edge enrichment p:",
    p_dist
)

print()

print(
    "Valid permutations:",
    N_PERM
)

print(
    "Random seed:",
    SEED
)

print(
    "Total replacement sampling operations:",
    replacement_operations
)

# ============================================================
# SAVE SUMMARY
# ============================================================

summary = pd.DataFrame([
    {
        "Test": "WES matched Edge10",
        "ClinVar_unique_coordinates": len(clin_u),
        "ClinVar_MANE_mapped": len(clin_spatial),
        "Matched_ClinVar_N": matched_n,
        "Matched_groups": len(group_plan),
        "Replacement_groups": replacement_groups,
        "Observed": obs_edge,
        "Null_mean": perm_edge.mean(),
        "Null_2.5pct": edge_low,
        "Null_97.5pct": edge_high,
        "P_value": p_edge,
        "N_permutations": N_PERM,
        "Seed": SEED
    },
    {
        "Test": "WES matched mean edge distance",
        "ClinVar_unique_coordinates": len(clin_u),
        "ClinVar_MANE_mapped": len(clin_spatial),
        "Matched_ClinVar_N": matched_n,
        "Matched_groups": len(group_plan),
        "Replacement_groups": replacement_groups,
        "Observed": obs_dist,
        "Null_mean": perm_dist.mean(),
        "Null_2.5pct": dist_low,
        "Null_97.5pct": dist_high,
        "P_value": p_dist,
        "N_permutations": N_PERM,
        "Seed": SEED
    }
])

summary.to_csv(
    OUT,
    sep="\t",
    index=False
)

print()
print(
    "Saved:",
    OUT
)

print("=" * 80)
print("DONE")
print("=" * 80)
