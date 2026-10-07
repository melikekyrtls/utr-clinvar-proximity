import pandas as pd
import numpy as np

# ============================================================
# FINAL GNOMAD WGS MATCHED SPATIAL CONTROL
# 5,000 permutations, seed=42
#
# Null:
# Within each Transcript + UTR_Type group, randomly sample
# the same number of gnomAD WGS AF>0 positions as observed
# ClinVar unique spatial coordinates.
#
# Tests:
# 1. Edge10 enrichment
# 2. Mean distance to nearest UTR edge
# ============================================================

SEED = 42
N_PERM = 5000

CLINVAR = r"..\data\canonical_utr_vus_with_alleles.tsv"
GNOMAD = r"..\data\gnomad_wgs_true_utr_variants.tsv"
MANE = r"..\data\MANE_UTR_coordinates_v2.tsv"

OUTPUT = r"..\data\gnomad_wgs_matched_permutation_FINAL.tsv"

rng = np.random.default_rng(SEED)

# ============================================================
# LOAD
# ============================================================

print("=" * 75)
print("FINAL GNOMAD WGS MATCHED SPATIAL CONTROL")
print("=" * 75)

v = pd.read_csv(CLINVAR, sep="\t", dtype=str)
g = pd.read_csv(GNOMAD, sep="\t", low_memory=False)
mane = pd.read_csv(MANE, sep="\t")

v["Position"] = pd.to_numeric(v["Position"], errors="coerce")

v["Chromosome"] = (
    v["Chromosome"]
    .astype(str)
    .str.replace("^chr", "", regex=True)
)

mane["Chromosome"] = (
    mane["Chromosome"]
    .astype(str)
    .str.replace("^chr", "", regex=True)
)

mane["Start"] = pd.to_numeric(mane["Start"], errors="coerce")
mane["Stop"] = pd.to_numeric(mane["Stop"], errors="coerce")

# ============================================================
# UNIQUE CLINVAR SPATIAL COORDINATES
# ============================================================

v = v.dropna(subset=["Position"])

v = v.drop_duplicates(
    subset=["Transcript", "UTR_Type", "Position"]
).copy()

print("\nUnique ClinVar spatial coordinates:", len(v))

# ============================================================
# MAP CLINVAR TO TRANSCRIPT-ORIENTED UTR SPACE
# ============================================================

rows = []

for (transcript, utr_type), sub_v in v.groupby(
    ["Transcript", "UTR_Type"],
    sort=False
):

    segs = mane[
        (mane["Transcript"] == transcript) &
        (mane["UTR_Type"] == utr_type)
    ].copy()

    if segs.empty:
        continue

    strand = segs["Strand"].iloc[0]

    if strand == "+":
        segs = segs.sort_values("Start")
    else:
        segs = segs.sort_values("Start", ascending=False)

    cumulative = 0
    segments = []

    for _, s in segs.iterrows():

        start = int(s["Start"])
        stop = int(s["Stop"])
        length = stop - start + 1

        segments.append(
            (start, stop, cumulative, length)
        )

        cumulative += length

    total_length = cumulative

    for _, r in sub_v.iterrows():

        pos = int(r["Position"])

        for start, stop, cum, length in segments:

            if start <= pos <= stop:

                if strand == "+":
                    offset = pos - start
                else:
                    offset = stop - pos

                relative = cum + offset

                rows.append({
                    "Transcript": transcript,
                    "UTR_Type": utr_type,
                    "Position": pos,
                    "Relative_Position": relative,
                    "UTR_Length": total_length,
                    "Edge_Distance": min(
                        relative,
                        total_length - 1 - relative
                    )
                })

                break

c = pd.DataFrame(rows)

c["Edge10"] = (
    (c["Relative_Position"] + 1) / c["UTR_Length"] <= 0.10
) | (
    (c["Relative_Position"] + 1) / c["UTR_Length"] >= 0.90
)

print("ClinVar coordinates mapped to MANE:", len(c))

# ============================================================
# GNOMAD WGS: AF > 0 ONLY
# ============================================================

g["AF"] = pd.to_numeric(g["AF"], errors="coerce")

g = g[g["AF"] > 0].copy()

print("gNOMAD WGS AF>0 rows:", len(g))

# Keep one spatial record per Transcript + UTR + Position
g = g.drop_duplicates(
    subset=["Transcript", "UTR_Type", "Position"]
).copy()

print("Unique gNOMAD spatial coordinates:", len(g))

# ============================================================
# PREPARE GROUP POOLS
# ============================================================

group_cols = ["Transcript", "UTR_Type"]

clinvar_groups = {
    key: group.copy()
    for key, group in c.groupby(group_cols, sort=False)
}

gnomad_groups = {
    key: group.copy()
    for key, group in g.groupby(group_cols, sort=False)
}

matched_groups = sorted(
    set(clinvar_groups) & set(gnomad_groups)
)

print("Matched groups:", len(matched_groups))

# Only retain groups where gnomAD has at least one usable position.
matched_c = pd.concat(
    [clinvar_groups[k] for k in matched_groups],
    ignore_index=True
)

print("Matched ClinVar coordinates:", len(matched_c))

# ============================================================
# OBSERVED STATISTICS
# ============================================================

observed_edge10 = matched_c["Edge10"].mean()
observed_mean_distance = matched_c["Edge_Distance"].mean()

print("\nOBSERVED")
print("Edge10:", observed_edge10)
print("Mean edge distance:", observed_mean_distance)

# ============================================================
# PERMUTATION NULL
# ============================================================

perm_edge10 = np.empty(N_PERM)
perm_distance = np.empty(N_PERM)

replacement_groups = 0

for i in range(N_PERM):

    edge_values = []
    distance_values = []

    for key in matched_groups:

        cv = clinvar_groups[key]
        pool = gnomad_groups[key]

        n = len(cv)

        replace = len(pool) < n

        if replace:
            replacement_groups += 1

        idx = rng.choice(
            len(pool),
            size=n,
            replace=replace
        )

        sampled = pool.iloc[idx]

        edge_values.extend(
            sampled["Edge10"].to_numpy()
        )

        distance_values.extend(
            sampled["Distance_to_Nearest_Edge"].to_numpy()
        )

    perm_edge10[i] = np.mean(edge_values)
    perm_distance[i] = np.mean(distance_values)

# ============================================================
# EMPIRICAL P VALUES
#
# Edge10:
# observed HIGHER than null = enrichment
#
# Distance:
# observed LOWER than null = enrichment toward edges
# ============================================================

p_edge10 = (
    np.sum(perm_edge10 >= observed_edge10) + 1
) / (N_PERM + 1)

p_distance = (
    np.sum(perm_distance <= observed_mean_distance) + 1
) / (N_PERM + 1)

# ============================================================
# 95% NULL INTERVALS
# ============================================================

edge_low, edge_high = np.percentile(
    perm_edge10,
    [2.5, 97.5]
)

dist_low, dist_high = np.percentile(
    perm_distance,
    [2.5, 97.5]
)

# ============================================================
# OUTPUT
# ============================================================

print("\n" + "=" * 75)
print("FINAL RESULTS")
print("=" * 75)

print("\nClinVar unique coordinates:", len(v))
print("ClinVar mapped:", len(c))
print("Matched ClinVar:", len(matched_c))
print("Matched groups:", len(matched_groups))

print("\nObserved ClinVar Edge10:")
print(observed_edge10)

print("\nWGS matched null mean Edge10:")
print(perm_edge10.mean())

print("\nWGS null 95% interval Edge10:")
print(edge_low, edge_high)

print("\nEdge10 enrichment empirical p:")
print(p_edge10)

print("\nObserved ClinVar mean edge distance:")
print(observed_mean_distance)

print("\nWGS matched null mean edge distance:")
print(perm_distance.mean())

print("\nWGS null 95% interval distance:")
print(dist_low, dist_high)

print("\nDistance-to-edge enrichment empirical p:")
print(p_distance)

print("\nValid permutations:", N_PERM)
print("Random seed:", SEED)

print("\nTotal group sampling operations using replacement:")
print(replacement_groups)

# ============================================================
# SAVE PERMUTATION SUMMARY
# ============================================================

summary = pd.DataFrame({
    "Statistic": [
        "ClinVar_unique_coordinates",
        "ClinVar_mapped_coordinates",
        "Matched_ClinVar_coordinates",
        "Matched_groups",
        "Observed_Edge10",
        "Null_mean_Edge10",
        "Null_95CI_low_Edge10",
        "Null_95CI_high_Edge10",
        "Empirical_p_Edge10_enrichment",
        "Observed_mean_edge_distance",
        "Null_mean_edge_distance",
        "Null_95CI_low_distance",
        "Null_95CI_high_distance",
        "Empirical_p_distance_enrichment",
        "N_permutations",
        "Seed",
        "Replacement_sampling_operations"
    ],
    "Value": [
        len(v),
        len(c),
        len(matched_c),
        len(matched_groups),
        observed_edge10,
        perm_edge10.mean(),
        edge_low,
        edge_high,
        p_edge10,
        observed_mean_distance,
        perm_distance.mean(),
        dist_low,
        dist_high,
        p_distance,
        N_PERM,
        SEED,
        replacement_groups
    ]
})

summary.to_csv(
    OUTPUT,
    sep="\t",
    index=False
)

print("\nSaved:")
print(OUTPUT)
