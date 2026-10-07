import csv
import random
import statistics

random.seed(42)

# ============================================================
# FILES
# ============================================================

CLINVAR = "ClinVar_UTR_GRCh38.tsv"
MANE = "MANE_UTR_coordinates_v2.tsv"
CANONICAL = "canonical_utr_vus_dataset.tsv"

N_PERMUTATIONS = 5000


# ============================================================
# LOAD MANE UTR INTERVALS
# ============================================================

segments = {}

with open(MANE, encoding="utf-8") as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        chrom = row["Chromosome"]
        start = int(row["Start"])
        stop = int(row["Stop"])
        transcript = row["Transcript"]
        utr_type = row["UTR_Type"]

        key = (chrom, transcript, utr_type)

        segments.setdefault(key, []).append((start, stop))


# ============================================================
# MERGE UTR SEGMENTS
# ============================================================

merged_segments = {}

for key, intervals in segments.items():

    intervals = sorted(intervals)

    merged = []

    for start, stop in intervals:

        if not merged:
            merged.append([start, stop])

        elif start <= merged[-1][1] + 1:
            merged[-1][1] = max(
                merged[-1][1],
                stop
            )

        else:
            merged.append([start, stop])

    merged_segments[key] = merged


# ============================================================
# GENERATE ALL UTR POSITIONS
# ============================================================

utr_positions = {}

for key, intervals in merged_segments.items():

    positions = []

    for start, stop in intervals:
        positions.extend(range(start, stop + 1))

    utr_positions[key] = positions


# ============================================================
# LOAD CANONICAL VUS DATA
# ============================================================

canonical = []

with open(CANONICAL, encoding="utf-8") as f:

    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:

        canonical.append({
            "gene": row["Gene"],
            "utr": row["UTR_Type"],
            "transcript": row["Transcript"],
            "position": int(row["Position"]),
            "distance": float(row["Distance"])
        })


# ============================================================
# LOAD P/LP POSITIONS FROM CLINVAR
# ============================================================

plp = {}

with open(CLINVAR, encoding="utf-8") as f:

    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:

        if row["Assembly"] != "GRCh38":
            continue

        significance = row["ClinicalSignificance"]

        if significance not in [
            "Pathogenic",
            "Likely pathogenic",
            "Pathogenic/Likely pathogenic",
            "Likely pathogenic/Pathogenic"
        ]:
            continue

        transcript = row["UTR_Transcript"]
        utr_type = row["UTR_Type"]

        if transcript == "na" or utr_type == "na":
            continue

        pos = row["PositionVCF"]

        if pos == "-1":
            pos = row["Start"]

        try:
            pos = int(pos)
        except:
            continue

        gene = row["GeneSymbol"]

        key = (
            gene,
            utr_type,
            transcript
        )

        plp.setdefault(key, set()).add(pos)


# ============================================================
# GROUP CANONICAL VUS
# ============================================================

groups = {}

for row in canonical:

    key = (
        row["gene"],
        row["utr"],
        row["transcript"]
    )

    groups.setdefault(key, []).append(
        row["position"]
    )


# ============================================================
# MATCH GROUPS TO MANE COORDINATES
# ============================================================

group_positions = {}

multiple_coordinate_candidates = 0

for key, vus_positions in groups.items():

    gene, utr_type, transcript = key

    candidates = []

    for mane_key, positions in utr_positions.items():

        chrom, tr, ut = mane_key

        if tr == transcript and ut == utr_type:
            candidates.append(positions)

    if len(candidates) > 1:
        multiple_coordinate_candidates += 1

    if not candidates:
        continue

    # If multiple chromosome entries exist,
    # use the longest matching UTR coordinate set.
    available_utr = max(
        candidates,
        key=len
    )

    plp_positions = plp.get(
        key,
        set()
    )

    if not plp_positions:
        continue

    # Count VUS observations at each unique coordinate.
    vus_counts = {}

    for position in vus_positions:

        vus_counts[position] = (
            vus_counts.get(position, 0) + 1
        )

    unique_vus_positions = list(
        vus_counts.keys()
    )

    # Exclude positions occupied by P/LP
    # from the randomization pool.
    random_pool = [
        p for p in available_utr
        if p not in plp_positions
    ]

    n_unique = len(
        unique_vus_positions
    )

    if len(random_pool) < n_unique:

        print(
            "SKIPPED GROUP:",
            key,
            "| VUS rows =",
            len(vus_positions),
            "| unique VUS coordinates =",
            n_unique,
            "| random pool =",
            len(random_pool)
        )

        continue

    group_positions[key] = {

        "vus_counts": vus_counts,

        "unique_vus_n": n_unique,

        "available": random_pool,

        "plp": sorted(plp_positions)
    }


# ============================================================
# NEAREST DISTANCE FUNCTION
# ============================================================

def nearest_distance(
    position,
    plp_positions
):

    best = None

    for p in plp_positions:

        d = abs(
            position - p
        )

        if best is None or d < best:
            best = d

    return best


# ============================================================
# OBSERVED DISTANCES
# ============================================================

observed_distances = []

for key, positions in groups.items():

    plp_positions = plp.get(
        key,
        set()
    )

    if not plp_positions:
        continue

    for position in positions:

        distances = [

            abs(position - p)

            for p in plp_positions

            if p != position

        ]

        if distances:

            observed_distances.append(
                min(distances)
            )


# Make sure every canonical VUS has
# an observed distance.

if len(observed_distances) != len(canonical):

    raise ValueError(

        "Observed distance count does not match "
        "canonical VUS count: "

        f"{len(observed_distances)} vs "
        f"{len(canonical)}"
    )


observed_mean = statistics.mean(
    observed_distances
)


# ============================================================
# PERMUTATION TEST
# ============================================================

null_means = []

valid_permutations = 0


for iteration in range(
    N_PERMUTATIONS
):

    permutation_distances = []


    for key, info in group_positions.items():

        vus_counts = info[
            "vus_counts"
        ]

        n_unique = info[
            "unique_vus_n"
        ]

        pool = info[
            "available"
        ]

        plp_positions = info[
            "plp"
        ]


        # Sample unique genomic coordinates.
        random_positions = random.sample(
            pool,
            n_unique
        )


        # Preserve the observed multiplicity
        # distribution, but randomly assign
        # multiplicities to sampled coordinates.
        multiplicities = list(
            vus_counts.values()
        )

        random.shuffle(
            multiplicities
        )


        for position, count in zip(
            random_positions,
            multiplicities
        ):

            d = nearest_distance(
                position,
                plp_positions
            )

            if d is not None:

                for _ in range(count):

                    permutation_distances.append(
                        d
                    )


    if not permutation_distances:
        continue


    null_means.append(
        statistics.mean(
            permutation_distances
        )
    )

    valid_permutations += 1


# ============================================================
# FINAL PERMUTATION STATISTICS
# ============================================================

if not null_means:

    raise ValueError(
        "No valid permutations were generated."
    )


null_mean = statistics.mean(
    null_means
)


sorted_null = sorted(
    null_means
)


low_index = int(
    0.025 * len(sorted_null)
)


high_index = int(
    0.975 * len(sorted_null)
)


high_index = min(
    high_index,
    len(sorted_null) - 1
)


null_low = sorted_null[
    low_index
]


null_high = sorted_null[
    high_index
]


# Empirical one-sided p-value.
# We test whether observed distance
# is smaller than randomized distance.

extreme = sum(
    x <= observed_mean
    for x in null_means
)


p_value = (
    extreme + 1
) / (
    len(null_means) + 1
)


ratio = (
    null_mean /
    observed_mean
)


# ============================================================
# FINAL OUTPUT
# ============================================================

print()

print("=" * 75)
print(
    "FINAL CANONICAL PROXIMITY PERMUTATION TEST"
)
print("=" * 75)

print()

print(
    "Canonical VUS:",
    len(canonical)
)

print(
    "Canonical VUS distances:",
    len(observed_distances)
)

print(
    "Total canonical groups:",
    len(groups)
)

print(
    "Usable permutation groups:",
    len(group_positions)
)

print(
    "Groups with multiple MANE coordinate candidates:",
    multiple_coordinate_candidates
)

print(
    "Valid permutations:",
    valid_permutations
)

print()

print(
    "Observed mean nearest P/LP distance:",
    round(observed_mean, 4)
)

print(
    "Permutation mean nearest P/LP distance:",
    round(null_mean, 4)
)

print(
    "Permutation 95% interval:",
    round(null_low, 4),
    "to",
    round(null_high, 4)
)

print()

print(
    "Random / observed ratio:",
    round(ratio, 4)
)

print(
    "Empirical one-sided p-value:",
    p_value
)

print()

if p_value < 0.05:

    print(
        "RESULT: VUS are significantly closer "
        "to known P/LP variants than expected "
        "under within-group randomization."
    )

else:

    print(
        "RESULT: Proximity difference is not "
        "statistically significant."
    )

print()

print("=" * 75)
print("END")
print("=" * 75)

