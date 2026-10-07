import csv
import random
import statistics

random.seed(42)

# ============================================================
# FILES
# ============================================================

CLINVAR = "ClinVar_UTR_GRCh38.tsv"
MANE = "MANE_UTR_coordinates_v2.tsv"
OVERLAP = "canonical_utr_overlapping_other_transcript_cds.tsv"
CANONICAL = "canonical_utr_vus_dataset.tsv"

N_PERMUTATIONS = 5000


# ============================================================
# LOAD CDS-OVERLAP INTERVALS
# ============================================================

overlap_intervals = {}

with open(OVERLAP, encoding="utf-8") as f:

    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:

        key = (
            row["Gene"],
            row["UTR_Type"],
            row["Transcript"]
        )

        interval = (
            int(row["Overlap_Start"]),
            int(row["Overlap_Stop"])
        )

        overlap_intervals.setdefault(
            key,
            []
        ).append(interval)


# ============================================================
# MERGE OVERLAP INTERVALS
# ============================================================

for key in overlap_intervals:

    intervals = sorted(
        overlap_intervals[key]
    )

    merged = []

    for start, stop in intervals:

        if not merged:

            merged.append(
                [start, stop]
            )

        elif start <= merged[-1][1] + 1:

            merged[-1][1] = max(
                merged[-1][1],
                stop
            )

        else:

            merged.append(
                [start, stop]
            )

    overlap_intervals[key] = merged


def is_in_overlap(key, position):

    for start, stop in overlap_intervals.get(
        key,
        []
    ):

        if start <= position <= stop:
            return True

    return False


# ============================================================
# LOAD CANONICAL VUS
# ============================================================

canonical = []
removed_vus = []

with open(CANONICAL, encoding="utf-8") as f:

    reader = csv.DictReader(
        f,
        delimiter="\t"
    )

    for row in reader:

        key = (
            row["Gene"],
            row["UTR_Type"],
            row["Transcript"]
        )

        position = int(
            row["Position"]
        )

        if is_in_overlap(
            key,
            position
        ):

            removed_vus.append(row)

        else:

            canonical.append({

                "gene": row["Gene"],

                "utr": row["UTR_Type"],

                "transcript": row["Transcript"],

                "position": position,

                "distance": float(
                    row["Distance"]
                )
            })


# ============================================================
# LOAD MANE UTR SEGMENTS
# ============================================================

segments = {}

with open(MANE, encoding="utf-8") as f:

    reader = csv.DictReader(
        f,
        delimiter="\t"
    )

    for row in reader:

        key = (
            row["Chromosome"],
            row["Transcript"],
            row["UTR_Type"]
        )

        segments.setdefault(
            key,
            []
        ).append(
            (
                int(row["Start"]),
                int(row["Stop"])
            )
        )


# ============================================================
# MERGE MANE UTR SEGMENTS
# ============================================================

merged_segments = {}

for key, intervals in segments.items():

    intervals = sorted(
        intervals
    )

    merged = []

    for start, stop in intervals:

        if not merged:

            merged.append(
                [start, stop]
            )

        elif start <= merged[-1][1] + 1:

            merged[-1][1] = max(
                merged[-1][1],
                stop
            )

        else:

            merged.append(
                [start, stop]
            )

    merged_segments[key] = merged


# ============================================================
# GENERATE UTR POSITIONS
# ============================================================

utr_positions = {}

for key, intervals in merged_segments.items():

    positions = []

    for start, stop in intervals:

        positions.extend(
            range(
                start,
                stop + 1
            )
        )

    utr_positions[key] = positions


# ============================================================
# LOAD P/LP POSITIONS
# ============================================================

plp = {}

with open(CLINVAR, encoding="utf-8") as f:

    reader = csv.DictReader(
        f,
        delimiter="\t"
    )

    for row in reader:

        if row["Assembly"] != "GRCh38":
            continue

        significance = row[
            "ClinicalSignificance"
        ]

        if significance not in [

            "Pathogenic",

            "Likely pathogenic",

            "Pathogenic/Likely pathogenic",

            "Likely pathogenic/Pathogenic"

        ]:

            continue

        transcript = row[
            "UTR_Transcript"
        ]

        utr_type = row[
            "UTR_Type"
        ]

        if (
            transcript == "na"
            or
            utr_type == "na"
        ):

            continue

        pos = row[
            "PositionVCF"
        ]

        if pos == "-1":

            pos = row["Start"]

        try:

            pos = int(pos)

        except:

            continue

        gene = row[
            "GeneSymbol"
        ]

        key = (
            gene,
            utr_type,
            transcript
        )

        # Exclude P/LP variants that themselves
        # fall inside another transcript's CDS.

        if is_in_overlap(
            key,
            pos
        ):

            continue

        plp.setdefault(
            key,
            set()
        ).add(pos)


# ============================================================
# GROUP REMAINING VUS
# ============================================================

groups = {}

for row in canonical:

    key = (
        row["gene"],
        row["utr"],
        row["transcript"]
    )

    groups.setdefault(
        key,
        []
    ).append(
        row["position"]
    )


# ============================================================
# BUILD USABLE PERMUTATION GROUPS
# ============================================================

group_positions = {}

skipped_no_plp = 0
skipped_no_mane = 0
skipped_small_pool = 0

for key, vus_positions in groups.items():

    gene, utr_type, transcript = key

    candidates = []

    for mane_key, positions in utr_positions.items():

        chrom, tr, ut = mane_key

        if (
            tr == transcript
            and
            ut == utr_type
        ):

            candidates.append(
                positions
            )

    if not candidates:

        skipped_no_mane += 1

        continue

    available_utr = max(
        candidates,
        key=len
    )

    plp_positions = plp.get(
        key,
        set()
    )

    if not plp_positions:

        skipped_no_plp += 1

        continue

    # Count VUS observations at each coordinate.

    vus_counts = {}

    for position in vus_positions:

        vus_counts[position] = (
            vus_counts.get(
                position,
                0
            ) + 1
        )

    unique_vus_positions = list(
        vus_counts.keys()
    )

    # Randomization pool:
    # - inside the same UTR
    # - not a P/LP coordinate
    # - not inside another transcript's CDS overlap

    random_pool = [

        p

        for p in available_utr

        if p not in plp_positions

        and not is_in_overlap(
            key,
            p
        )

    ]

    n_unique = len(
        unique_vus_positions
    )

    if len(random_pool) < n_unique:

        skipped_small_pool += 1

        continue

    group_positions[key] = {

        "vus_counts": vus_counts,

        "unique_vus_n": n_unique,

        "available": random_pool,

        "plp": sorted(
            plp_positions
        )
    }


# ============================================================
# DEFINE ANALYSIS POPULATION
# ============================================================

analysis_canonical = []

excluded_no_plp = []

for row in canonical:

    key = (
        row["gene"],
        row["utr"],
        row["transcript"]
    )

    if key in group_positions:

        analysis_canonical.append(
            row
        )

    else:

        excluded_no_plp.append(
            row
        )


# ============================================================
# OBSERVED DISTANCES
# ============================================================

observed_distances = []

for row in analysis_canonical:

    key = (
        row["gene"],
        row["utr"],
        row["transcript"]
    )

    position = row["position"]

    plp_positions = plp[key]

    distances = [

        abs(
            position - p
        )

        for p in plp_positions

        if p != position

    ]

    if distances:

        observed_distances.append(
            min(distances)
        )


# ============================================================
# SAFETY CHECK
# ============================================================

if len(observed_distances) != len(
    analysis_canonical
):

    raise ValueError(

        "Analysis distance count mismatch: "

        f"{len(observed_distances)} vs "

        f"{len(analysis_canonical)}"
    )


observed_mean = statistics.mean(
    observed_distances
)


# ============================================================
# PERMUTATION TEST
# ============================================================

null_means = []

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

        # Sample unique positions.

        random_positions = random.sample(

            pool,

            n_unique
        )

        # Preserve observed multiplicity
        # distribution.

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

            distances = [

                abs(
                    position - p
                )

                for p in plp_positions

            ]

            if distances:

                d = min(
                    distances
                )

                for _ in range(
                    count
                ):

                    permutation_distances.append(
                        d
                    )

    if permutation_distances:

        null_means.append(

            statistics.mean(
                permutation_distances
            )
        )


# ============================================================
# FINAL STATISTICS
# ============================================================

if not null_means:

    raise ValueError(
        "No valid permutations."
    )


null_mean = statistics.mean(
    null_means
)


sorted_null = sorted(
    null_means
)


low_index = int(
    0.025 *
    len(sorted_null)
)


high_index = min(

    int(
        0.975 *
        len(sorted_null)
    ),

    len(sorted_null) - 1
)


null_low = sorted_null[
    low_index
]


null_high = sorted_null[
    high_index
]


# One-sided test:
# observed distance <= randomized distance

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
    "CORRECTED CDS-OVERLAP "
    "SENSITIVITY ANALYSIS"
)

print("=" * 75)

print()

print(
    "Original canonical VUS:",
    len(canonical)
)

print(
    "VUS removed inside CDS overlap:",
    len(removed_vus)
)

print(
    "VUS remaining after positional filter:",
    len(canonical)
)

print(
    "VUS retained for analysis:",
    len(analysis_canonical)
)

print(
    "VUS excluded because no usable P/LP remained:",
    len(excluded_no_plp)
)

print()

print(
    "Removed VUS %:",
    round(
        100 *
        len(removed_vus) /
        (
            len(canonical)
            +
            len(removed_vus)
        ),
        4
    )
)

print()

print(
    "Original groups:",
    len(groups)
)

print(
    "Usable permutation groups:",
    len(group_positions)
)

print(
    "Groups skipped - no MANE:",
    skipped_no_mane
)

print(
    "Groups skipped - no remaining P/LP:",
    skipped_no_plp
)

print(
    "Groups skipped - insufficient random pool:",
    skipped_small_pool
)

print()

print(
    "P/LP groups after overlap filtering:",
    len(plp)
)

print(
    "Valid permutations:",
    len(null_means)
)

print()

print(
    "Observed mean nearest P/LP distance:",
    round(
        observed_mean,
        4
    )
)

print(
    "Permutation mean:",
    round(
        null_mean,
        4
    )
)

print(
    "Permutation 95% interval:",
    round(
        null_low,
        4
    ),
    "to",
    round(
        null_high,
        4
    )
)

print()

print(
    "Random / observed ratio:",
    round(
        ratio,
        4
    )
)

print(
    "Empirical one-sided p-value:",
    p_value
)

print()

if p_value < 0.05:

    print(
        "RESULT: Proximity remains significant "
        "after excluding UTR positions overlapping "
        "another transcript's CDS."
    )

else:

    print(
        "RESULT: Proximity is not significant "
        "after excluding UTR positions overlapping "
        "another transcript's CDS."
    )

print()

print("=" * 75)

print("END")

print("=" * 75)
