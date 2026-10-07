import csv
import random
import math

random.seed(42)

# ============================================================
# LOAD CANONICAL DATASET
# ============================================================

rows = []

with open("canonical_utr_vus_dataset.tsv", encoding="utf-8") as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        rows.append({
            "gene": row["Gene"],
            "utr": row["UTR_Type"],
            "transcript": row["Transcript"],
            "distance": float(row["Distance"]),
            "bin": row["Distance_Bin"],
            "review": row["Review_Category"]
        })

# ============================================================
# DISTANCE ORDER
# ============================================================

bin_order = [
    "0-5",
    "6-10",
    "11-25",
    "26-50",
    "51-100",
    ">100"
]

# ============================================================
# CALCULATE TREND Z
# ============================================================

def trend_z(data):

    scores = {
        "0-5": 1,
        "6-10": 2,
        "11-25": 3,
        "26-50": 4,
        "51-100": 5,
        ">100": 6
    }

    # Multiple submitter = success
    total = 0
    multiple = 0

    group_counts = {}

    for row in data:

        if row["bin"] not in scores:
            continue

        if row["review"] not in [
            "Multiple submitters",
            "Single submitter",
            "Expert panel",
            "Practice guideline"
        ]:
            continue

        b = row["bin"]

        if b not in group_counts:
            group_counts[b] = {
                "n": 0,
                "multiple": 0
            }

        group_counts[b]["n"] += 1

        if row["review"] == "Multiple submitters":
            group_counts[b]["multiple"] += 1

        total += 1

    if total == 0:
        return None

    # Cochran-Armitage style ordered score test
    weighted = 0

    for b in bin_order:
        if b in group_counts:
            weighted += (
                scores[b] *
                group_counts[b]["multiple"]
            )

    p = sum(
        group_counts[b]["multiple"]
        for b in group_counts
    )

    n = total

    mean_score = sum(
        scores[b] * group_counts[b]["n"]
        for b in group_counts
    ) / n

    variance_score = sum(
        group_counts[b]["n"] *
        (scores[b] - mean_score) ** 2
        for b in group_counts
    ) / n

    if variance_score == 0:
        return 0

    observed_score = weighted - p * mean_score

    denominator = math.sqrt(
        p * (1 - p / n) * variance_score
    )

    if denominator == 0:
        return 0

    return observed_score / denominator

# ============================================================
# OBSERVED
# ============================================================

observed_z = trend_z(rows)

# ============================================================
# GROUP DATA
# ============================================================

groups = {}

for row in rows:

    key = (
        row["gene"],
        row["utr"],
        row["transcript"]
    )

    if key not in groups:
        groups[key] = []

    groups[key].append(row)

# ============================================================
# WITHIN-GROUP PERMUTATION
# ============================================================

N_PERMUTATIONS = 5000

permutation_z = []

for iteration in range(N_PERMUTATIONS):

    permuted = []

    for key, group_rows in groups.items():

        # Preserve the exact review labels/counts
        labels = [
            row["review"]
            for row in group_rows
        ]

        random.shuffle(labels)

        for row, label in zip(group_rows, labels):

            permuted.append({
                "bin": row["bin"],
                "review": label
            })

    z = trend_z(permuted)

    if z is not None:
        permutation_z.append(z)

# ============================================================
# EMPIRICAL P-VALUE
# ============================================================

# We test the observed negative trend.
extreme = sum(
    z <= observed_z
    for z in permutation_z
)

p_value = (
    extreme + 1
) / (
    len(permutation_z) + 1
)

# ============================================================
# NULL DISTRIBUTION SUMMARY
# ============================================================

mean_null = (
    sum(permutation_z) /
    len(permutation_z)
)

sorted_z = sorted(permutation_z)

low = sorted_z[
    int(0.025 * len(sorted_z))
]

high = sorted_z[
    int(0.975 * len(sorted_z))
]

# ============================================================
# OUTPUT
# ============================================================

print()
print("=" * 75)
print("WITHIN-GROUP REVIEW-STATUS PERMUTATION TEST")
print("=" * 75)

print()
print("VUS:", len(rows))
print("Groups:", len(groups))
print("Permutations:", len(permutation_z))

print()
print("Observed trend Z:", round(observed_z, 6))

print("Null mean Z:", round(mean_null, 6))

print(
    "Null 95% interval:",
    round(low, 6),
    "to",
    round(high, 6)
)

print()
print("Empirical one-sided p-value:", p_value)

print()

if observed_z < 0:
    print("Observed direction: NEGATIVE")
else:
    print("Observed direction: POSITIVE")

if p_value < 0.05:
    print()
    print(
        "RESULT: The negative review-status trend "
        "persists after within-group permutation."
    )
else:
    print()
    print(
        "RESULT: The trend is not significant "
        "after within-group permutation."
    )

