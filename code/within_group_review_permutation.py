import pandas as pd
import numpy as np

BASE = r"..\data"

N_PERM = 5000
SEED = 42

df = pd.read_csv(
    BASE + r"\canonical_utr_vus_dataset.tsv",
    sep="\t"
)

# ------------------------------------------------------------
# DATA CLEANING
# ------------------------------------------------------------

df["Distance"] = pd.to_numeric(
    df["Distance"],
    errors="coerce"
)

df = df.dropna(
    subset=[
        "Distance",
        "Review_Category",
        "Gene",
        "UTR_Type",
        "Transcript"
    ]
).copy()

# Same transformation used in the regression
df["log_distance"] = np.log10(
    df["Distance"] + 1
)

# Multiple submitter = 1
df["Multiple"] = (
    df["Review_Category"]
    .str.lower()
    == "multiple submitters"
).astype(int)

group_cols = [
    "Gene",
    "UTR_Type",
    "Transcript"
]

# ------------------------------------------------------------
# OBSERVED STATISTIC
#
# Difference:
# mean log-distance among Multiple submitter
# minus
# mean log-distance among Single submitter
#
# Negative = Multiple submitter variants are closer.
# ------------------------------------------------------------

x = df["log_distance"].to_numpy(dtype=float)
y = df["Multiple"].to_numpy(dtype=int)

observed_multiple_mean = x[y == 1].mean()
observed_single_mean = x[y == 0].mean()

observed_difference = (
    observed_multiple_mean
    - observed_single_mean
)

print("=" * 75)
print("WITHIN-GROUP REVIEW-STATUS PERMUTATION TEST")
print("=" * 75)

print("Rows:", len(df))
print(
    "Multiple submitter:",
    int(y.sum())
)
print(
    "Single submitter:",
    int((y == 0).sum())
)

print(
    "\nObserved mean log-distance, Multiple:",
    observed_multiple_mean
)

print(
    "Observed mean log-distance, Single:",
    observed_single_mean
)

print(
    "Observed difference (Multiple - Single):",
    observed_difference
)

print(
    "\nNegative difference means Multiple-submitter variants "
    "are closer to P/LP variants."
)

# ------------------------------------------------------------
# PREPARE GROUPS
#
# Review-status labels will be shuffled ONLY within the
# same gene-UTR-transcript group.
# Thus each group's number of Multiple and Single records
# remains exactly unchanged.
# ------------------------------------------------------------

groups = []

for key, group in df.groupby(group_cols):

    gx = group["log_distance"].to_numpy(dtype=float)
    gy = group["Multiple"].to_numpy(dtype=int)

    n_multiple = int(gy.sum())
    n_total = len(gy)

    # Groups with no variation cannot contribute to a
    # within-group comparison.
    if n_multiple == 0 or n_multiple == n_total:
        continue

    groups.append(
        (
            gx,
            n_multiple,
            n_total
        )
    )

print(
    "Groups with both review categories:",
    len(groups)
)

# ------------------------------------------------------------
# PERMUTATION
# ------------------------------------------------------------

rng = np.random.default_rng(SEED)

null_statistics = np.empty(
    N_PERM,
    dtype=float
)

for perm in range(N_PERM):

    total_multiple_sum = 0.0
    total_multiple_n = 0
    total_single_sum = 0.0
    total_single_n = 0

    for gx, n_multiple, n_total in groups:

        # Randomly choose which positions receive the
        # fixed number of Multiple-submitter labels.
        chosen = rng.choice(
            n_total,
            size=n_multiple,
            replace=False
        )

        mask = np.zeros(
            n_total,
            dtype=bool
        )

        mask[chosen] = True

        total_multiple_sum += gx[mask].sum()
        total_multiple_n += n_multiple

        total_single_sum += gx[~mask].sum()
        total_single_n += n_total - n_multiple

    perm_multiple_mean = (
        total_multiple_sum
        / total_multiple_n
    )

    perm_single_mean = (
        total_single_sum
        / total_single_n
    )

    null_statistics[perm] = (
        perm_multiple_mean
        - perm_single_mean
    )

# ------------------------------------------------------------
# RESULTS
# ------------------------------------------------------------

null_mean = null_statistics.mean()

null_low = np.percentile(
    null_statistics,
    2.5
)

null_high = np.percentile(
    null_statistics,
    97.5
)

null_sd = null_statistics.std(
    ddof=1
)

# One-sided test:
# H1 = observed difference is MORE NEGATIVE
# than expected under within-group label randomization.

extreme = np.sum(
    null_statistics <= observed_difference
)

p_value = (
    extreme + 1
) / (
    N_PERM + 1
)

permutation_z = (
    observed_difference
    - null_mean
) / null_sd

print("\n" + "=" * 75)
print("RESULTS")
print("=" * 75)

print(
    "Observed difference:",
    f"{observed_difference:.8f}"
)

print(
    "Permutation mean:",
    f"{null_mean:.8f}"
)

print(
    "Permutation 95% interval:",
    f"{null_low:.8f} to {null_high:.8f}"
)

print(
    "Permutation SD:",
    f"{null_sd:.8f}"
)

print(
    "Permutation Z:",
    f"{permutation_z:.4f}"
)

print(
    "Extreme permutations:",
    int(extreme)
)

print(
    "One-sided empirical p-value:",
    p_value
)

print(
    "\nInterpretation:"
)

if p_value < 0.05:
    print(
        "The observed review-status/distance association "
        "is stronger than expected after preserving "
        "group-specific review-status counts."
    )
else:
    print(
        "The observed review-status/distance association "
        "is not stronger than expected after preserving "
        "group-specific review-status counts."
    )

print("\n" + "=" * 75)
print("DONE")
print("=" * 75)
