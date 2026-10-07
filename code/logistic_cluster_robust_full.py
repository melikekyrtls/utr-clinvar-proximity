import csv
import math

# ============================================================
# LOAD MANE UTR LENGTHS
# ============================================================

utr_lengths = {}

with open("MANE_UTR_coordinates_v2.tsv", encoding="utf-8") as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        key = (row["Transcript"], row["UTR_Type"])
        start = int(row["Start"])
        stop = int(row["Stop"])
        length = stop - start + 1

        utr_lengths[key] = utr_lengths.get(key, 0) + length


# ============================================================
# LOAD CANONICAL DATASET
# ============================================================

data = []

with open("canonical_utr_vus_dataset.tsv", encoding="utf-8") as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:

        if row["Review_Category"] not in [
            "Multiple submitters",
            "Single submitter",
            "Expert panel",
            "Practice guideline"
        ]:
            continue

        key = (row["Transcript"], row["UTR_Type"])

        if key not in utr_lengths:
            continue

        distance = float(row["Distance"])
        length = utr_lengths[key]

        x_distance = math.log10(distance + 1)
        x_length = length / 1000.0

        # 3' UTR = reference
        # 5' UTR = 1
        x_type = 1.0 if row["UTR_Type"] == "5prime_UTR" else 0.0

        y = 1 if row["Review_Category"] == "Multiple submitters" else 0

        cluster = (
            row["Gene"],
            row["UTR_Type"],
            row["Transcript"]
        )

        data.append(
            (x_distance, x_length, x_type, y, cluster)
        )


# ============================================================
# MATRIX UTILITIES
# ============================================================

def invert_matrix(M):

    n = len(M)

    A = [
        M[i][:] + [
            1.0 if i == j else 0.0
            for j in range(n)
        ]
        for i in range(n)
    ]

    for col in range(n):

        pivot = max(
            range(col, n),
            key=lambda r: abs(A[r][col])
        )

        A[col], A[pivot] = A[pivot], A[col]

        if abs(A[col][col]) < 1e-15:
            raise RuntimeError("Singular matrix.")

        pv = A[col][col]

        for j in range(2 * n):
            A[col][j] /= pv

        for r in range(n):

            if r == col:
                continue

            factor = A[r][col]

            for j in range(2 * n):
                A[r][j] -= factor * A[col][j]

    return [
        A[i][n:]
        for i in range(n)
    ]


def mat_vec(M, v):

    return [
        sum(M[i][j] * v[j] for j in range(len(v)))
        for i in range(len(M))
    ]


# ============================================================
# LOGISTIC REGRESSION
# ============================================================

beta = [0.0, 0.0, 0.0, 0.0]

for iteration in range(100):

    g = [0.0] * 4

    H = [
        [0.0] * 4
        for _ in range(4)
    ]

    for xd, xl, xt, y, cluster in data:

        X = [1.0, xd, xl, xt]

        z = sum(
            beta[i] * X[i]
            for i in range(4)
        )

        if z >= 0:
            p = 1.0 / (1.0 + math.exp(-z))
        else:
            ez = math.exp(z)
            p = ez / (1.0 + ez)

        w = p * (1.0 - p)

        for i in range(4):

            g[i] += X[i] * (y - p)

            for j in range(4):
                H[i][j] += w * X[i] * X[j]

    H_inv = invert_matrix(H)

    step = mat_vec(H_inv, g)

    beta_new = [
        beta[i] + step[i]
        for i in range(4)
    ]

    beta = beta_new

    if max(abs(x) for x in step) < 1e-10:
        break


# ============================================================
# MODEL-BASED INFORMATION MATRIX
# ============================================================

H = [
    [0.0] * 4
    for _ in range(4)
]

# Cluster score vectors
cluster_scores = {}

for xd, xl, xt, y, cluster in data:

    X = [1.0, xd, xl, xt]

    z = sum(
        beta[i] * X[i]
        for i in range(4)
    )

    if z >= 0:
        p = 1.0 / (1.0 + math.exp(-z))
    else:
        ez = math.exp(z)
        p = ez / (1.0 + ez)

    w = p * (1.0 - p)

    for i in range(4):
        for j in range(4):
            H[i][j] += w * X[i] * X[j]

    score = [
        X[i] * (y - p)
        for i in range(4)
    ]

    if cluster not in cluster_scores:
        cluster_scores[cluster] = [0.0] * 4

    for i in range(4):
        cluster_scores[cluster][i] += score[i]


# ============================================================
# CLUSTER-ROBUST SANDWICH VARIANCE
# ============================================================

bread = invert_matrix(H)

meat = [
    [0.0] * 4
    for _ in range(4)
]

for score in cluster_scores.values():

    for i in range(4):
        for j in range(4):
            meat[i][j] += score[i] * score[j]


# V = H^-1 * meat * H^-1

temp = [
    [
        sum(
            bread[i][k] * meat[k][j]
            for k in range(4)
        )
        for j in range(4)
    ]
    for i in range(4)
]

V = [
    [
        sum(
            temp[i][k] * bread[k][j]
            for k in range(4)
        )
        for j in range(4)
    ]
    for i in range(4)
]


se = [
    math.sqrt(max(V[i][i], 0.0))
    for i in range(4)
]


# ============================================================
# OUTPUT
# ============================================================

names = [
    "Intercept",
    "Distance",
    "UTR_length_per_1000bp",
    "5prime_UTR_vs_3prime"
]

print()
print("=" * 80)
print("CLUSTER-ROBUST LOGISTIC REGRESSION")
print("=" * 80)

print()
print(
    "Model:"
)

print(
    "Multiple-submitter review ~ "
    "log10(distance + 1) + UTR length + UTR type"
)

print()
print("N:", len(data))
print("Clusters:", len(cluster_scores))
print(
    "Multiple-submitter:",
    sum(row[3] for row in data)
)

print()

for i in range(4):

    z = beta[i] / se[i]

    p_value = math.erfc(
        abs(z) / math.sqrt(2)
    )

    odds_ratio = math.exp(beta[i])

    ci_low = math.exp(
        beta[i] - 1.96 * se[i]
    )

    ci_high = math.exp(
        beta[i] + 1.96 * se[i]
    )

    print(names[i])

    print(
        "  coefficient:",
        round(beta[i], 6)
    )

    print(
        "  cluster-robust SE:",
        round(se[i], 6)
    )

    print(
        "  Wald Z:",
        round(z, 6)
    )

    print(
        "  p:",
        p_value
    )

    print(
        "  OR:",
        round(odds_ratio, 6)
    )

    print(
        "  95% CI:",
        round(ci_low, 6),
        "to",
        round(ci_high, 6)
    )

    print()


print("=" * 80)

if beta[1] < 0:
    print(
        "RESULT: The negative distance association persists "
        "after adjustment for UTR length and UTR type, "
        "with gene-UTR-transcript cluster-robust SEs."
    )
else:
    print(
        "RESULT: The distance association is not negative "
        "after full adjustment."
    )

print("=" * 80)

