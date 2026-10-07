import csv
import math

# ============================================================
# LOAD MANE UTR LENGTHS
# ============================================================

utr_lengths = {}

with open("MANE_UTR_coordinates_v2.tsv", encoding="utf-8") as f:

    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:

        key = (
            row["Transcript"],
            row["UTR_Type"]
        )

        start = int(row["Start"])
        stop = int(row["Stop"])

        length = stop - start + 1

        utr_lengths[key] = (
            utr_lengths.get(key, 0) + length
        )


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

        key = (
            row["Transcript"],
            row["UTR_Type"]
        )

        if key not in utr_lengths:
            continue

        distance = float(row["Distance"])
        utr_length = utr_lengths[key]

        x_distance = math.log10(distance + 1)

        # Per 1000 bp
        x_length = utr_length / 1000.0

        # 3' UTR = reference (0)
        # 5' UTR = 1
        x_type = (
            1.0
            if row["UTR_Type"] == "5prime_UTR"
            else 0.0
        )

        y = (
            1
            if row["Review_Category"] == "Multiple submitters"
            else 0
        )

        data.append(
            (x_distance, x_length, x_type, y)
        )


# ============================================================
# LOGISTIC REGRESSION
#
# y ~ intercept
#     + distance
#     + UTR length
#     + UTR type
# ============================================================

def sigmoid(z):

    if z >= 0:
        return 1 / (1 + math.exp(-z))

    ez = math.exp(z)

    return ez / (1 + ez)


beta = [0.0, 0.0, 0.0, 0.0]


# ============================================================
# NEWTON-RAPHSON
# ============================================================

for iteration in range(100):

    g = [0.0] * 4

    H = [
        [0.0] * 4
        for _ in range(4)
    ]

    for xd, xl, xt, y in data:

        X = [
            1.0,
            xd,
            xl,
            xt
        ]

        z = sum(
            beta[i] * X[i]
            for i in range(4)
        )

        p = sigmoid(z)

        w = p * (1 - p)

        for i in range(4):

            g[i] += X[i] * (y - p)

            for j in range(4):

                H[i][j] += (
                    w * X[i] * X[j]
                )


    # --------------------------------------------------------
    # Solve H * step = gradient
    # Gaussian elimination
    # --------------------------------------------------------

    A = [
        H[i][:] + [g[i]]
        for i in range(4)
    ]

    for col in range(4):

        pivot = max(
            range(col, 4),
            key=lambda r: abs(A[r][col])
        )

        A[col], A[pivot] = (
            A[pivot],
            A[col]
        )

        if abs(A[col][col]) < 1e-15:
            raise RuntimeError(
                "Singular Hessian matrix."
            )

        pivot_value = A[col][col]

        for j in range(col, 5):
            A[col][j] /= pivot_value

        for r in range(4):

            if r == col:
                continue

            factor = A[r][col]

            for j in range(col, 5):
                A[r][j] -= (
                    factor * A[col][j]
                )

    step = [
        A[i][4]
        for i in range(4)
    ]

    beta = [
        beta[i] + step[i]
        for i in range(4)
    ]

    if max(abs(x) for x in step) < 1e-10:
        break


# ============================================================
# INVERT FINAL HESSIAN
# ============================================================

H = [
    [0.0] * 4
    for _ in range(4)
]

for xd, xl, xt, y in data:

    X = [
        1.0,
        xd,
        xl,
        xt
    ]

    z = sum(
        beta[i] * X[i]
        for i in range(4)
    )

    p = sigmoid(z)

    w = p * (1 - p)

    for i in range(4):

        for j in range(4):

            H[i][j] += (
                w * X[i] * X[j]
            )


# Gauss-Jordan inverse

A = []

for i in range(4):

    A.append(
        H[i][:]
        + [
            1.0 if i == j else 0.0
            for j in range(4)
        ]
    )


for col in range(4):

    pivot = max(
        range(col, 4),
        key=lambda r: abs(A[r][col])
    )

    A[col], A[pivot] = (
        A[pivot],
        A[col]
    )

    pivot_value = A[col][col]

    for j in range(8):
        A[col][j] /= pivot_value

    for r in range(4):

        if r == col:
            continue

        factor = A[r][col]

        for j in range(8):

            A[r][j] -= (
                factor * A[col][j]
            )


cov = [
    A[i][4:8]
    for i in range(4)
]

se = [
    math.sqrt(cov[i][i])
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
print("=" * 75)
print("LOGISTIC REGRESSION WITH UTR TYPE")
print("=" * 75)

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
        "  SE:",
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


print("=" * 75)

if beta[1] < 0:

    print(
        "RESULT: Distance remains negatively associated "
        "with multiple-submitter review after adjusting "
        "for UTR length and UTR type."
    )

else:

    print(
        "RESULT: Distance is not negatively associated "
        "after adjustment."
    )

print("=" * 75)

