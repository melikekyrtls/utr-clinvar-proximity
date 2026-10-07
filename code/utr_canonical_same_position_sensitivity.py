import csv
import math
from collections import defaultdict

CLINVAR = "ClinVar_UTR_GRCh38.tsv"
MANE = "MANE_UTR_coordinates_v2.tsv"
CANONICAL = "canonical_utr_vus_dataset.tsv"


def is_vus(sig):
    return sig.lower().strip() == "uncertain significance"


def is_plp(sig):
    return sig.lower().strip() in {
        "pathogenic",
        "likely pathogenic",
        "pathogenic/likely pathogenic",
        "likely pathogenic/pathogenic"
    }


def get_position(row):
    try:
        pos = int(row["PositionVCF"])
        if pos >= 0:
            return pos
    except:
        pass

    try:
        return int(row["Start"])
    except:
        return None


# ------------------------------------------------------------
# 1. Load canonical VUS population
# ------------------------------------------------------------

canonical_keys = set()

with open(CANONICAL, "r", encoding="utf-8", newline="") as f:

    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:

        key = (
            row["Gene"],
            row["UTR_Type"],
            row["Transcript"],
            int(row["Position"])
        )

        canonical_keys.add(key)


print("=" * 70)
print("CANONICAL SAME-POSITION SENSITIVITY ANALYSIS")
print("=" * 70)

print("Canonical unique VUS coordinates:", len(canonical_keys))


# ------------------------------------------------------------
# 2. Load MANE UTR intervals
# ------------------------------------------------------------

mane = defaultdict(list)

with open(MANE, "r", encoding="utf-8", newline="") as f:

    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:

        chrom = row["Chromosome"].strip()

        if not chrom.startswith("chr"):
            chrom = "chr" + chrom

        key = (
            chrom,
            row["Transcript"].strip(),
            row["UTR_Type"].strip()
        )

        mane[key].append(
            (
                int(row["Start"]),
                int(row["Stop"])
            )
        )


# ------------------------------------------------------------
# 3. Load ClinVar records
# ------------------------------------------------------------

records = defaultdict(lambda: {"vus": [], "plp": []})

with open(CLINVAR, "r", encoding="utf-8", newline="") as f:

    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:

        if row.get("Assembly", "").strip() != "GRCh38":
            continue

        sig = row.get("ClinicalSignificance", "").strip()

        if not (is_vus(sig) or is_plp(sig)):
            continue

        gene = row.get("GeneSymbol", "").strip()
        transcript = row.get("UTR_Transcript", "").strip()
        utr_type = row.get("UTR_Type", "").strip()

        if not gene or not transcript or not utr_type:
            continue

        position = get_position(row)

        if position is None:
            continue

        chrom = row.get("Chromosome", "").strip()

        if not chrom.startswith("chr"):
            chrom = "chr" + chrom

        # Must map to MANE exactly like canonical analysis
        mane_key = (
            chrom,
            transcript,
            utr_type
        )

        matched = False

        for start, stop in mane.get(mane_key, []):

            if start <= position <= stop:
                matched = True
                break

        if not matched:
            continue

        record = {
            "position": position,
            "ref": row.get("ReferenceAlleleVCF", "").strip(),
            "alt": row.get("AlternateAlleleVCF", "").strip(),
            "variation_id": row.get("VariationID", "").strip(),
            "review": row.get("ReviewStatus", "").lower()
        }

        key = (
            gene,
            utr_type,
            transcript
        )

        if is_vus(sig):
            records[key]["vus"].append(record)

        else:
            records[key]["plp"].append(record)


# ------------------------------------------------------------
# 4. Keep ONLY canonical VUS coordinates
# ------------------------------------------------------------

results = []

same_position_cases = 0

for key, data in records.items():

    gene, utr_type, transcript = key

    if not data["vus"] or not data["plp"]:
        continue

    for vus in data["vus"]:

        canonical_key = (
            gene,
            utr_type,
            transcript,
            vus["position"]
        )

        if canonical_key not in canonical_keys:
            continue

        distances = []

        for plp in data["plp"]:

            vus_allele = (
                vus["ref"],
                vus["alt"]
            )

            plp_allele = (
                plp["ref"],
                plp["alt"]
            )

            # Same exact allele:
            # do not compare it against itself.
            if (
                vus["position"] == plp["position"]
                and vus_allele == plp_allele
            ):
                continue

            # Different allele at same genomic position:
            # sensitivity analysis treats this as 0 bp.
            if (
                vus["position"] == plp["position"]
                and vus_allele != plp_allele
            ):
                distance = 0
                same_position_cases += 1

            else:
                distance = abs(
                    vus["position"] - plp["position"]
                )

            distances.append(distance)

        if not distances:
            continue

        nearest = min(distances)

        multiple = (
            "multiple submitters" in vus["review"]
        )

        results.append({
            "distance": nearest,
            "multiple": multiple
        })


# ------------------------------------------------------------
# 5. Ordered trend test
# ------------------------------------------------------------

bins = [
    (0, 5),
    (6, 10),
    (11, 25),
    (26, 50),
    (51, 100),
    (101, 10**18)
]

scores = [1, 2, 3, 4, 5, 6]

totals = [0] * 6
multiples = [0] * 6

for r in results:

    d = r["distance"]

    for i, (lo, hi) in enumerate(bins):

        if lo <= d <= hi:

            totals[i] += 1

            if r["multiple"]:
                multiples[i] += 1

            break


n = sum(totals)
total_multiple = sum(multiples)

overall_rate = total_multiple / n

weighted_score = sum(
    scores[i] * totals[i]
    for i in range(6)
)

weighted_multiple = sum(
    scores[i] * multiples[i]
    for i in range(6)
)

numerator = (
    weighted_multiple
    - overall_rate * weighted_score
)

mean_score = weighted_score / n

score_variance = sum(
    totals[i] * (scores[i] - mean_score) ** 2
    for i in range(6)
)

variance = (
    overall_rate
    * (1 - overall_rate)
    * score_variance
)

z = numerator / math.sqrt(variance)

p = math.erfc(abs(z) / math.sqrt(2))


# ------------------------------------------------------------
# 6. Output
# ------------------------------------------------------------

print()
print("Canonical population retained:", len(canonical_keys))
print("VUS analyzed:", n)

print(
    "Same-position different-allele cases:",
    same_position_cases
)

print()

print(
    "Multiple-submitter:",
    total_multiple,
    "/",
    n,
    "=",
    round(100 * total_multiple / n, 4),
    "%"
)

print()
print("DISTANCE BINS")
print("-" * 70)

labels = [
    "0-5",
    "6-10",
    "11-25",
    "26-50",
    "51-100",
    ">100"
]

for label, total, mult in zip(
    labels,
    totals,
    multiples
):

    rate = (
        100 * mult / total
        if total
        else 0
    )

    print(
        f"{label:8s} "
        f"{mult:5d}/{total:<5d} "
        f"= {rate:.2f}%"
    )

print()
print("ORDERED TREND TEST")
print("-" * 70)

print("Z =", round(z, 6))
print("Two-sided p =", f"{p:.6e}")

if z < 0:
    print(
        "RESULT: Negative trend preserved."
    )
else:
    print(
        "RESULT: Negative trend NOT preserved."
    )

print()
print("=" * 70)
print("AUDIT COMPLETE")
print("=" * 70)
