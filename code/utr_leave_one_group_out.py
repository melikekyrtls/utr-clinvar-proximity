import csv
from collections import defaultdict

# ============================================================
# LOAD CLINVAR
# ============================================================

data = []

with open("ClinVar_UTR_GRCh38.tsv", encoding="utf-8") as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        cs = row["ClinicalSignificance"].strip()

        if cs in [
            "Uncertain significance",
            "Uncertain significance/Uncertain risk allele"
        ]:
            category = "VUS"

        elif cs in [
            "Pathogenic",
            "Likely pathogenic",
            "Pathogenic/Likely pathogenic"
        ]:
            category = "P/LP"

        else:
            continue

        try:
            pos = int(row["PositionVCF"])
        except:
            continue

        data.append({
            "gene": row["GeneSymbol"],
            "utr": row["UTR_Type"],
            "pos": pos,
            "category": category,
            "transcript": row["UTR_Transcript"]
        })


# ============================================================
# LOAD MANE UTR
# ============================================================

segments = defaultdict(list)

with open("MANE_UTR_coordinates_v2.tsv", encoding="utf-8") as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        try:
            start = int(row["Start"])
            stop = int(row["Stop"])
        except:
            continue

        key = (row["Transcript"], row["UTR_Type"])

        segments[key].append({
            "start": start,
            "stop": stop,
            "strand": row["Strand"]
        })


# ============================================================
# BUILD TRANSCRIPT COORDINATE MAP
# ============================================================

def build_coordinate_map(segment_list):

    strand = segment_list[0]["strand"]

    if strand == "+":
        ordered = sorted(segment_list, key=lambda x: x["start"])
    else:
        ordered = sorted(
            segment_list,
            key=lambda x: x["start"],
            reverse=True
        )

    coordinate_map = {}
    transcript_pos = 1

    for seg in ordered:

        if strand == "+":
            positions = range(seg["start"], seg["stop"] + 1)
        else:
            positions = range(seg["stop"], seg["start"] - 1, -1)

        for genomic_pos in positions:
            coordinate_map[genomic_pos] = transcript_pos
            transcript_pos += 1

    return coordinate_map


coordinate_maps = {
    key: build_coordinate_map(segs)
    for key, segs in segments.items()
}


# ============================================================
# GROUP VARIANTS
# ============================================================

groups = defaultdict(lambda: {
    "VUS": [],
    "P/LP": []
})

for x in data:

    key = (
        x["gene"],
        x["utr"],
        x["transcript"]
    )

    groups[key][x["category"]].append(x)


# ============================================================
# BUILD USABLE GROUPS
# ============================================================

usable_groups = {}

for key, values in groups.items():

    gene, utr, transcript = key
    coord_map = coordinate_maps.get((transcript, utr))

    if not coord_map:
        continue

    vus_coords = [
        coord_map[v["pos"]]
        for v in values["VUS"]
        if v["pos"] in coord_map
    ]

    plp_coords = [
        coord_map[p["pos"]]
        for p in values["P/LP"]
        if p["pos"] in coord_map
    ]

    if not vus_coords or not plp_coords:
        continue

    usable_vus = []

    for v in vus_coords:

        if any(p != v for p in plp_coords):
            usable_vus.append(v)

    if usable_vus:
        usable_groups[key] = {
            "vus": usable_vus,
            "plp": plp_coords
        }


# ============================================================
# CALCULATE GROUP CONTRIBUTIONS
# ============================================================

def nearest_distance(v, plp):
    return min(
        abs(v - p)
        for p in plp
        if p != v
    )


group_stats = {}

for key, values in usable_groups.items():

    distances = []

    for v in values["vus"]:
        distances.append(
            nearest_distance(v, values["plp"])
        )

    group_stats[key] = {
        "n_vus": len(distances),
        "mean_distance": sum(distances) / len(distances)
    }


# ============================================================
# OVERALL OBSERVED
# ============================================================

all_distances = []

for key, values in usable_groups.items():

    for v in values["vus"]:
        all_distances.append(
            nearest_distance(v, values["plp"])
        )

overall_mean = sum(all_distances) / len(all_distances)


print()
print("=" * 75)
print("LEAVE-ONE-GROUP-OUT ROBUSTNESS")
print("=" * 75)

print()
print("Groups:", len(usable_groups))
print("VUS:", len(all_distances))
print("Overall observed mean nearest distance:",
      round(overall_mean, 2))

# ============================================================
# REMOVE EACH GROUP ONE AT A TIME
# ============================================================

results = []

for removed_key in usable_groups:

    remaining_distances = []

    for key, values in usable_groups.items():

        if key == removed_key:
            continue

        for v in values["vus"]:
            remaining_distances.append(
                nearest_distance(v, values["plp"])
            )

    mean_without = (
        sum(remaining_distances) /
        len(remaining_distances)
    )

    results.append({
        "key": removed_key,
        "mean": mean_without,
        "change": mean_without - overall_mean
    })


# ============================================================
# EXTREMES
# ============================================================

lowest = min(results, key=lambda x: x["mean"])
highest = max(results, key=lambda x: x["mean"])

print()
print("LOWEST mean after removing one group")
print("-" * 75)

print(
    lowest["key"],
    "| mean =",
    round(lowest["mean"], 2),
    "| change =",
    round(lowest["change"], 2)
)

print()
print("HIGHEST mean after removing one group")
print("-" * 75)

print(
    highest["key"],
    "| mean =",
    round(highest["mean"], 2),
    "| change =",
    round(highest["change"], 2)
)


# ============================================================
# TOP 10 MOST INFLUENTIAL
# ============================================================

print()
print("10 MOST INFLUENTIAL GROUPS")
print("-" * 75)

influential = sorted(
    results,
    key=lambda x: abs(x["change"]),
    reverse=True
)[:10]

for r in influential:

    gene, utr, transcript = r["key"]

    print(
        f"{gene:12s} "
        f"{utr:15s} "
        f"VUS={group_stats[r['key']]['n_vus']:4d} "
        f"group_mean={group_stats[r['key']]['mean_distance']:8.2f} "
        f"without={r['mean']:8.2f} "
        f"change={r['change']:8.2f}"
    )


# ============================================================
# RANGE
# ============================================================

print()
print("ROBUSTNESS RANGE")
print("-" * 75)

print(
    "Minimum mean after removing a group:",
    round(lowest["mean"], 2)
)

print(
    "Maximum mean after removing a group:",
    round(highest["mean"], 2)
)

print(
    "Overall mean:",
    round(overall_mean, 2)
)

print()
print("=" * 75)
print("AUDIT COMPLETE")
print("=" * 75)
