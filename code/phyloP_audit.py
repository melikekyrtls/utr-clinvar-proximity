import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu

BASE = r"..\data"
FILE = BASE + r"\canonical_utr_vus_gnomad_annotations_clean.tsv"

df = pd.read_csv(FILE, sep="\t", low_memory=False)
df["phyloP"] = pd.to_numeric(df["phyloP"], errors="coerce")

sub = df.dropna(subset=["phyloP", "Review_Category"]).copy()

sub["Multiple"] = (
    sub["Review_Category"].str.lower() == "multiple submitters"
).astype(int)

single = sub.loc[sub["Multiple"] == 0, "phyloP"].to_numpy()
multiple = sub.loc[sub["Multiple"] == 1, "phyloP"].to_numpy()

print("=" * 70)
print("phyloP AUDIT")
print("=" * 70)

print("Total rows:", len(df))
print("phyloP non-missing:", len(sub))
print("Single:", len(single))
print("Multiple:", len(multiple))

print("\nSINGLE")
print("Mean:", np.mean(single))
print("Median:", np.median(single))

print("\nMULTIPLE")
print("Mean:", np.mean(multiple))
print("Median:", np.median(multiple))

u, p = mannwhitneyu(multiple, single, alternative="two-sided")

print("\nMANN-WHITNEY")
print("U:", u)
print("P:", p)

rng = np.random.default_rng(42)
observed_diff = np.mean(multiple) - np.mean(single)

combined = np.concatenate([multiple, single])
n_multiple = len(multiple)

n_perm = 5000
perm_diffs = np.empty(n_perm)

for i in range(n_perm):
    shuffled = rng.permutation(combined)
    perm_multiple = shuffled[:n_multiple]
    perm_single = shuffled[n_multiple:]
    perm_diffs[i] = np.mean(perm_multiple) - np.mean(perm_single)

p_perm = (
    np.sum(np.abs(perm_diffs) >= abs(observed_diff)) + 1
) / (n_perm + 1)

print("\nPERMUTATION")
print("Observed mean difference:", observed_diff)
print("Permutation mean difference:", np.mean(perm_diffs))
print("Permutation 95% interval:",
      np.percentile(perm_diffs, 2.5),
      np.percentile(perm_diffs, 97.5))
print("Permutation p:", p_perm)

print("\nPHYLOP RANGE")
print("Min:", np.min(sub["phyloP"]))
print("Max:", np.max(sub["phyloP"]))
print("Missing:", df["phyloP"].isna().sum())
