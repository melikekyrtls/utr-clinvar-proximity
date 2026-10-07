import pandas as pd
import numpy as np
import statsmodels.formula.api as smf

BASE = r"..\data"

df = pd.read_csv(
    BASE + r"\canonical_utr_vus_gnomad_annotations_clean.tsv",
    sep="\t"
)

mane = pd.read_csv(
    BASE + r"\MANE_UTR_coordinates_v2.tsv",
    sep="\t"
)

# UTR length
lengths = (
    mane.groupby(["Transcript", "UTR_Type"])
    .apply(lambda x: x["Stop"].max() - x["Start"].min() + 1)
    .reset_index(name="UTR_Length")
)

df = df.merge(
    lengths,
    on=["Transcript", "UTR_Type"],
    how="left"
)

# WES AF mevcut olan satırlar
df["WES_AF"] = pd.to_numeric(df["WES_AF"], errors="coerce")
df["Distance"] = pd.to_numeric(df["Distance"], errors="coerce")

sub = df.dropna(
    subset=["WES_AF", "Distance", "UTR_Length", "Review_Category"]
).copy()

# Outcome
sub["Multiple"] = (
    sub["Review_Category"].str.lower() == "multiple submitters"
).astype(int)

# Distance transform
sub["log_distance"] = np.log10(sub["Distance"] + 1)

# AF transform
sub["log_wes_af"] = np.log10(sub["WES_AF"] + 1e-8)

print("=" * 70)
print("WES AF REGRESSION")
print("=" * 70)
print("N:", len(sub))
print("Multiple:", sub["Multiple"].sum())
print("Single:", (sub["Multiple"] == 0).sum())

# Model A: same model, no AF
model_a = smf.logit(
    "Multiple ~ log_distance + UTR_Length + C(UTR_Type)",
    data=sub
).fit(disp=False)

# Model B: add WES AF
model_b = smf.logit(
    "Multiple ~ log_distance + UTR_Length + C(UTR_Type) + log_wes_af",
    data=sub
).fit(disp=False)

def report(model, name):
    print("\n" + "-" * 70)
    print(name)
    print("-" * 70)

    for term in ["log_distance", "UTR_Length", "C(UTR_Type)[T.5prime_UTR]", "log_wes_af"]:
        if term in model.params.index:
            beta = model.params[term]
            orr = np.exp(beta)
            ci = model.conf_int().loc[term]
            print(
                f"{term}: "
                f"OR={orr:.6f}, "
                f"95% CI=({np.exp(ci[0]):.6f}, {np.exp(ci[1]):.6f}), "
                f"p={model.pvalues[term]:.6g}"
            )

report(model_a, "MODEL A: WITHOUT WES AF")
report(model_b, "MODEL B: WITH WES AF")

print("\n" + "=" * 70)
print("DISTANCE COMPARISON")
print("=" * 70)

print("Distance OR without AF:", np.exp(model_a.params["log_distance"]))
print("Distance OR with AF:   ", np.exp(model_b.params["log_distance"]))

change = (
    (np.exp(model_b.params["log_distance"]) -
     np.exp(model_a.params["log_distance"]))
    / np.exp(model_a.params["log_distance"])
    * 100
)

print("OR change (%):", change)
