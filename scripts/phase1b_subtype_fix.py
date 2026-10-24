# Phase 1b: corrected subtype statistics. Loads the NMF scores already saved by phase1_rerun_cells.py,
# so it runs in seconds. Fixes the bug where empty categorical levels turned every Kruskal-Wallis into NaN.
# EDIT the ATLAS and STUDY_COL lines to match what you used in phase1_rerun_cells.py.
import warnings
warnings.filterwarnings("ignore")
from pathlib import Path
from itertools import combinations
import numpy as np
import pandas as pd
import anndata as ad
from scipy.stats import kruskal, mannwhitneyu, rankdata
from statsmodels.stats.multitest import multipletests
import statsmodels.formula.api as smf
import statsmodels.api as sm

pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 50)
pd.set_option("display.max_colwidth", 45)

OUT = Path("data/outputs")
ATLAS = Path("data/IntegratedAtlas.h5ad")   # <-- same path you used before
STUDY_COL = "batch"                          # <-- same column you used before

W_df = pd.read_csv(OUT / "weighted_patient_program_scores_v2.csv", index_col=0)
prog_cols = list(W_df.columns)

# ---- patient labels, forced to plain strings (this is the fix) ----
adata = ad.read_h5ad(ATLAS, backed="r")
obs = adata.obs[["donor_id", "disease", STUDY_COL]].copy()
for c in obs.columns:
    obs[c] = obs[c].astype(str)

# sanity check: does any donor carry more than one label?
multi = obs.groupby("donor_id").nunique()
print("donors with >1 disease label:", int((multi["disease"] > 1).sum()),
      "| donors with >1 study:", int((multi[STUDY_COL] > 1).sum()))

mode = lambda s: s.mode().iloc[0]
labels = obs.groupby("donor_id").agg(disease=("disease", mode), study=(STUDY_COL, mode))
patient_df = W_df.join(labels, how="inner")

# ---- keep disease groups with >= 5 donors ----
counts = patient_df["disease"].value_counts()
keep = counts[counts >= 5].index
pdf = patient_df[patient_df["disease"].isin(keep)].copy()
k_groups, N = pdf["disease"].nunique(), len(pdf)
print(f"\nN = {N} patients, {k_groups} groups")
print(pdf["disease"].value_counts(), "\n")

print("DISEASE x STUDY (shows how tangled the labels are with study):")
print(pd.crosstab(pdf["disease"], pdf["study"]), "\n")

# ---- Kruskal-Wallis across disease labels ----
rows = []
for p in prog_cols:
    groups = [g[p].values for _, g in pdf.groupby("disease")]
    H, pval = kruskal(*groups)
    rows.append({"program": p, "H": H, "df": k_groups - 1, "p": pval,
                 "eta2": (H - k_groups + 1) / (N - k_groups)})
res = pd.DataFrame(rows)
res["FDR"] = multipletests(res["p"], method="fdr_bh")[1]
print("KRUSKAL-WALLIS BY DISEASE LABEL")
print(res.round(6).to_string(index=False), "\n")
res.to_csv(OUT / "program_subtype_association_v2.csv", index=False)

print("MEDIAN PROGRAM LOADING BY GROUP")
print(pdf.groupby("disease")[prog_cols].median().round(3).T, "\n")

# ---- same test, but study is also a factor (rank-based ANOVA) ----
adj = []
for p in prog_cols:
    d = pdf[[p, "disease", "study"]].rename(columns={p: "val"}).copy()
    d["y"] = rankdata(d["val"])
    try:
        tab = sm.stats.anova_lm(smf.ols("y ~ C(disease) + C(study)", data=d).fit(), typ=2)
        adj.append({"program": p, "p_disease_adj_study": tab.loc["C(disease)", "PR(>F)"],
                    "p_study": tab.loc["C(study)", "PR(>F)"]})
    except Exception:
        adj.append({"program": p, "p_disease_adj_study": np.nan, "p_study": np.nan})
adj = pd.DataFrame(adj)
adj["FDR_disease_adj"] = multipletests(adj["p_disease_adj_study"].fillna(1), method="fdr_bh")[1]
print("DISEASE EFFECT ADJUSTED FOR STUDY  |  STUDY EFFECT")
print(adj.round(6).to_string(index=False), "\n")
adj.to_csv(OUT / "program_subtype_adjusted_for_study_v2.csv", index=False)

# ---- strongest control for batch: compare disease groups WITHIN a single study ----
wrows = []
for study, sdf in pdf.groupby("study"):
    gs = [g for g, n in sdf["disease"].value_counts().items() if n >= 5]
    for a, b in combinations(sorted(gs), 2):
        for p in prog_cols:
            x = sdf.loc[sdf["disease"] == a, p]
            y = sdf.loc[sdf["disease"] == b, p]
            U, pv = mannwhitneyu(x, y, alternative="two-sided")
            wrows.append({"study": study, "group_a": a, "group_b": b, "n_a": len(x), "n_b": len(y),
                          "program": p, "p": pv,
                          "rank_biserial_a_vs_b": 2 * U / (len(x) * len(y)) - 1})
if wrows:
    within = pd.DataFrame(wrows)
    within["FDR"] = multipletests(within["p"], method="fdr_bh")[1]
    within = within.sort_values("p")
    print("WITHIN-STUDY COMPARISONS (disease groups with >=5 donors in the same study)")
    print(within.head(15).round(4).to_string(index=False))
    within.to_csv(OUT / "within_study_comparisons_v2.csv", index=False)
else:
    print("No study contains two disease groups with >= 5 donors, so no within-study test is possible.")
