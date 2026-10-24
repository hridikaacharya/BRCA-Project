# Phase 1 rerun: weighted NMF, stability, subtype stats, surrogate models
# Paste each "# %%" block into its own notebook cell (or open this file in VS Code / Jupyter as-is).
# Run from the same folder your earlier notebooks ran from (so "data/outputs/..." resolves).
# All new files get a "_v2" suffix so nothing from the earlier run is overwritten.

# %% [Cell 0] Imports and paths  -- EDIT ATLAS and STUDY_COL
import json, warnings
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.decomposition import NMF
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.model_selection import train_test_split, KFold
from sklearn.metrics import r2_score
from sklearn.ensemble import RandomForestRegressor, ExtraTreesRegressor
from scipy.optimize import linear_sum_assignment
from scipy.stats import kruskal, rankdata
from statsmodels.stats.multitest import multipletests
import statsmodels.formula.api as smf
import statsmodels.api as sm

warnings.filterwarnings("ignore", category=UserWarning)

OUT = Path("data/outputs")
WEIGHTED_CSV = OUT / "weighted_patient_communication_matrix.csv"
ATLAS = Path("data/IntegratedAtlas.h5ad")   # <-- EDIT: wherever IntegratedAtlas.h5ad lives
STUDY_COL = "batch"                          # <-- CHECK: the adata.obs column that identifies the study
K = 10
SEED = 42

# %% [Cell 1] Load WEIGHTED matrix, drop zero-variance features, save the retained list
weighted = pd.read_csv(WEIGHTED_CSV, index_col=0)
print("weighted matrix:", weighted.shape)            # expect (138, 3489)

feature_var = weighted.var(axis=0)
retained = feature_var[feature_var > 0].index.tolist()
X = weighted[retained].copy()
print("after zero-variance filter:", X.shape)        # may differ from 3465 (different matrix)

with open(OUT / "retained_features_v2.json", "w") as f:
    json.dump(retained, f)

scaler = MinMaxScaler()
X_scaled = pd.DataFrame(scaler.fit_transform(X), index=X.index, columns=X.columns)

# %% [Cell 2] NMF on the weighted matrix; save W and H
nmf = NMF(n_components=K, init="nndsvda", random_state=SEED, max_iter=2000)
W = nmf.fit_transform(X_scaled.values)
H = nmf.components_
print("W", W.shape, "H", H.shape, "reconstruction error", round(nmf.reconstruction_err_, 4),
      "iterations", nmf.n_iter_)

prog_cols = [f"Program_{i}" for i in range(K)]
W_df = pd.DataFrame(W, index=X_scaled.index, columns=prog_cols)
H_df = pd.DataFrame(H, index=prog_cols, columns=X_scaled.columns)
W_df.to_csv(OUT / "weighted_patient_program_scores_v2.csv")
H_df.to_csv(OUT / "weighted_program_feature_weights_v2.csv")

# %% [Cell 3] REAL stability test: random init x 50 seeds, plus 80% patient subsampling
def match_to_reference(H_ref, H_new):
    """Hungarian-match factors by cosine similarity; return matched cosines in reference order."""
    sim = cosine_similarity(H_ref, H_new)
    r, c = linear_sum_assignment(-sim)
    return sim[r, c]

N_RUNS = 50
rng = np.random.default_rng(SEED)

seed_sims, boot_sims = [], []
for s in range(N_RUNS):
    m = NMF(n_components=K, init="random", random_state=1000 + s, max_iter=2000)
    m.fit(X_scaled.values)
    seed_sims.append(match_to_reference(H, m.components_))

    idx = rng.choice(len(X_scaled), size=int(0.8 * len(X_scaled)), replace=False)
    mb = NMF(n_components=K, init="random", random_state=2000 + s, max_iter=2000)
    mb.fit(X_scaled.values[idx])
    boot_sims.append(match_to_reference(H, mb.components_))

stab = pd.DataFrame({
    "random_init_mean_cos": np.mean(seed_sims, axis=0),
    "random_init_min_cos": np.min(seed_sims, axis=0),
    "subsample80_mean_cos": np.mean(boot_sims, axis=0),
    "subsample80_min_cos": np.min(boot_sims, axis=0),
}, index=prog_cols).round(3)
print(stab)
stab.to_csv(OUT / "nmf_stability_v2.csv")
# Report these numbers as-is. Programs with low matched cosine are not reproducible factors.

# %% [Cell 4] Patient labels from the atlas (obs only; does not load expression into memory)
import anndata as ad
adata = ad.read_h5ad(ATLAS, backed="r")
obs = adata.obs[["donor_id", "disease", STUDY_COL]].copy()
mode = lambda s: s.astype(str).mode().iloc[0]
labels = obs.groupby("donor_id", observed=True).agg(disease=("disease", mode),
                                                    study=(STUDY_COL, mode))
missing = set(W_df.index) - set(labels.index)
print("patients in matrix without labels:", len(missing))   # should be 0; if not, fix index format
patient_df = W_df.join(labels, how="inner")
print(patient_df.shape)
print(pd.crosstab(patient_df["disease"], patient_df["study"]))   # shows disease/study confounding

# %% [Cell 5] Kruskal-Wallis by disease label: print EVERY program, with H, df, p, FDR, eta^2
counts = patient_df["disease"].value_counts()
valid = counts[counts >= 5].index
pdf = patient_df[patient_df["disease"].isin(valid)].copy()
k_groups, N = pdf["disease"].nunique(), len(pdf)
print(f"N = {N} patients, {k_groups} groups:\n{pdf['disease'].value_counts()}\n")

rows = []
for p in prog_cols:
    groups = [g[p].values for _, g in pdf.groupby("disease")]
    H_stat, pval = kruskal(*groups)
    eta2 = (H_stat - k_groups + 1) / (N - k_groups)
    rows.append({"program": p, "H": H_stat, "df": k_groups - 1, "p": pval, "eta2": eta2})
res = pd.DataFrame(rows)
res["FDR"] = multipletests(res["p"], method="fdr_bh")[1]
print(res.round(6).to_string(index=False))
res.to_csv(OUT / "program_subtype_association_v2.csv", index=False)

# %% [Cell 6] Sensitivity: disease effect adjusted for study (rank-based ANOVA)
# Caution: disease and study are partly confounded, so some models may be rank-deficient.
adj_rows = []
for p in prog_cols:
    d = pdf[[p, "disease", "study"]].rename(columns={p: "val"}).copy()
    d["y"] = rankdata(d["val"])
    try:
        fit = smf.ols("y ~ C(disease) + C(study)", data=d).fit()
        tab = sm.stats.anova_lm(fit, typ=2)
        adj_rows.append({"program": p,
                         "p_disease_adj_study": tab.loc["C(disease)", "PR(>F)"],
                         "p_study": tab.loc["C(study)", "PR(>F)"]})
    except Exception as e:
        adj_rows.append({"program": p, "p_disease_adj_study": np.nan, "p_study": np.nan})
adj = pd.DataFrame(adj_rows)
adj["FDR_disease_adj"] = multipletests(adj["p_disease_adj_study"].fillna(1), method="fdr_bh")[1]
print(adj.round(6).to_string(index=False))
adj.to_csv(OUT / "program_subtype_adjusted_for_study_v2.csv", index=False)

# %% [Cell 7] Surrogate models on the SAME 3,4xx retained features and the NEW weighted targets
# Holdout (same 80/20 seed as before) AND 5-fold CV. Note: NMF targets were fit on all patients,
# so this measures surrogate fidelity, not independent predictive validity. Say so in the paper.
from xgboost import XGBRegressor

Xs = X.values            # unscaled weighted features, as in the original Notebook 14
Y = W_df.values

def make_models():
    return {
        "RandomForest": lambda: RandomForestRegressor(n_estimators=100, random_state=SEED, n_jobs=-1),
        "ExtraTrees": lambda: ExtraTreesRegressor(n_estimators=100, random_state=SEED, n_jobs=-1),
        "XGBoost": lambda: XGBRegressor(n_estimators=200, learning_rate=0.05, max_depth=4,
                                        random_state=SEED, n_jobs=-1),
    }

# --- holdout
Xtr, Xte, Ytr, Yte = train_test_split(Xs, Y, test_size=0.20, random_state=SEED)
hold = {}
for name, mk in make_models().items():
    hold[name] = [r2_score(Yte[:, j], mk().fit(Xtr, Ytr[:, j]).predict(Xte)) for j in range(K)]
hold_df = pd.DataFrame(hold, index=prog_cols).round(3)
print("HOLDOUT R2 (n_test =", len(Xte), ")\n", hold_df)

# --- 5-fold CV (mean and SD of per-fold R2)
kf = KFold(n_splits=5, shuffle=True, random_state=SEED)
cv_mean, cv_sd = {}, {}
for name, mk in make_models().items():
    scores = np.zeros((5, K))
    for f, (tr, te) in enumerate(kf.split(Xs)):
        for j in range(K):
            scores[f, j] = r2_score(Y[te, j], mk().fit(Xs[tr], Y[tr, j]).predict(Xs[te]))
    cv_mean[name], cv_sd[name] = scores.mean(0), scores.std(0)
cv_df = pd.concat({"mean": pd.DataFrame(cv_mean, index=prog_cols),
                   "sd": pd.DataFrame(cv_sd, index=prog_cols)}, axis=1).round(3)
print("5-FOLD CV R2\n", cv_df)

hold_df.to_csv(OUT / "surrogate_holdout_r2_v2.csv")
cv_df.to_csv(OUT / "surrogate_cv_r2_v2.csv")
# Decide which programs count as "high fidelity" from CV mean AND sd, not from the holdout alone.

# %% [Cell 8] Software versions for the Methods section
import sys, sklearn, xgboost, scipy, statsmodels
print("python", sys.version.split()[0])
for mod in ["numpy", "pandas", "scipy", "sklearn", "xgboost", "statsmodels", "scanpy", "anndata", "liana", "shap"]:
    try:
        m = __import__(mod)
        print(mod, getattr(m, "__version__", "?"))
    except Exception:
        print(mod, "NOT INSTALLED in this environment")
