# Phase 1c: TreeSHAP on the surrogate models, restricted to programs that passed BOTH criteria:
#   (1) factor reproducible across random inits / subsamples, (2) surrogate CV R2 consistently high.
# Requires: pip install shap
import json, warnings
warnings.filterwarnings("ignore")
from pathlib import Path
import numpy as np
import pandas as pd
import shap
from scipy.stats import spearmanr
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.model_selection import KFold

pd.set_option("display.width", 250)
pd.set_option("display.max_colwidth", 60)

OUT = Path("data/outputs")
PROGRAMS = ["Program_1", "Program_3", "Program_6", "Program_9"]
TOPN = 20
SEED = 42

weighted = pd.read_csv(OUT / "weighted_patient_communication_matrix.csv", index_col=0)
retained = json.load(open(OUT / "retained_features_v2.json"))
X = weighted[retained]
W_df = pd.read_csv(OUT / "weighted_patient_program_scores_v2.csv", index_col=0).loc[X.index]
H_df = pd.read_csv(OUT / "weighted_program_feature_weights_v2.csv", index_col=0)
assert list(H_df.columns) == retained, "H columns do not match retained features"

def mean_abs_shap(model, data):
    sv = shap.TreeExplainer(model).shap_values(data)
    return np.abs(sv).mean(axis=0)

def make_model():
    return ExtraTreesRegressor(n_estimators=200, random_state=SEED, n_jobs=-1)

top_rows, summary = [], []
kf = KFold(n_splits=5, shuffle=True, random_state=SEED)

for prog in PROGRAMS:
    y = W_df[prog].values
    full = make_model().fit(X.values, y)
    imp = pd.Series(mean_abs_shap(full, X.values), index=X.columns)
    top = imp.sort_values(ascending=False).head(TOPN)
    nmf_w = H_df.loc[prog]
    nmf_top = set(nmf_w.sort_values(ascending=False).head(TOPN).index)

    # stability of the top-20 list across 5 training folds (Jaccard vs full-data list)
    jacc = []
    for tr, _ in kf.split(X.values):
        m = make_model().fit(X.values[tr], y[tr])
        fold_top = set(pd.Series(mean_abs_shap(m, X.values[tr]), index=X.columns)
                       .sort_values(ascending=False).head(TOPN).index)
        jacc.append(len(fold_top & set(top.index)) / len(fold_top | set(top.index)))

    rho = spearmanr(imp.values, nmf_w.values).correlation
    summary.append({"program": prog,
                    "top20_overlap_with_NMF_top20": len(set(top.index) & nmf_top),
                    "spearman_shap_vs_NMF_weight": round(rho, 3),
                    "top20_fold_jaccard_mean": round(float(np.mean(jacc)), 3),
                    "top20_fold_jaccard_min": round(float(np.min(jacc)), 3)})
    for feat, val in top.items():
        top_rows.append({"program": prog, "feature": feat, "mean_abs_shap": val,
                         "nmf_weight": nmf_w[feat], "in_nmf_top20": feat in nmf_top})
    print(f"\n=== {prog}: top {TOPN} features by mean |SHAP| ===")
    print(pd.DataFrame({"mean_abs_shap": top.round(5),
                        "nmf_weight": nmf_w[top.index].round(4)}))

summary = pd.DataFrame(summary)
print("\nSUMMARY (SHAP explains the surrogate model, not the biology)")
print(summary.to_string(index=False))
pd.DataFrame(top_rows).to_csv(OUT / "shap_top_features_v2.csv", index=False)
summary.to_csv(OUT / "shap_summary_v2.csv", index=False)
