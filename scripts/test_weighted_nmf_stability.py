
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.decomposition import NMF
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics.pairwise import cosine_similarity
from scipy.optimize import linear_sum_assignment

OUT = Path("data/outputs")
X_df = pd.read_csv(OUT / "weighted_matrix_3465.csv", index_col=0)
X = MinMaxScaler().fit_transform(X_df.to_numpy(dtype=float))

def fit_nmf(seed):
    model = NMF(
        n_components=10,
        init="random",
        solver="cd",
        max_iter=2000,
        tol=1e-4,
        random_state=seed
    )
    W = model.fit_transform(X)
    return W, model.components_, model

# Use the existing NNDSVDA weighted solution as the reference.
reference_df = pd.read_csv(
    OUT / "program_feature_weights_weighted.csv",
    index_col=0
)
reference_H = reference_df.to_numpy(dtype=float)

rows = []
for seed in range(1, 21):
    W, H, model = fit_nmf(seed)
    similarity = cosine_similarity(reference_H, H)
    ref_idx, run_idx = linear_sum_assignment(-similarity)

    for i, j in zip(ref_idx, run_idx):
        rows.append({
            "seed": seed,
            "reference_program": i,
            "matched_run_program": j,
            "cosine_similarity": float(similarity[i, j]),
            "relative_reconstruction_error": float(
                np.linalg.norm(X - W @ H, "fro") /
                max(np.linalg.norm(X, "fro"), 1e-12)
            )
        })

detail = pd.DataFrame(rows)
detail.to_csv(OUT / "weighted_nmf_component_stability.csv", index=False)

summary = (
    detail.groupby("reference_program")["cosine_similarity"]
    .agg(["mean", "median", "min", "max", "std"])
    .reset_index()
)
summary.to_csv(OUT / "weighted_nmf_component_stability_summary.csv", index=False)

print("Component-level stability against the NNDSVDA reference:")
print(summary.to_string(index=False))
print("\nLowest similarities by run:")
print(
    detail.sort_values("cosine_similarity")
    .groupby("seed", group_keys=False)
    .head(1)
    .sort_values("cosine_similarity")
    .to_string(index=False)
)
print("\nSaved component-level stability tables.")
