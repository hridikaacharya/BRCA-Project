
from pathlib import Path
from itertools import combinations

import numpy as np
import pandas as pd
from sklearn.decomposition import NMF
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics.pairwise import cosine_similarity
from scipy.optimize import linear_sum_assignment

OUT = Path("data/outputs")
INPUT = OUT / "weighted_matrix_3465.csv"

N_SEEDS = 10
K_VALUES = range(5, 13)

df = pd.read_csv(INPUT, index_col=0)
X = MinMaxScaler().fit_transform(df.to_numpy(dtype=float))
norm_x = max(np.linalg.norm(X, ord="fro"), 1e-12)

summary_rows = []
pair_rows = []
component_rows = []

for k in K_VALUES:
    print(f"\nTesting k={k}...", flush=True)
    fitted = []

    for seed in range(N_SEEDS):
        model = NMF(
            n_components=k,
            init="random",
            solver="cd",
            max_iter=2000,
            tol=1e-4,
            random_state=seed,
        )

        W = model.fit_transform(X)
        H = model.components_

        rel_error = (
            np.linalg.norm(X - W @ H, ord="fro") / norm_x
        )

        fitted.append({
            "seed": seed,
            "H": H,
            "error": float(rel_error),
            "iterations": int(model.n_iter_),
            "converged": model.n_iter_ < 2000,
        })

    all_pair_means = []
    all_pair_mins = []
    all_pair_p10 = []

    # Compare every unique pair of runs.
    for run_a, run_b in combinations(fitted, 2):
        sim = cosine_similarity(run_a["H"], run_b["H"])
        rows, cols = linear_sum_assignment(-sim)
        matched = sim[rows, cols]

        all_pair_means.append(float(matched.mean()))
        all_pair_mins.append(float(matched.min()))
        all_pair_p10.append(float(np.quantile(matched, 0.10)))

        pair_rows.append({
            "k": k,
            "seed_a": run_a["seed"],
            "seed_b": run_b["seed"],
            "mean_matched_cosine": float(matched.mean()),
            "median_matched_cosine": float(np.median(matched)),
            "p10_matched_cosine": float(np.quantile(matched, 0.10)),
            "minimum_matched_cosine": float(matched.min()),
            "relative_error_a": run_a["error"],
            "relative_error_b": run_b["error"],
        })

        for a, b in zip(rows, cols):
            component_rows.append({
                "k": k,
                "seed_a": run_a["seed"],
                "seed_b": run_b["seed"],
                "component_a": int(a),
                "component_b": int(b),
                "cosine_similarity": float(sim[a, b]),
            })

    errors = [run["error"] for run in fitted]
    iterations = [run["iterations"] for run in fitted]
    n_pairs = len(all_pair_means)

    summary_rows.append({
        "k": k,
        "n_seeds": N_SEEDS,
        "n_pairs": n_pairs,
        "best_relative_error": min(errors),
        "mean_relative_error": float(np.mean(errors)),
        "worst_relative_error": max(errors),
        "mean_pairwise_matched_cosine": float(np.mean(all_pair_means)),
        "median_pairwise_matched_cosine": float(np.median(all_pair_means)),
        "p10_pairwise_mean_cosine": float(np.quantile(all_pair_means, 0.10)),
        "mean_pairwise_worst_component_cosine": float(np.mean(all_pair_mins)),
        "worst_component_cosine_across_pairs": float(min(all_pair_mins)),
        "mean_pairwise_p10_component_cosine": float(np.mean(all_pair_p10)),
        "runs_reaching_iteration_limit": sum(
            i >= 2000 for i in iterations
        ),
    })

summary = pd.DataFrame(summary_rows)
pairs = pd.DataFrame(pair_rows)
components = pd.DataFrame(component_rows)

summary.to_csv(OUT / "weighted_nmf_k_allpairs_summary.csv", index=False)
pairs.to_csv(OUT / "weighted_nmf_k_allpairs_pairs.csv", index=False)
components.to_csv(
    OUT / "weighted_nmf_k_allpairs_components.csv", index=False
)

print("\n=== ALL-PAIRS RANK SUMMARY ===")
print(summary.round(4).to_string(index=False))
print("\nSaved:")
print(OUT / "weighted_nmf_k_allpairs_summary.csv")
print(OUT / "weighted_nmf_k_allpairs_pairs.csv")
print(OUT / "weighted_nmf_k_allpairs_components.csv")
