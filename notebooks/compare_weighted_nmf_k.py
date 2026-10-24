
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.decomposition import NMF
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics.pairwise import cosine_similarity
from scipy.optimize import linear_sum_assignment

OUT = Path("data/outputs")
df = pd.read_csv(OUT / "weighted_matrix_3465.csv", index_col=0)
X = MinMaxScaler().fit_transform(df.to_numpy(dtype=float))
norm_x = max(np.linalg.norm(X, ord="fro"), 1e-12)

results = []

for k in range(5, 13):
    fitted = []

    for seed in range(5):
        model = NMF(
            n_components=k,
            init="random",
            solver="cd",
            max_iter=2000,
            tol=1e-4,
            random_state=seed
        )
        W = model.fit_transform(X)
        H = model.components_
        rel_error = np.linalg.norm(X - W @ H, ord="fro") / norm_x
        fitted.append((H, rel_error, model.n_iter_))

    # Compare each run's components with run 0 after optimal matching.
    reference_H = fitted[0][0]
    pairwise_means = []
    pairwise_mins = []

    for H, _, _ in fitted[1:]:
        sim = cosine_similarity(reference_H, H)
        r, c = linear_sum_assignment(-sim)
        matched = sim[r, c]
        pairwise_means.append(float(matched.mean()))
        pairwise_mins.append(float(matched.min()))

    errors = [x[1] for x in fitted]
    iterations = [x[2] for x in fitted]

    results.append({
        "k": k,
        "best_relative_error": min(errors),
        "mean_relative_error": np.mean(errors),
        "worst_relative_error": max(errors),
        "mean_matched_cosine_vs_seed0": np.mean(pairwise_means),
        "worst_matched_cosine_vs_seed0": min(pairwise_mins),
        "runs_reaching_iteration_limit": sum(i >= 2000 for i in iterations)
    })

summary = pd.DataFrame(results)
summary.to_csv(OUT / "weighted_nmf_k_comparison.csv", index=False)
print(summary.to_string(index=False))
print("\nSaved weighted_nmf_k_comparison.csv")
