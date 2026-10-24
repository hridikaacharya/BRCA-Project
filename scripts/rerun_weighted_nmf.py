
from pathlib import Path
import json
import numpy as np
import pandas as pd
from sklearn.decomposition import NMF
from sklearn.preprocessing import MinMaxScaler

# ---------- Paths ----------
BASE = Path(".")
INPUT = BASE / "data" / "outputs" / "weighted_patient_communication_matrix.csv"
OUT = BASE / "data" / "outputs"
OUT.mkdir(parents=True, exist_ok=True)

# ---------- Load and validate ----------
if not INPUT.exists():
    raise FileNotFoundError(f"Missing weighted matrix: {INPUT.resolve()}")

df = pd.read_csv(INPUT, index_col=0)

print("Input shape:", df.shape)
print("First patient IDs:", df.index[:5].tolist())

if df.shape[0] != 138:
    raise ValueError(f"Expected 138 patients, got {df.shape[0]}")

if df.index.has_duplicates:
    raise ValueError("Duplicate patient IDs found.")

if df.columns.has_duplicates:
    raise ValueError("Duplicate feature names found.")

X_df = df.apply(pd.to_numeric, errors="raise")
X = X_df.to_numpy(dtype=float)

if not np.isfinite(X).all():
    raise ValueError("Matrix contains NaN or infinite values.")

if (X < 0).any():
    raise ValueError(
        "NMF requires nonnegative inputs. Inspect the matrix; do not "
        "silently clip negative values."
    )

# Remove features that have zero variance across patients.
variances = X_df.var(axis=0, ddof=0)
keep = variances > 0
removed = X_df.columns[~keep].tolist()
X_df = X_df.loc[:, keep]

print("Features before filtering:", len(variances))
print("Zero-variance features removed:", len(removed))
print("Features retained:", X_df.shape[1])

if X_df.shape[1] == 0:
    raise ValueError("No nonconstant features remain.")

# Save the exact retained feature list for downstream alignment.
X_df.to_csv(OUT / "weighted_matrix_3465.csv")
pd.Series(removed, name="feature").to_csv(
    OUT / "weighted_nmf_removed_features.csv", index=False
)
pd.Series(X_df.columns, name="feature").to_csv(
    OUT / "weighted_nmf_feature_order.csv", index=False
)

# Match the previous pipeline's scaling, if it used MinMaxScaler.
# This makes the input nonnegative and scales each feature to [0, 1].
scaler = MinMaxScaler()
X_scaled = scaler.fit_transform(X_df.to_numpy(dtype=float))

# ---------- NMF ----------
# NNDSVDa gives a deterministic initialization, so do not call its
# repeated identical runs a stability analysis.
model = NMF(
    n_components=10,
    init="nndsvda",
    solver="cd",
    max_iter=2000,
    tol=1e-4,
    random_state=42
)

W = model.fit_transform(X_scaled)
H = model.components_

patient_scores = pd.DataFrame(
    W,
    index=X_df.index,
    columns=[f"Program_{i}" for i in range(10)]
)
feature_weights = pd.DataFrame(
    H,
    index=[f"Program_{i}" for i in range(10)],
    columns=X_df.columns
)

patient_scores.index.name = "patient_id"
feature_weights.index.name = "program"

patient_scores.to_csv(OUT / "patient_program_scores_weighted.csv")
feature_weights.to_csv(OUT / "program_feature_weights_weighted.csv")

# Record reconstruction metrics and settings.
reconstruction = W @ H
fro_error = float(np.linalg.norm(X_scaled - reconstruction, ord="fro"))
relative_error = float(
    fro_error / max(np.linalg.norm(X_scaled, ord="fro"), 1e-12)
)

metadata = {
    "input_file": str(INPUT),
    "n_patients": int(X_df.shape[0]),
    "n_features_retained": int(X_df.shape[1]),
    "n_features_removed": int(len(removed)),
    "n_components": 10,
    "init": "nndsvda",
    "solver": "cd",
    "max_iter": 2000,
    "tol": 1e-4,
    "random_state": 42,
    "scaling": "MinMaxScaler fitted across all patients",
    "frobenius_reconstruction_error": fro_error,
    "relative_reconstruction_error": relative_error,
    "n_iter": int(model.n_iter_),
    "converged_before_max_iter": bool(model.n_iter_ < 2000)
}

with open(OUT / "weighted_nmf_metadata.json", "w", encoding="utf-8") as f:
    json.dump(metadata, f, indent=2)

print("\nNMF completed.")
print(json.dumps(metadata, indent=2))
print("Saved weighted patient scores and feature weights.")
