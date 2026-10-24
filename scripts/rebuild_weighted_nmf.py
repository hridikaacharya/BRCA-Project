
from pathlib import Path
from datetime import datetime, timezone
import json
import warnings

import numpy as np
import pandas as pd
import sklearn
from sklearn.preprocessing import MinMaxScaler
from sklearn.decomposition import NMF

# --------------------------------------------------
# Configuration
# --------------------------------------------------
ROOT = Path(__file__).resolve().parent
OUT = ROOT / "data" / "outputs"
INPUT = OUT / "weighted_patient_communication_matrix.csv"

K = 10
SEED = 42
MAX_ITER = 2000

# Use new filenames for now; do not overwrite ambiguous
# historical results or label this fit "validated".
W_PATH = OUT / "weighted_nmf_patient_program_scores.csv"
H_PATH = OUT / "weighted_nmf_program_feature_weights.csv"
FEATURES_PATH = OUT / "weighted_nmf_retained_features.csv"
REPORT_PATH = OUT / "weighted_nmf_provenance.json"

# --------------------------------------------------
# Load and validate input
# --------------------------------------------------
if not INPUT.exists():
    raise FileNotFoundError(f"Input not found: {INPUT}")

X_df = pd.read_csv(INPUT, index_col=0)

if X_df.empty:
    raise ValueError("Input matrix is empty.")
if not X_df.index.is_unique:
    raise ValueError("Patient IDs are not unique.")
if not X_df.columns.is_unique:
    raise ValueError("Feature names are not unique.")

X_df = X_df.apply(pd.to_numeric, errors="raise")
X = X_df.to_numpy(dtype=np.float64)

if not np.isfinite(X).all():
    raise ValueError("Input contains NaN or infinite values.")
if (X < 0).any():
    raise ValueError("NMF requires non-negative input values.")

# --------------------------------------------------
# Remove constant features and scale each feature
# across patients to [0, 1].
# --------------------------------------------------
variances = X_df.var(axis=0, ddof=0)
keep = variances > 0

X_df = X_df.loc[:, keep]
if X_df.shape[1] == 0:
    raise ValueError("No non-constant features remain.")

scaler = MinMaxScaler()
X_scaled = scaler.fit_transform(X_df.to_numpy(dtype=np.float64))

if not np.isfinite(X_scaled).all() or (X_scaled < 0).any():
    raise ValueError("Preprocessed matrix failed validation.")

# --------------------------------------------------
# Fit canonical weighted NMF
# --------------------------------------------------
model = NMF(
    n_components=K,
    init="nndsvda",
    random_state=SEED,
    max_iter=MAX_ITER,
    solver="cd",
    tol=1e-4,
)

with warnings.catch_warnings(record=True) as caught:
    warnings.simplefilter("always")
    W = model.fit_transform(X_scaled)
    H = model.components_

# --------------------------------------------------
# Check reconstruction
# --------------------------------------------------
X_hat = W @ H
residual_norm = float(np.linalg.norm(X_scaled - X_hat))
input_norm = float(np.linalg.norm(X_scaled))
relative_error = residual_norm / input_norm if input_norm else None

program_names = [f"Program_{i + 1}" for i in range(K)]

W_df = pd.DataFrame(
    W,
    index=X_df.index,
    columns=program_names,
)
W_df.index.name = "patient_id"

H_df = pd.DataFrame(
    H,
    index=program_names,
    columns=X_df.columns,
)
H_df.index.name = "program"

# --------------------------------------------------
# Save outputs and provenance
# --------------------------------------------------
OUT.mkdir(parents=True, exist_ok=True)

W_df.to_csv(W_PATH)
H_df.to_csv(H_PATH)

pd.DataFrame({
    "feature": X_df.columns,
    "variance_before_scaling": variances.loc[X_df.columns].to_numpy(),
}).to_csv(FEATURES_PATH, index=False)

report = {
    "run_utc": datetime.now(timezone.utc).isoformat(),
    "input_file": str(INPUT.relative_to(ROOT)),
    "input_shape": list(X_df.shape),
    "original_input_shape": list(X_df.shape[:1]) + [int(len(keep))],
    "patients": int(X_df.shape[0]),
    "features_before_filtering": int(len(keep)),
    "features_removed_zero_variance": int((~keep).sum()),
    "features_retained": int(X_df.shape[1]),
    "preprocessing": (
        "Removed zero-variance features; MinMax-scaled each feature "
        "across patients to [0, 1]."
    ),
    "model": "sklearn.decomposition.NMF",
    "n_components": K,
    "init": "nndsvda",
    "random_state": SEED,
    "solver": "cd",
    "max_iter": MAX_ITER,
    "tolerance": 1e-4,
    "iterations_used": int(model.n_iter_),
    "converged_without_convergence_warning": not any(
        "Maximum number of iterations" in str(w.message)
        for w in caught
    ),
    "reconstruction_error_frobenius": float(model.reconstruction_err_),
    "relative_reconstruction_error": relative_error,
    "sklearn_version": sklearn.__version__,
    "outputs": {
        "patient_scores": str(W_PATH.relative_to(ROOT)),
        "feature_weights": str(H_PATH.relative_to(ROOT)),
        "retained_features": str(FEATURES_PATH.relative_to(ROOT)),
    },
    "warnings": [str(w.message) for w in caught],
}

REPORT_PATH.write_text(json.dumps(report, indent=2), encoding="utf-8")

# --------------------------------------------------
# Report
# --------------------------------------------------
print("\nWEIGHTED NMF REBUILD COMPLETE")
print("-" * 40)
print("Input:", INPUT)
print("Patients:", W_df.shape[0])
print("Features retained:", H_df.shape[1])
print("Programs:", H_df.shape[0])
print("W shape:", W_df.shape)
print("H shape:", H_df.shape)
print("Iterations:", model.n_iter_)
print("Reconstruction error:", model.reconstruction_err_)
print("Relative reconstruction error:", relative_error)
print("Convergence warning:", not report["converged_without_convergence_warning"])
print("\nSaved:")
print(W_PATH)
print(H_PATH)
print(FEATURES_PATH)
print(REPORT_PATH)

if caught:
    print("\nWarnings:")
    for warning in caught:
        print("-", warning.message)
