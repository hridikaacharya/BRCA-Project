"""
BRCA exploratory figure generator
Project root: C:\\BRCA\\BRCA Project

Reads existing project outputs and creates figures/reports under:
    figures/generated/

Does not modify input data. Disease association tables are labelled exploratory
until independently reproduced from patient-level metadata and program scores.
"""

from __future__ import annotations

import json
import hashlib
import warnings
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

ROOT = Path(r"C:\BRCA\BRCA Project")
OUT = ROOT / "data" / "outputs"
FIGROOT = ROOT / "figures" / "generated"
PANELS = FIGROOT / "individual"
COMBINED = FIGROOT / "combined"
PNG = FIGROOT / "png_600dpi"
SVG = FIGROOT / "svg"
CAPTIONS = FIGROOT / "captions"

for folder in (PANELS, COMBINED, PNG, SVG, CAPTIONS):
    folder.mkdir(parents=True, exist_ok=True)

# A restrained, consistent palette; plots remain readable when printed.
INK = "#263238"
MUTED = "#607D8B"
GRID = "#D9E0E4"
PALETTE = [
    "#456990", "#5B8E7D", "#E0A458", "#C06C84", "#7B6DAB",
    "#8C9A7A", "#4D9DE0", "#D17A22", "#6C757D", "#B56576",
]
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "axes.titlesize": 11,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "text.color": INK,
    "axes.labelcolor": INK,
    "axes.edgecolor": GRID,
    "xtick.color": INK,
    "ytick.color": INK,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
})

REPORT = []
CAPTION_TEXT = {}

def log(msg=""):
    print(msg)
    REPORT.append(str(msg))

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

def locate(filename: str, folders=None):
    folders = folders or [OUT]
    for folder in folders:
        p = folder / filename
        if p.exists():
            return p
    return None

def read_csv_required(filename: str, folders=None):
    path = locate(filename, folders)
    if path is None:
        log(f"SKIP: required file not found: {filename}")
        return None
    try:
        df = pd.read_csv(path)
        log(f"READ: {path.relative_to(ROOT)} | shape={df.shape}")
        return df
    except Exception as e:
        log(f"SKIP: could not read {filename}: {type(e).__name__}: {e}")
        return None

def clean_program_score_file(filename):
    path = locate(filename)
    if path is None:
        return None
    df = pd.read_csv(path)
    # Accommodate files whose index was saved as an unnamed first column.
    if "patient_id" not in df.columns:
        first = df.columns[0]
        if str(first).startswith("Unnamed"):
            df = df.rename(columns={first: "patient_id"})
        else:
            df = df.rename(columns={first: "patient_id"})
    df["patient_id"] = df["patient_id"].astype(str)
    program_cols = [c for c in df.columns if str(c).startswith("Program_")]
    if len(program_cols) != 10:
        raise ValueError(f"{filename}: expected 10 program columns, found {len(program_cols)}")
    return df[["patient_id"] + program_cols]

def canonicalize_program_labels(df):
    """Map one-based Program_1..10 to canonical Program_0..9 only if needed."""
    if df is None:
        return None
    cols = [c for c in df.columns if c.startswith("Program_")]
    if set(cols) == {f"Program_{i}" for i in range(1, 11)}:
        df = df.rename(columns={f"Program_{i}": f"Program_{i-1}" for i in range(1, 11)})
    expected = [f"Program_{i}" for i in range(10)]
    if not all(c in df.columns for c in expected):
        raise ValueError("Could not normalize program labels to Program_0..Program_9")
    return df[["patient_id"] + expected]

def save_figure(fig, stem, caption):
    """Save the full multi-panel figure plus cropped individual axes."""
    fig.canvas.draw()
    full_outputs = [
        (COMBINED / f"{stem}.pdf", "pdf"),
        (PNG / f"{stem}.png", "png"),
        (SVG / f"{stem}.svg", "svg"),
    ]
    for path, fmt in full_outputs:
        kwargs = {"dpi": 600, "bbox_inches": "tight"} if fmt == "png" else {"bbox_inches": "tight"}
        fig.savefig(path, format=fmt, **kwargs)

    # Export each real plotting axis separately. Colorbar axes are excluded.
    renderer = fig.canvas.get_renderer()
    letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    panel_number = 0
    for ax in fig.axes:
        if not ax.axison or ax.get_label() == "<colorbar>":
            continue
        try:
            bbox = ax.get_tightbbox(renderer)
            if bbox is None:
                continue
            bbox = bbox.transformed(fig.dpi_scale_trans.inverted()).expanded(1.08, 1.12)
            panel_stem = f"{stem}_panel_{letters[panel_number % len(letters)]}"
            for ext, fmt in [("pdf", "pdf"), ("svg", "svg"), ("png", "png")]:
                path = PANELS / f"{panel_stem}.{ext}"
                kwargs = {"dpi": 600, "bbox_inches": bbox} if fmt == "png" else {"bbox_inches": bbox}
                fig.savefig(path, format=fmt, **kwargs)
            CAPTION_TEXT[panel_stem] = caption
            panel_number += 1
        except Exception as e:
            log(f"WARNING: could not export cropped panel from {stem}: {type(e).__name__}: {e}")

    CAPTION_TEXT[stem] = caption
    plt.close(fig)
    log(f"FIGURE: {stem} -> combined PDF/SVG/600-dpi PNG plus {panel_number} cropped panels")

def finish_composite(fig, stem, caption):
    for path, fmt in [
        (COMBINED / f"{stem}.pdf", "pdf"),
        (COMBINED / f"{stem}.svg", "svg"),
        (COMBINED / f"{stem}.png", "png"),
    ]:
        kwargs = {"dpi": 600, "bbox_inches": "tight"} if fmt == "png" else {"bbox_inches": "tight"}
        fig.savefig(path, format=fmt, **kwargs)
    CAPTION_TEXT[stem] = caption
    plt.close(fig)
    log(f"COMPOSITE: {stem} -> PDF, SVG, 600-dpi PNG")

def style_axes(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color(GRID)
    ax.spines["bottom"].set_color(GRID)
    ax.grid(axis="y", color=GRID, linewidth=0.6, alpha=0.65)
    ax.set_axisbelow(True)

def get_atlas_metadata():
    atlas_path = ROOT / "data" / "IntegratedAtlas.h5ad"
    if not atlas_path.exists():
        log("SKIP cohort metadata: IntegratedAtlas.h5ad not found")
        return None
    try:
        import anndata as ad
        atlas = ad.read_h5ad(atlas_path, backed="r")
        cols = [c for c in ["donor_id", "batch", "disease", "grade", "cell_type"]
                if c in atlas.obs.columns]
        obs = atlas.obs[cols].copy()
        atlas.file.close()
        if "donor_id" not in obs:
            log("SKIP cohort metadata: atlas lacks donor_id")
            return None
        # The atlas is cell-level. Reduce to one row per donor and check consistency.
        donor_meta = obs.drop_duplicates("donor_id").copy()
        n_pairs = obs[["donor_id", "batch"]].drop_duplicates().shape[0] if "batch" in obs else None
        log(f"Atlas cells: {len(obs):,}; unique donors: {obs['donor_id'].nunique()}")
        if n_pairs is not None:
            log(f"Unique donor-batch pairs: {n_pairs}")
            if n_pairs != obs["donor_id"].nunique():
                log("WARNING: at least one donor occurs in multiple batches; inspect before donor-level plots.")
        return obs, donor_meta
    except Exception as e:
        log(f"SKIP cohort metadata: {type(e).__name__}: {e}")
        return None

def make_fig1_cohort():
    got = get_atlas_metadata()
    if got is None:
        return None
    obs, donor_meta = got
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.4), constrained_layout=True)

    if "batch" in obs:
        batch = obs["batch"].astype("string").fillna("Missing")
        axes[0].barh(batch.index.astype(str), batch.values, color=MUTED)
        axes[0].set_title("A. Cells contributed by study")
        axes[0].set_xlabel("Number of cells")
        style_axes(axes[0])
    else:
        axes[0].text(0.5, 0.5, "Study/batch metadata unavailable", ha="center", va="center")
        axes[0].set_axis_off()

    if "disease" in obs:
        # Cell-level counts are explicitly labelled as such, not as patient counts.
        disease = obs["disease"].astype("string").fillna("Missing")
        axes[1].barh(disease.index.astype(str), disease.values, color=PALETTE[1])
        axes[1].set_title("B. Atlas cells by recorded disease label")
        axes[1].set_xlabel("Number of cells (not patients)")
        axes[1].tick_params(axis="y", labelsize=7)
        style_axes(axes[1])
    else:
        axes[1].text(0.5, 0.5, "Disease metadata unavailable", ha="center", va="center")
        axes[1].set_axis_off()

    caption = (
        "Figure 1. Cohort composition in the integrated single-cell atlas. "
        "Panel A shows cell counts by source study/batch; panel B shows cell counts "
        "by the atlas's recorded disease label. Counts are cells, not independent patients. "
        "Disease labels and cell numbers vary across source studies, so study composition "
        "must be considered when interpreting downstream associations."
    )
    save_figure(fig, "Figure1_cohort_composition", caption)
    return fig

def make_fig2_scores_and_features(scores):
    if scores is None:
        log("SKIP Figure 2: program scores unavailable")
        return None
    prog_cols = [f"Program_{i}" for i in range(10)]
    x = scores[prog_cols].apply(pd.to_numeric, errors="coerce")
    if x.isna().any().any():
        log("SKIP Figure 2 score heatmap: nonnumeric/missing score values detected")
        return None

    # Cluster ordering is intentionally not performed: this is a reproducible display
    # in source patient order, not an inferred clustering result.
    fig, axes = plt.subplots(1, 2, figsize=(15, 7), constrained_layout=True,
                             gridspec_kw={"width_ratios": [1.4, 1]})
    im = axes[0].imshow(x.to_numpy(), aspect="auto", interpolation="nearest",
                        cmap="viridis")
    axes[0].set_title("A. Patient-by-program scores")
    axes[0].set_xlabel("NMF program")
    axes[0].set_ylabel("Patients (source file order)")
    axes[0].set_xticks(range(10))
    axes[0].set_xticklabels([f"P{i}" for i in range(10)])
    cb = fig.colorbar(im, ax=axes[0], fraction=0.025, pad=0.02)
    cb.set_label("NMF coefficient (saved scale)")

    weights = read_csv_required("weighted_program_feature_weights_v2.csv")
    if weights is None:
        weights = read_csv_required("weighted_nmf_program_feature_weights.csv")
        if weights is not None and weights.columns[0].startswith("Unnamed"):
            weights = weights.rename(columns={weights.columns[0]: "program"})
        if weights is not None and "program" in weights.columns:
            weights["program"] = weights["program"].astype(str)
            weights["program"] = weights["program"].map(
                lambda s: f"Program_{int(s.split('_')[-1])-1}"
                if s.startswith("Program_") and s.split("_")[-1].isdigit()
                else s
            )
    else:
        if weights.columns[0].startswith("Unnamed"):
            weights = weights.rename(columns={weights.columns[0]: "program"})

    if weights is not None and "program" in weights.columns:
        feature_cols = [c for c in weights.columns if c != "program"]
        # Only plot the 8 highest-weight features per program to keep labels legible.
        records = []
        for _, row in weights.iterrows():
            prog = str(row["program"])
            vals = pd.to_numeric(row[feature_cols], errors="coerce")
            vals = vals.dropna().sort_values(ascending=False).head(8)
            for feat, val in vals.items():
                records.append((prog, str(feat), float(val)))
        if records:
            top = pd.DataFrame(records, columns=["program", "feature", "weight"])
            top["label"] = top["program"] + " | " + top["feature"]
            top = top.sort_values(["program", "weight"], ascending=[True, True])
            axes[1].barh(top["label"], top["weight"], color=PALETTE[0])
            axes[1].set_title("B. Top eight feature weights per program")
            axes[1].set_xlabel("NMF feature weight")
            axes[1].tick_params(axis="y", labelsize=5.5)
            style_axes(axes[1])
        else:
            axes[1].text(0.5, 0.5, "No usable feature weights", ha="center", va="center")
            axes[1].set_axis_off()
    else:
        axes[1].text(0.5, 0.5, "Feature-weight file unavailable", ha="center", va="center")
        axes[1].set_axis_off()

    caption = (
        "Figure 2. Exploratory NMF representation of patient-level abundance-weighted "
        "ligand–receptor expression-product scores. Panel A displays the ten saved NMF "
        "coefficients for 138 patients in the input-file order; it is not a patient-cluster "
        "analysis. Panel B displays the eight largest saved feature weights per program. "
        "These weights describe the fitted factorization and do not by themselves establish "
        "directional cell-to-cell signaling or functional communication."
    )
    save_figure(fig, "Figure2_program_scores_and_features", caption)
    return fig

def make_fig3_stability():
    summary = read_csv_required("weighted_nmf_k_allpairs_summary.csv")
    init = read_csv_required("weighted_nmf_initialization_stability.csv")
    if summary is None:
        log("SKIP Figure 3: rank stability summary missing")
        return None

    # Detect expected columns flexibly, but never silently substitute a different measure.
    cols = {str(c).lower(): c for c in summary.columns}
    def find_col(candidates):
        for c in candidates:
            if c in cols:
                return cols[c]
        return None
    kcol = find_col(["k", "n_components", "rank"])
    errcol = find_col(["mean_relative_error", "relative_reconstruction_error",
                       "mean_relative_reconstruction_error"])
    coscol = find_col(["mean_pairwise_matched_cosine", "mean_matched_cosine",
                       "mean_pairwise_cosine", "mean_cosine_similarity"])
    if kcol is None:
        log(f"SKIP Figure 3: cannot identify rank column in {list(summary.columns)}")
        return None

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    if errcol:
        xx = pd.to_numeric(summary[kcol], errors="coerce")
        yy = pd.to_numeric(summary[errcol], errors="coerce")
        valid = xx.notna() & yy.notna()
        axes[0].plot(xx[valid], yy[valid], marker="o", linewidth=1.8, color=PALETTE[0])
        axes[0].set_title("A. Reconstruction error by NMF rank")
        axes[0].set_xlabel("NMF rank (k)")
        axes[0].set_ylabel(errcol.replace("_", " "))
        style_axes(axes[0])
    else:
        axes[0].text(0.5, 0.5, f"Error column not identified\nAvailable: {list(summary.columns)}",
                     ha="center", va="center", wrap=True)
        axes[0].set_axis_off()

    if coscol:
        xx = pd.to_numeric(summary[kcol], errors="coerce")
        yy = pd.to_numeric(summary[coscol], errors="coerce")
        valid = xx.notna() & yy.notna()
        axes[1].plot(xx[valid], yy[valid], marker="o", linewidth=1.8, color=PALETTE[1])
        axes[1].set_title("B. Mean matched-component similarity")
        axes[1].set_xlabel("NMF rank (k)")
        axes[1].set_ylabel("Mean matched cosine similarity")
        axes[1].set_ylim(0, 1.02)
        style_axes(axes[1])
    else:
        axes[1].text(0.5, 0.5, f"Cosine column not identified\nAvailable: {list(summary.columns)}",
                     ha="center", va="center", wrap=True)
        axes[1].set_axis_off()

    if init is not None:
        log(f"Initialization stability table columns: {list(init.columns)}")
        log("Initialization stability values are included in the audit report; "
            "a plot is only added if the table's schema is interpretable.")
    caption = (
        "Figure 3. Exploratory NMF rank sensitivity. Panel A shows the saved reconstruction-error "
        "summary across tested ranks; panel B shows the saved mean matched-component cosine "
        "similarity where that metric is available. Lower reconstruction error alone does not "
        "identify a uniquely optimal rank, and mean similarity can conceal individual unstable "
        "components. The tested ranks and initialization settings are exploratory."
    )
    save_figure(fig, "Figure3_nmf_rank_stability", caption)
    return fig

def make_fig4_association_summaries():
    unadj = read_csv_required("program_subtype_association_v2.csv")
    adj = read_csv_required("program_subtype_adjusted_for_study_v2.csv")
    within = read_csv_required("within_study_comparisons_v2.csv")
    if unadj is None and adj is None and within is None:
        log("SKIP Figure 4: no saved association summary tables found")
        return None

    fig, axes = plt.subplots(1, 3, figsize=(15, 6), constrained_layout=True)
    panels = [
        (unadj, "FDR", "eta2", "A. Unadjusted association summaries",
         "Saved eta-squared estimate"),
        (adj, "FDR_disease_adj", "p_disease_adj_study",
         "B. Study-adjusted association summaries", "Saved disease p-value"),
        (within, "FDR", "rank_biserial_a_vs_b",
         "C. Within-study effect sizes", "Saved rank-biserial effect size"),
    ]
    for ax, (df, qcol, valcol, title, xlabel) in zip(axes, panels):
        if df is None or qcol not in df.columns or valcol not in df.columns:
            ax.text(0.5, 0.5, "Required columns unavailable", ha="center", va="center")
            ax.set_axis_off()
            continue
        d = df.copy()
        if "program" not in d.columns and "program" in [str(c).lower() for c in d.columns]:
            pass
        if "program" not in d.columns:
            # Within-study table contains multiple rows per program; select the
            # smallest saved FDR per program only for display, not as a new test.
            if "program" in d.columns:
                pass
            else:
                ax.text(0.5, 0.5, "Program labels unavailable", ha="center", va="center")
                ax.set_axis_off()
                continue
        d[qcol] = pd.to_numeric(d[qcol], errors="coerce")
        d[valcol] = pd.to_numeric(d[valcol], errors="coerce")
        d = d.dropna(subset=[qcol, valcol])
        if "program" in d.columns:
            d = d.sort_values(qcol).drop_duplicates("program", keep="first")
            labels = d["program"].astype(str)
        else:
            labels = d.index.astype(str)
        # Plot effect sizes where available; for adjusted p values use -log10(p).
        values = d[valcol].to_numpy(dtype=float)
        if "p-value" in xlabel.lower() or "p-value" in valcol.lower() or valcol == "p_disease_adj_study":
            values = -np.log10(np.clip(values, 1e-300, 1))
            xlab = r"$-\log_{10}$(saved p-value)"
        else:
            xlab = xlabel
        order = np.argsort(values)
        ax.barh(np.asarray(labels)[order], values[order],
                color=[PALETTE[i % len(PALETTE)] for i in range(len(order))])
        ax.set_title(title + "\n(saved results; not independently reproduced)")
        ax.set_xlabel(xlab)
        ax.tick_params(axis="y", labelsize=7)
        style_axes(ax)

    caption = (
        "Figure 4. Saved exploratory disease-association summaries. These panels visualize "
        "previously generated summary tables and are not an independent reanalysis. Panel A "
        "shows the saved unadjusted effect-size estimates; panel B displays saved disease "
        "p-values from the study-adjusted analysis on a negative-log scale; panel C shows "
        "saved within-study rank-biserial effect sizes where available. Confirm the underlying "
        "patient-level metadata, statistical model, multiple-testing correction, and direction "
        "of contrasts before treating any association as a manuscript result."
    )
    save_figure(fig, "Figure4_saved_association_summaries_exploratory", caption)
    return fig

def save_captions():
    md = [
        "# Figure captions",
        "",
        "Generated automatically by `generate_brca_figures.py`.",
        "Saved association summaries remain exploratory until independently reproduced.",
        "",
    ]
    for name, caption in CAPTION_TEXT.items():
        md.extend([f"## {name.replace('_', ' ')}", "", caption, ""])
    (CAPTIONS / "figure_captions.md").write_text("\n".join(md), encoding="utf-8")
    for name, caption in CAPTION_TEXT.items():
        (CAPTIONS / f"{name}.txt").write_text(caption + "\n", encoding="utf-8")

def main():
    log("BRCA FIGURE GENERATION")
    log(f"Started: {datetime.now().isoformat(timespec='seconds')}")
    log(f"Project root: {ROOT}")
    log(f"Output root: {FIGROOT}")
    log("Input data will not be modified.")
    if not ROOT.exists():
        raise FileNotFoundError(f"Project root does not exist: {ROOT}")
    if not OUT.exists():
        raise FileNotFoundError(f"Expected output folder does not exist: {OUT}")

    # Verify canonical scores against both saved representations.
    try:
        canonical = canonicalize_program_labels(
            clean_program_score_file("weighted_nmf_patient_program_scores.csv")
        )
        newer = canonicalize_program_labels(
            clean_program_score_file("weighted_patient_program_scores_v2.csv")
        )
        if canonical is not None and newer is not None:
            if canonical["patient_id"].tolist() != newer["patient_id"].tolist():
                log("WARNING: score files have different patient order.")
            cols = [f"Program_{i}" for i in range(10)]
            a = canonical[cols].to_numpy(dtype=float)
            b = newer[cols].to_numpy(dtype=float)
            log(f"Score file alignment max abs diff: {np.max(np.abs(a-b)):.12g}")
            log(f"Score file alignment allclose: {np.allclose(a, b, rtol=1e-7, atol=1e-10)}")
            if not np.allclose(a, b, rtol=1e-7, atol=1e-10):
                log("WARNING: score files differ; using weighted_nmf_patient_program_scores.csv.")
        scores = canonical
    except Exception as e:
        log(f"WARNING: score verification failed: {type(e).__name__}: {e}")
        scores = None

    # Cohort composition
    make_fig1_cohort()

    # Program scores / feature weights
    make_fig2_scores_and_features(scores)

    # Rank and initialization stability
    make_fig3_stability()

    # Existing association summaries, explicitly flagged exploratory
    make_fig4_association_summaries()

    save_captions()

    # Manifest hashes all files that this script generated.
    manifest = []
    for p in sorted(FIGROOT.rglob("*")):
        if p.is_file() and p.name not in {"file_manifest.json", "audit_report.txt"}:
            manifest.append({
                "path": str(p.relative_to(ROOT)),
                "size_bytes": p.stat().st_size,
                "sha256": sha256(p),
            })
    (FIGROOT / "file_manifest.json").write_text(
        json.dumps({
            "generated_at": datetime.now().isoformat(timespec="seconds"),
            "project_root": str(ROOT),
            "script": "generate_brca_figures.py",
            "files": manifest,
        }, indent=2),
        encoding="utf-8",
    )

    log("")
    log("=== FINAL OUTPUT SUMMARY ===")
    log(f"Individual panels: {PANELS}")
    log(f"Combined figures: {COMBINED}")
    log(f"600-dpi PNGs: {PNG}")
    log(f"SVG files: {SVG}")
    log(f"Captions: {CAPTIONS}")
    log("Association panels use saved exploratory tables and are not independently reproduced.")
    log("No source data were modified.")
    (FIGROOT / "audit_report.txt").write_text("\n".join(REPORT) + "\n", encoding="utf-8")
    print("\nDone. Review figures/generated/audit_report.txt before using any panel in a manuscript.")

if __name__ == "__main__":
    main()
