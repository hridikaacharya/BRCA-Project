# BRCA Project

## Overview

This repository contains the analysis pipeline for a study of patient level cell–cell communication in breast cancer.

The project uses an integrated single cell RNA-seq atlas to identify recurrent communication programs, characterize their biological properties, and interpret them using surrogate machine learning models and SHAP. The final stage of the project will evaluate these findings using independent spatial transcriptomic data.

---

## Research Question

Can latent communication programs be identified from patient specific ligand receptor interactions, and can those programs be interpreted in terms of the signaling events that define them?

---

## Dataset

Primary dataset

* Integrated Breast Cancer Atlas
* 138 patients
* around 621,000 cells
* Multiple breast cancer subtypes
* Patient level clinical metadata

Large input datasets are excluded from this repository because of their size.

---

## Analysis Workflow

```
Integrated atlas
    ↓
Cell–cell communication inference
    ↓
Patient × ligand to receptor matrix
    ↓
Network weighting
    ↓
Non negative Matrix Factorization
    ↓
Communication programs
    ↓
Program validation
    ↓
Biological characterization
    ↓
Surrogate machine learning models
    ↓
SHAP interpretation
    ↓
Prioritization of ligand to receptor interactions
    ↓
Spatial validation
```

---

# Repository Structure

```
BRCA Project/

data/
    outputs/
    supplementary/

notebooks/

figures/

README.md
requirements.txt
```

Large datasets (`.h5ad`, spatial data, raw sequencing files) are ignored using `.gitignore`.

---

# Notebook Guide

## Notebooks 1–10

These notebooks prepare the data used throughout the project.

Tasks include:

* loading the integrated atlas
* preprocessing
* communication inference
* construction of patient-level communication matrices
* network weighting
* export of intermediate datasets

Main outputs:

* communication matrix
* weighted communication matrix

---

## Notebook 11: Communication Program Discovery

Purpose

Identify recurrent communication programs shared across patients.

Methods

* feature filtering
* MinMax scaling
* Non negative Matrix Factorization (NMF)

Outputs

* patient program scores
* program feature weights

Additional analyses

* reconstruction error
* stability across random seeds
* cosine similarity

Result

Repeated NMF runs converged to the same solution under different random seeds, indicating stable communication programs.

---

## Notebook 12: Program Validation

Purpose

Evaluate the robustness of the NMF decomposition.

Analyses

* reconstruction error across different values of *k*
* elbow analysis
* stability across repeated runs
* feature overlap
* sparsity of program weights
* PCA of patient program activities
* correlation between communication programs

Outputs

* validated patient program matrix
* validated program feature matrix

---

## Notebook 13: Biological Characterization

Purpose

Determine whether communication programs are associated with disease groups.

Methods

* Kruskal Wallis test
* Benjamini Hochberg FDR correction
* eta-squared effect size

Outputs

* program disease association table

This notebook is intended as biological characterization rather than predictive modelling.

---

## Notebook 14: Surrogate Models

Purpose

Evaluate whether communication programs can be approximated by supervised models and interpreted with SHAP.

Models evaluated

* Random Forest
* Extra Trees
* XGBoost
* LightGBM

Evaluation metrics

* R²
* RMSE
* MAE
* cross validation

Current observations

* XGBoost achieved the highest average performance across communication programs.
* Several programs were reconstructed with high fidelity.
* Other programs remained difficult to model, suggesting differences in complexity or signal strength.

SHAP values were computed using the best-performing model for each communication program.

Outputs

* model comparison
* program specific model performance
* SHAP importance scores

---

## Notebook 15

Current status: in progress

Planned analyses

* combine NMF weights with SHAP importance
* identify ligand to receptor interactions supported by both approaches
* compare programs
* prepare candidate interactions for downstream validation

---

## Planned Work

* spatial transcriptomic validation
* pathway enrichment
* network visualization
* manuscript figures
* manuscript preparation

---

# Software

Python packages used in this project include

* Scanpy
* NumPy
* Pandas
* SciPy
* scikit-learn
* XGBoost
* LightGBM
* SHAP
* Matplotlib

Package versions are listed in `requirements.txt`.

---

# Current Status

Completed

* communication inference
* patient communication matrix
* communication program discovery
* program validation
* biological characterization
* surrogate modelling
* SHAP analysis

In progress

* integrated prioritization of ligand–receptor interactions

Planned

* spatial validation
* manuscript figures
* publication
