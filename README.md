# Interpretable Single-Cell Analysis of Breast Cancer Cell–Cell Signaling

## Overview

In this project, I investigated cell-cell signaling patterns in breast cancer using single-cell RNA sequencing data. Analysis combined ligand–receptor features, patient level analysis, non-negative matrix factorization (NMF), and explainable machine learning to examine patterns of signaling activity and their relationship with tumor characteristics.

The analysis identified signaling programs that occur consistently across patients and investigated which features contribute to those programs and their computational prediction.

## Research questions

What recurring programs of ligand–receptor expression features can be identified across patients with breast cancer?

How do the activity scores of the identified programs vary across available tumor characteristics and clinical categories, and how robust are these associations across the eight studies contributing to the atlas?

How accurately can interpretable surrogate models(XGBoost, LightGBM, Random Foresst, Extra Trees) reproduce NMF derived program scores in patients, and which ligand receptor features contribute most to their predictions, assessed using SHAP?

## Data

I used Integrated Breast Cancer Atlas by Chen et al.(2026) This is a published single-cell RNA-sequencing atlas of human breast cancer. The atlas contains 621,200 cells from 138 patients with untreated primary tumors, assembled from eight studies.

Primary publication: Chen et al. (2026), A highly resolved integrated single-cell atlas of human breast cancers. 

Dataset access: CELLxGENE collection
Atlas publication: https://doi.org/10.1093/nargab/lqaf217
Atlas authors' code: https://github.com/montilab/chen_et_al_brca_atlas

The primary input for my study is IntegratedAtlas.h5ad which contains the integrated single-cell expression matrix and the cell level metadata associated with it. I used supplemental files to connect donor identifiers with available clinical and tumour annotations.

I derived pateint-level ligand receptor expression features from the single cell data and used these features for a non-negative matrix factorization (NMF), program characterization, and exploratory analyses of associations with available tumor characteristics.

The original dataset source file is not included in not available in this repository. Users must obtain the source data separately and follow the directory structure to reproduce the results. The tables and figures in this repository are the computational outputs of the project which were not a part of the original atlas. 

## Methods

The analysis pipeline includes:

1. **Dataset exploration and metadata integration** Examine the atlas, patient annotations, and study level metadata.
2. **Ligand receptor feature construction** Derive patient level features from cell population abundance and ligand receptor information.
3. **Feature weighting and filtering** Construct an abundance weighted representation and remove features with zero variance where appropriate.
4. **Signaling program discovery** Apply NMF to identify recurring patterns in the patient level feature matrix.
5. **Program characterization** Examine associations between the resulting programs and tumor characteristics.
6. **Explainable machine learning** Train surrogate models to approximate program scores or assignments and investigate feature contributions.
7. **Evaluation and sensitivity analysis** Assess model performance, program stability, and potential sources of confounding.

The abundance weighted features represent estimated signaling capacity based on the expression derived ligand receptor features and cell population abundance. They are not direct measurements of physical interactions, directional communication, or functional signaling.

## Repository structure

| Directory | Contents |
|---|---|
| `notebooks/` | Exploratory analysis and analysis workflow |
| `data/` | Data files and intermediate analysis outputs |
| `figures/` | Generated figures and visualizations |
| `scripts/` | Standalone analysis and utility scripts |
| `results/` | Empty for now|


## Reproducibility

Reproduction instructions are still being finalized.

They will specify:

- The Python version and required packages.
- How to obtain the source data.
- The order in which notebooks and scripts should be run.
- Which intermediate files are generated at each stage.
- Any computational or memory requirements.


## Limitations

This is an exploratory computational analysis. Any associations found between signaling programs and tumor characteristics do not establish causal relationship.

Model performance on derived program scores also does not establish clinical predictive validity.

The current analysis does not include validated patient survival or recurrence outcomes. Therefore, the results should not be interpreted as a validated breast cancer prognostic model.

## References

1. Chen, A., Kroehling, L., Ennis, C. S., Denis, G. V., & Monti, S. (2026). A highly resolved integrated single-cell atlas of human breast cancers. *NAR Genomics and Bioinformatics, 8*(1), lqaf217. https://doi.org/10.1093/nargab/lqaf217

2. Dimitrov, D., Schäfer, P. S. L., Farr, E., Rodriguez-Mier, P., Lobentanzer, S., Badia-i-Mompel, P., Dugourd, A., Tanevski, J., Ramirez Flores, R. O., & Saez-Rodriguez, J. (2024). LIANA+ provides an all-in-one framework for cell–cell communication inference. *Nature Cell Biology, 26*, 1613–1622. https://doi.org/10.1038/s41556-024-01469-w

3. Lee, D. D., & Seung, H. S. (1999). Learning the parts of objects by non-negative matrix factorization. *Nature, 401*, 788–791. https://doi.org/10.1038/44565

4. Lundberg, S. M., & Lee, S.-I. (2017). A unified approach to interpreting model predictions. In *Advances in Neural Information Processing Systems 30* (pp. 4765–4774). https://papers.nips.cc/paper/7062-a-unified-approach-to-interpreting-model-predictions

5. Virshup, I., Rybakov, S., Theis, F. J., Angerer, P., & Wolf, F. A. (2024). anndata: Access and store annotated data matrices. *Journal of Open Source Software, 9*(101), 4371. https://doi.org/10.21105/joss.04371

6. Wolf, F. A., Angerer, P., & Theis, F. J. (2018). SCANPY: large scale single-cell gene expression data analysis. *Genome Biology, 19*, 15. https://doi.org/10.1186/s13059-017-1382-0

## Project status

Analysis, quality control, and repository preparation are ongoing. Methods, results, and reproduction instructions will be updated as the workflow is audited.
