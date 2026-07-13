
# ======================================================
# Notebook 02: Metadata Integration
#
# Goal:
# Audit and integrate patientlevel clinical metadata
# with the integrated breast cancer single-cell atlas.
#
# Output:
# Master patient metadata table for downstream analyses.
# ======================================================




import scanpy as sc
import pandas as pd
import numpy as np

adata = sc.read_h5ad(
    "data/IntegratedAtlas.h5ad",
    backed="r"
)

print(adata)


    AnnData object with n_obs × n_vars = 621200 × 37389 backed at 'data\\IntegratedAtlas.h5ad'
        obs: 'tissue_ontology_term_id', 'tissue_type', 'assay_ontology_term_id', 'disease_ontology_term_id', 'cell_type_ontology_term_id', 'self_reported_ethnicity_ontology_term_id', 'development_stage_ontology_term_id', 'sex_ontology_term_id', 'donor_id', 'suspension_type', 'grade', 'author_cell_type', 'batch', 'is_primary_data', 'cell_type', 'assay', 'disease', 'sex', 'tissue', 'self_reported_ethnicity', 'development_stage', 'observation_joinid'
        var: 'feature_is_filtered', 'feature_name', 'feature_reference', 'feature_biotype', 'feature_length', 'feature_type'
        uns: 'batch_condition', 'citation', 'default_embedding', 'is_pre_analysis', 'organism', 'organism_ontology_term_id', 'schema_reference', 'schema_version', 'title'
        obsm: 'X_rpca', 'X_umap'
    



import os

os.listdir("data/integratedatlas_supplemental_files")





    ['SuppTable1.xlsx',
     'SuppTable2 (epi).xlsx',
     'SuppTable3 (imm).xlsx',
     'SuppTable4 (strom).xlsx',
     'supp_figures.pdf']





patient_metadata = pd.read_excel(
    "data/integratedatlas_supplemental_files/SuppTable1.xlsx"
)

patient_metadata.head()





<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>donor_id</th>
      <th>Basal</th>
      <th>Her2</th>
      <th>LumA</th>
      <th>LumB</th>
      <th>Normal</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>BC17086-12-Tumor</td>
      <td>0.254329</td>
      <td>0.429654</td>
      <td>0.0</td>
      <td>0.316017</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>1</th>
      <td>BC17086-24-Tumor</td>
      <td>0.479841</td>
      <td>0.333790</td>
      <td>0.0</td>
      <td>0.186369</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>2</th>
      <td>BC17086-25-Tumor</td>
      <td>0.044725</td>
      <td>0.482418</td>
      <td>0.0</td>
      <td>0.472856</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>3</th>
      <td>BC17086-35-Tumor</td>
      <td>0.000000</td>
      <td>0.567416</td>
      <td>0.0</td>
      <td>0.432584</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>4</th>
      <td>BC258-Tumor</td>
      <td>0.000000</td>
      <td>0.487568</td>
      <td>0.0</td>
      <td>0.512432</td>
      <td>0.0</td>
    </tr>
  </tbody>
</table>
</div>





print(patient_metadata.shape)
patient_metadata.columns.tolist()


    (138, 6)
    




    ['donor_id', 'Basal', 'Her2', 'LumA', 'LumB', 'Normal']





patient_metadata.isna().sum().sort_values(ascending=False)





    donor_id    0
    Basal       0
    Her2        0
    LumA        0
    LumB        0
    Normal      0
    dtype: int64





patient_metadata.head(20)





<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>donor_id</th>
      <th>Basal</th>
      <th>Her2</th>
      <th>LumA</th>
      <th>LumB</th>
      <th>Normal</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>BC17086-12-Tumor</td>
      <td>0.254329</td>
      <td>0.429654</td>
      <td>0.000000</td>
      <td>0.316017</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>1</th>
      <td>BC17086-24-Tumor</td>
      <td>0.479841</td>
      <td>0.333790</td>
      <td>0.000000</td>
      <td>0.186369</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>2</th>
      <td>BC17086-25-Tumor</td>
      <td>0.044725</td>
      <td>0.482418</td>
      <td>0.000000</td>
      <td>0.472856</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>3</th>
      <td>BC17086-35-Tumor</td>
      <td>0.000000</td>
      <td>0.567416</td>
      <td>0.000000</td>
      <td>0.432584</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>4</th>
      <td>BC258-Tumor</td>
      <td>0.000000</td>
      <td>0.487568</td>
      <td>0.000000</td>
      <td>0.512432</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>5</th>
      <td>BC302-Tumor</td>
      <td>0.000000</td>
      <td>0.341736</td>
      <td>0.000000</td>
      <td>0.658264</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>6</th>
      <td>BC389-Tumor</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.708314</td>
      <td>0.117440</td>
      <td>0.174246</td>
    </tr>
    <tr>
      <th>7</th>
      <td>BC392-Tumor</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.581742</td>
      <td>0.000000</td>
      <td>0.418258</td>
    </tr>
    <tr>
      <th>8</th>
      <td>BC393-Tumor</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.555948</td>
      <td>0.444052</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>9</th>
      <td>BC394-Tumor</td>
      <td>0.075120</td>
      <td>0.493940</td>
      <td>0.000000</td>
      <td>0.430940</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>10</th>
      <td>BC397-Tumor</td>
      <td>0.000000</td>
      <td>0.268810</td>
      <td>0.000000</td>
      <td>0.731190</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>11</th>
      <td>BC401-Tumor</td>
      <td>0.000000</td>
      <td>0.569668</td>
      <td>0.000000</td>
      <td>0.430332</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>12</th>
      <td>BC419-Tumor</td>
      <td>0.080594</td>
      <td>0.453399</td>
      <td>0.000000</td>
      <td>0.466006</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>13</th>
      <td>BC428-Tumor</td>
      <td>0.000000</td>
      <td>0.453507</td>
      <td>0.000000</td>
      <td>0.546493</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>14</th>
      <td>BIOKEY-1</td>
      <td>1.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>15</th>
      <td>BIOKEY-10</td>
      <td>0.089028</td>
      <td>0.570509</td>
      <td>0.000000</td>
      <td>0.340463</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>16</th>
      <td>BIOKEY-11</td>
      <td>0.000000</td>
      <td>0.632520</td>
      <td>0.000000</td>
      <td>0.367480</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>17</th>
      <td>BIOKEY-12</td>
      <td>0.000000</td>
      <td>0.180664</td>
      <td>0.000000</td>
      <td>0.819336</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>18</th>
      <td>BIOKEY-13</td>
      <td>0.148432</td>
      <td>0.851568</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>19</th>
      <td>BIOKEY-14</td>
      <td>0.113012</td>
      <td>0.000000</td>
      <td>0.156638</td>
      <td>0.000000</td>
      <td>0.730350</td>
    </tr>
  </tbody>
</table>
</div>

# Predicted subtype for each patient
patient_metadata["Predicted_Subtype"] = patient_metadata[
    ["Basal", "Her2", "LumA", "LumB", "Normal"]
].idxmax(axis=1)

patient_metadata.head()
<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>donor_id</th>
      <th>Basal</th>
      <th>Her2</th>
      <th>LumA</th>
      <th>LumB</th>
      <th>Normal</th>
      <th>Predicted_Subtype</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>BC17086-12-Tumor</td>
      <td>0.254329</td>
      <td>0.429654</td>
      <td>0.0</td>
      <td>0.316017</td>
      <td>0.0</td>
      <td>Her2</td>
    </tr>
    <tr>
      <th>1</th>
      <td>BC17086-24-Tumor</td>
      <td>0.479841</td>
      <td>0.333790</td>
      <td>0.0</td>
      <td>0.186369</td>
      <td>0.0</td>
      <td>Basal</td>
    </tr>
    <tr>
      <th>2</th>
      <td>BC17086-25-Tumor</td>
      <td>0.044725</td>
      <td>0.482418</td>
      <td>0.0</td>
      <td>0.472856</td>
      <td>0.0</td>
      <td>Her2</td>
    </tr>
    <tr>
      <th>3</th>
      <td>BC17086-35-Tumor</td>
      <td>0.000000</td>
      <td>0.567416</td>
      <td>0.0</td>
      <td>0.432584</td>
      <td>0.0</td>
      <td>Her2</td>
    </tr>
    <tr>
      <th>4</th>
      <td>BC258-Tumor</td>
      <td>0.000000</td>
      <td>0.487568</td>
      <td>0.0</td>
      <td>0.512432</td>
      <td>0.0</td>
      <td>LumB</td>
    </tr>
  </tbody>
</table>
</div>





patient_metadata["Predicted_Subtype"].value_counts()






    Predicted_Subtype
    LumB      40
    Basal     37
    LumA      29
    Her2      21
    Normal    11
    Name: count, dtype: int64





for f in [
    "SuppTable2 (epi).xlsx",
    "SuppTable3 (imm).xlsx",
    "SuppTable4 (strom).xlsx",
]:
    df = pd.read_excel(f"data/integratedatlas_supplemental_files/{f}")
    print(f)
    print(df.shape)
    print(df.columns.tolist()[:20])
    print("-"*80)


    SuppTable2 (epi).xlsx
    (100, 7)
    ['gene', 'coef', 'mean', 'pval', 'fdr', 'group', 'cluster']
    
    SuppTable3 (imm).xlsx
    (48, 6)
    ['cluster', 'cell_type', 'sub_type', 'Cell Ontology', 'notes', 'action']
    
    SuppTable4 (strom).xlsx
    (17, 6)
    ['cluster', 'cell_type', 'sub_type', 'Cluster', 'cell ontology', 'notes']    