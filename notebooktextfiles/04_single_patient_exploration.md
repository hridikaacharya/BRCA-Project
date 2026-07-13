
!pip install liana scanpy anndata numpy pandas matplotlib


    ran successfully


import scanpy as sc
import liana as li
import pandas as pd
import numpy as np





import anndata as ad

f = ad.read_h5ad(
    r"C:\BRCA\BRCA Project\data\IntegratedAtlas.h5ad",
    backed="r"
)

f





    AnnData object with n_obs × n_vars = 621200 × 37389 backed at 'C:\\BRCA\\BRCA Project\\data\\IntegratedAtlas.h5ad'
        obs: 'tissue_ontology_term_id', 'tissue_type', 'assay_ontology_term_id', 'disease_ontology_term_id', 'cell_type_ontology_term_id', 'self_reported_ethnicity_ontology_term_id', 'development_stage_ontology_term_id', 'sex_ontology_term_id', 'donor_id', 'suspension_type', 'grade', 'author_cell_type', 'batch', 'is_primary_data', 'cell_type', 'assay', 'disease', 'sex', 'tissue', 'self_reported_ethnicity', 'development_stage', 'observation_joinid'
        var: 'feature_is_filtered', 'feature_name', 'feature_reference', 'feature_biotype', 'feature_length', 'feature_type'
        uns: 'batch_condition', 'citation', 'default_embedding', 'is_pre_analysis', 'organism', 'organism_ontology_term_id', 'schema_reference', 'schema_version', 'title'
        obsm: 'X_rpca', 'X_umap'





adata.obs["donor_id"].value_counts().head()





    donor_id
    pal_Patient 0177    18812
    pal_Patient 0176    17155
    pal_Patient 0135    14709
    TBB129              13340
    TBB338              12893
    Name: count, dtype: int64





adata.obs["cell_type"].value_counts().head(20)





    cell_type
    malignant cell                                                              242613
    fibroblast                                                                   45702
    macrophage                                                                   44437
    effector memory CD8-positive, alpha-beta T cell                              29592
    effector memory CD4-positive, alpha-beta T cell                              23696
    natural killer cell                                                          23508
    memory B cell                                                                21190
    CD4-positive, CD25-positive, CCR4-positive, alpha-beta regulatory T cell     19802
    naive T cell                                                                 18117
    IgG plasma cell                                                              17824
    exhausted T cell                                                             17359
    CD8-positive, alpha-beta regulatory T cell                                   14182
    conventional dendritic cell                                                  12855
    endothelial cell                                                             11927
    T follicular helper cell                                                     11838
    vein endothelial cell                                                        10318
    vascular associated smooth muscle cell                                        9172
    endothelial cell of vascular tree                                             8749
    cycling T cell                                                                6729
    myofibroblast cell                                                            6689
    Name: count, dtype: int64

donor = adata.obs["donor_id"].unique()[0]

adata_p = adata[adata.obs["donor_id"] == donor].to_memory()

adata_p





    AnnData object with n_obs × n_vars = 6176 × 37389
        obs: 'tissue_ontology_term_id', 'tissue_type', 'assay_ontology_term_id', 'disease_ontology_term_id', 'cell_type_ontology_term_id', 'self_reported_ethnicity_ontology_term_id', 'development_stage_ontology_term_id', 'sex_ontology_term_id', 'donor_id', 'suspension_type', 'grade', 'author_cell_type', 'batch', 'is_primary_data', 'cell_type', 'assay', 'disease', 'sex', 'tissue', 'self_reported_ethnicity', 'development_stage', 'observation_joinid'
        var: 'feature_is_filtered', 'feature_name', 'feature_reference', 'feature_biotype', 'feature_length', 'feature_type'
        uns: 'batch_condition', 'citation', 'default_embedding', 'is_pre_analysis', 'organism', 'organism_ontology_term_id', 'schema_reference', 'schema_version', 'title'
        obsm: 'X_rpca', 'X_umap'





adata_p.obs["cell_type"].value_counts()





    cell_type
    effector memory CD4-positive, alpha-beta T cell                             1310
    naive T cell                                                                1269
    effector memory CD8-positive, alpha-beta T cell                              875
    malignant cell                                                               706
    natural killer cell                                                          360
    memory B cell                                                                322
    T follicular helper cell                                                     313
    CD4-positive, CD25-positive, CCR4-positive, alpha-beta regulatory T cell     310
    fibroblast                                                                   186
    macrophage                                                                   101
    CD8-positive, alpha-beta regulatory T cell                                    89
    endothelial cell of vascular tree                                             77
    cycling T cell                                                                40
    conventional dendritic cell                                                   39
    plasmacytoid dendritic cell                                                   39
    exhausted T cell                                                              36
    endothelial cell                                                              35
    vein endothelial cell                                                         25
    endothelial cell of artery                                                    18
    vascular associated smooth muscle cell                                         8
    myofibroblast cell                                                             6
    cycling macrophage                                                             6
    pericyte                                                                       5
    myeloid dendritic cell                                                         1
    Name: count, dtype: int64


import liana as li

adata_p.var.head()

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
      <th>feature_is_filtered</th>
      <th>feature_name</th>
      <th>feature_reference</th>
      <th>feature_biotype</th>
      <th>feature_length</th>
      <th>feature_type</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>ENSG00000286448</th>
      <td>False</td>
      <td>ENSG00000286448</td>
      <td>NCBITaxon:9606</td>
      <td>gene</td>
      <td>728</td>
      <td>lncRNA</td>
    </tr>
    <tr>
      <th>ENSG00000225880</th>
      <td>False</td>
      <td>LINC00115</td>
      <td>NCBITaxon:9606</td>
      <td>gene</td>
      <td>874</td>
      <td>lncRNA</td>
    </tr>
    <tr>
      <th>ENSG00000230368</th>
      <td>False</td>
      <td>FAM41C</td>
      <td>NCBITaxon:9606</td>
      <td>gene</td>
      <td>919</td>
      <td>lncRNA</td>
    </tr>
    <tr>
      <th>ENSG00000187634</th>
      <td>False</td>
      <td>SAMD11</td>
      <td>NCBITaxon:9606</td>
      <td>gene</td>
      <td>1731</td>
      <td>protein_coding</td>
    </tr>
    <tr>
      <th>ENSG00000188976</th>
      <td>False</td>
      <td>NOC2L</td>
      <td>NCBITaxon:9606</td>
      <td>gene</td>
      <td>1244</td>
      <td>protein_coding</td>
    </tr>
  </tbody>
</table>
</div>

adata_p.var["feature_name"].head()





    ENSG00000286448    ENSG00000286448
    ENSG00000225880          LINC00115
    ENSG00000230368             FAM41C
    ENSG00000187634             SAMD11
    ENSG00000188976              NOC2L
    Name: feature_name, dtype: category
    Categories (37361, object): ['A1BG', 'A1BG-AS1', 'A1CF', 'A2M', ..., 'ZZEF1', 'ZZZ3', 'hsa-mir-1253', 'hsa-mir-423']

adata_p.var_names = adata_p.var["feature_name"]

adata_p.X





    <Compressed Sparse Row sparse matrix of dtype 'float32'
    	with 5871242 stored elements and shape (6176, 37389)>


adata_p.var_names = adata_p.var["feature_name"]




adata_p.var_names[:5]





    Index(['ENSG00000286448', 'LINC00115', 'FAM41C', 'SAMD11', 'NOC2L'], dtype='object', name='feature_name')


import scanpy as sc

sc.pp.normalize_total(adata_p, target_sum=1e4)
sc.pp.log1p(adata_p)

adata_p.raw = adata_p

import scanpy as sc

adata_p.X = adata_p.X.copy()

adata_p.raw = adata_p

import liana.resource as lr

lr_db = lr.select_resource("consensus")
lr_db.head()





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
      <th>ligand</th>
      <th>receptor</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>LGALS9</td>
      <td>PTPRC</td>
    </tr>
    <tr>
      <th>1</th>
      <td>LGALS9</td>
      <td>MET</td>
    </tr>
    <tr>
      <th>2</th>
      <td>LGALS9</td>
      <td>CD44</td>
    </tr>
    <tr>
      <th>3</th>
      <td>LGALS9</td>
      <td>LRP1</td>
    </tr>
    <tr>
      <th>4</th>
      <td>LGALS9</td>
      <td>CD47</td>
    </tr>
  </tbody>
</table>
</div>





import pandas as pd

expr = adata_p.to_df()
expr["cell_type"] = adata_p.obs["cell_type"].values

mean_expr = expr.groupby("cell_type").mean()
mean_expr.shape


    C:\Users\Lenovo\AppData\Local\Temp\ipykernel_24128\4244700143.py:6: FutureWarning: The default of observed=False is deprecated and will be changed to True in a future version of pandas. Pass observed=False to retain current behavior or observed=True to adopt the future default and silence this warning.
    




    (24, 37389)

cell_types = mean_expr.index.tolist()

mean_expr_aligned = mean_expr.copy()
mean_expr_aligned = mean_expr_aligned.loc[cell_types]

scores = []

for _, row in lig_rec.iterrows():
    lig = row["ligand"]
    rec = row["receptor"]

    if lig in mean_expr_aligned.columns and rec in mean_expr_aligned.columns:

        lig_vec = mean_expr_aligned[lig].values.astype(float)
        rec_vec = mean_expr_aligned[rec].values.astype(float)

        if lig_vec.shape == rec_vec.shape:

            score = float((lig_vec * rec_vec).sum())
            scores.append([lig, rec, score])

ccc_df = pd.DataFrame(scores, columns=["ligand", "receptor", "score"])
ccc_df.sort_values("score", ascending=False).head(10)


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
      <th>ligand</th>
      <th>receptor</th>
      <th>score</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>2773</th>
      <td>VIM</td>
      <td>CD44</td>
      <td>58.513568</td>
    </tr>
    <tr>
      <th>78</th>
      <td>HLA-B</td>
      <td>CD3D</td>
      <td>50.588818</td>
    </tr>
    <tr>
      <th>3032</th>
      <td>HMGB1</td>
      <td>CXCR4</td>
      <td>45.188863</td>
    </tr>
    <tr>
      <th>119</th>
      <td>APP</td>
      <td>CD74</td>
      <td>35.022281</td>
    </tr>
    <tr>
      <th>1432</th>
      <td>LGALS1</td>
      <td>PTPRC</td>
      <td>26.619255</td>
    </tr>
    <tr>
      <th>2514</th>
      <td>TIMP1</td>
      <td>CD63</td>
      <td>26.476228</td>
    </tr>
    <tr>
      <th>1431</th>
      <td>LGALS1</td>
      <td>ITGB1</td>
      <td>25.974649</td>
    </tr>
    <tr>
      <th>1623</th>
      <td>MIF</td>
      <td>TNFRSF14</td>
      <td>24.714781</td>
    </tr>
    <tr>
      <th>77</th>
      <td>HLA-B</td>
      <td>CD8A</td>
      <td>23.816556</td>
    </tr>
    <tr>
      <th>170</th>
      <td>HLA-C</td>
      <td>CD8A</td>
      <td>22.594716</td>
    </tr>
  </tbody>
</table>
</div>