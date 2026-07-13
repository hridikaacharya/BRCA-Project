
import sys
print(sys.executable)
```

    C:\Users\Lenovo\miniconda3\envs\brca\python.exe
    



import scanpy as sc

print(sc.__version__)

```

    1.11.5
    

    C:\Users\Lenovo\AppData\Local\Temp\ipykernel_5272\1266513164.py:3: FutureWarning: `__version__` is deprecated, use `importlib.metadata.version('scanpy')` instead
      print(sc.__version__)
    



import os

print(os.getcwd())
```

    C:\BRCA\BRCA Project
    



import os

print(os.listdir("data"))
```

    ['IntegratedAtlas.h5ad', 'integratedatlas_supplemental_files', 'wu_spatial']
    



adata = sc.read_h5ad("data/IntegratedAtlas.h5ad", backed="r")

print(adata)
```

    AnnData object with n_obs × n_vars = 621200 × 37389 backed at 'data\\IntegratedAtlas.h5ad'
        obs: 'tissue_ontology_term_id', 'tissue_type', 'assay_ontology_term_id', 'disease_ontology_term_id', 'cell_type_ontology_term_id', 'self_reported_ethnicity_ontology_term_id', 'development_stage_ontology_term_id', 'sex_ontology_term_id', 'donor_id', 'suspension_type', 'grade', 'author_cell_type', 'batch', 'is_primary_data', 'cell_type', 'assay', 'disease', 'sex', 'tissue', 'self_reported_ethnicity', 'development_stage', 'observation_joinid'
        var: 'feature_is_filtered', 'feature_name', 'feature_reference', 'feature_biotype', 'feature_length', 'feature_type'
        uns: 'batch_condition', 'citation', 'default_embedding', 'is_pre_analysis', 'organism', 'organism_ontology_term_id', 'schema_reference', 'schema_version', 'title'
        obsm: 'X_rpca', 'X_umap'
    



adata.obs.columns.tolist()
```




    ['tissue_ontology_term_id',
     'tissue_type',
     'assay_ontology_term_id',
     'disease_ontology_term_id',
     'cell_type_ontology_term_id',
     'self_reported_ethnicity_ontology_term_id',
     'development_stage_ontology_term_id',
     'sex_ontology_term_id',
     'donor_id',
     'suspension_type',
     'grade',
     'author_cell_type',
     'batch',
     'is_primary_data',
     'cell_type',
     'assay',
     'disease',
     'sex',
     'tissue',
     'self_reported_ethnicity',
     'development_stage',
     'observation_joinid']





adata.obs["donor_id"].nunique()
```




    138





adata.obs["batch"].value_counts()
```




    batch
    pal_2021          224823
    wu_natgen_2021     99876
    tietscher_2023     97301
    bassez_2021        78549
    wang_2024          48553
    qian_2020          40665
    liu_2023           27490
    gao_2021            3943
    Name: count, dtype: int64





# Cells per patient
cells_per_patient = adata.obs["donor_id"].value_counts()
print(cells_per_patient.describe())
```

    count      138.000000
    mean      4501.449275
    std       3535.196537
    min        165.000000
    25%       1774.500000
    50%       3789.500000
    75%       6045.500000
    max      18812.000000
    Name: count, dtype: float64
    



# Top 20 patients by cell count
cells_per_patient.head(20)


```




    donor_id
    pal_Patient 0177      18812
    pal_Patient 0176      17155
    pal_Patient 0135      14709
    TBB129                13340
    TBB338                12893
    pal_Patient 0178      12779
    pal_Patient 0167      12213
    pal_Patient 0337      11046
    BC392_Tumor           10786
    pal_Patient 0114       9851
    pal_Patient 0025       9763
    Patient4_T7            9752
    TBB011                 9476
    pal_Patient 0173       9338
    BC393_Tumor            9292
    TBB214                 9008
    pal_Patient 0554       8891
    wu_natgen_CID4471      8608
    TBB035                 8512
    wu_natgen_CID44971     7980
    Name: count, dtype: int64





# Cell type distribution (top 30)
adata.obs["cell_type"].value_counts().head(30)
```




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
    pericyte                                                                      5870
    cycling macrophage                                                            4811
    mast cell                                                                     4272
    endothelial cell of artery                                                    4108
    plasmacytoid dendritic cell                                                   2590
    myeloid dendritic cell                                                        1334
    capillary endothelial cell                                                     916
    cycling stromal cell                                                           581
    CD4-positive helper T cell                                                     419
    Name: count, dtype: int64





# Author annotations (top 30)
adata.obs["author_cell_type"].value_counts().head(30)

```




    author_cell_type
    epi_0              51241
    epi_1              40441
    epi_2              38378
    CD8_Tem            29592
    epi_low_ent        29062
    CD4_Tem            23696
    epi_3              20006
    CD4_Treg           19802
    CD4_Naive          18117
    Plasma_IgG         17824
    epi_4              17743
    CD8_Tex            17359
    CAFs (COL11A1+)    14727
    CAFs (LAMP5)       14387
    CD8_Isg            14182
    NK                 13821
    Mac/Mono           13488
    cDC2               12011
    CD4_Tfh            11838
    Endo Imm.          11158
    Endo Ven.          10318
    epi_5              10120
    VSMC                9172
    Endo Capil.         8749
    CAFs (PI16+)        8425
    NK_rNK              7065
    epi_6               7036
    Mac_RTM_LA          6776
    T_Prolif.           6729
    Myofibroblast       6689
    Name: count, dtype: int64





# Grade
adata.obs["grade"].value_counts(dropna=False)

```




    grade
    3.0    314868
    nan    171051
    2.0    103007
    1.0     32274
    Name: count, dtype: int64





# Missing values
adata.obs.isna().sum()

```




    tissue_ontology_term_id                     0
    tissue_type                                 0
    assay_ontology_term_id                      0
    disease_ontology_term_id                    0
    cell_type_ontology_term_id                  0
    self_reported_ethnicity_ontology_term_id    0
    development_stage_ontology_term_id          0
    sex_ontology_term_id                        0
    donor_id                                    0
    suspension_type                             0
    grade                                       0
    author_cell_type                            0
    batch                                       0
    is_primary_data                             0
    cell_type                                   0
    assay                                       0
    disease                                     0
    sex                                         0
    tissue                                      0
    self_reported_ethnicity                     0
    development_stage                           0
    observation_joinid                          0
    dtype: int64





# Number of unique patients per study
adata.obs.groupby("batch")["donor_id"].nunique()

```

    C:\Users\Lenovo\AppData\Local\Temp\ipykernel_5272\2308953735.py:2: FutureWarning: The default of observed=False is deprecated and will be changed to True in a future version of pandas. Pass observed=False to retain current behavior or observed=True to adopt the future default and silence this warning.
      adata.obs.groupby("batch")["donor_id"].nunique()
    



batch
    bassez_2021       31
    gao_2021           4
    liu_2023           4
    pal_2021          33
    qian_2020         14
    tietscher_2023    12
    wang_2024         14
    wu_natgen_2021    26
    Name: donor_id, dtype: int64