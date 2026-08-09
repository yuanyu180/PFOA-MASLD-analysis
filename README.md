# PFOA-MASLD Analysis Code

Analysis code for the manuscript:

**"Environmental PFOA Exposure and the Risk of Metabolic Dysfunction-Associated Steatotic Liver Disease: An Integrated Computational Toxicology and Multi-Omics Study"**

*Manuscript ID: PONE-D-26-24445 (PLOS ONE)*

This repository contains the R scripts used for all bioinformatics analyses reported in the manuscript: computational toxicology target prediction, transcriptomic data integration and machine learning, single-cell RNA-sequencing (scRNA-seq) analysis, and molecular docking/dynamics data processing.

## Repository Structure

```
├── 01_target_prediction/       # PFOA targets & NAFLD/MASLD disease genes
│   ├── 01_extract_targets_chembl.R    # Retrieve PFOA-related targets from ChEMBL (also STITCH & SwissTargetPrediction)
│   ├── 02_target_venn.R               # Venn diagram of PFOA targets from 3 databases (Fig 1B)
│   ├── 03_disease_gene_venn.R         # Venn diagram of disease genes from GEO/OMIM/GeneCards/TTD (Fig 1F)
│   ├── 04_disease_venn.R              # Disease gene collection
│   ├── 05_intersect_venn.R            # Intersection of PFOA targets & disease genes → 17 core genes (Fig 1G)
│   ├── 06_GO_enrichment.R             # GO enrichment of core genes (Fig 4B annotations)
│   └── 07_KEGG_enrichment.R           # KEGG pathway enrichment
├── 02_geo_differential/
│   └── 01_differential_expression.R   # limma differential expression analysis of GEO training/validation sets
├── 03_machine_learning/
│   ├── 01_ML_functions.R              # Core ML functions (Enet, Lasso, Ridge, Stepglm, SVM, RF, XGBoost, GBM, etc.)
│   ├── 02_run_ML_models.R             # Run 113 algorithm combinations + multi-logistic final model (Fig 2A)
│   ├── refer.methodLists.txt          # List of 113 model combinations
│   ├── data.train.txt                 # Training expression matrix (GSE66676 + GSE89632 + GSE164760, ComBat-corrected)
│   ├── data.test.txt                  # Validation expression matrices (GSE135251, GSE63067)
│   └── results/
│       ├── model.AUCmatrix.txt        # AUC of all 113 models across training & validation cohorts
│       └── model.genes.txt            # Genes selected by each model combination
├── 04_scRNA-seq/
│   └── 01_scRNA_seq_full_pipeline.R   # Full scRNA-seq pipeline (Figs 5-6): QC → Harmony integration → clustering → cell-type annotation → cell proportion → DEG → monocle3 trajectory
└── 05_molecular_docking/
    └── README.md                      # Molecular docking & MD simulation workflow notes
```

## Correspondence Between Code and Manuscript Figures

| Manuscript figure | Analysis | Code |
|---|---|---|
| Fig 1B, 1F, 1G | Target/disease gene Venn diagrams | `01_target_prediction/02,03,05` |
| Fig 1C-D | Batch effect correction (PCA before/after ComBat) | `03_machine_learning` (data preparation) |
| Fig 1E, 2H | Differential expression (volcano plots) | `02_geo_differential/01` |
| Fig 2A-G | 113-model ML framework, ROC, confusion matrices | `03_machine_learning/02_run_ML_models.R` |
| Fig 3A-B | SHAP feature importance & beeswarm | `03_machine_learning/02_run_ML_models.R` |
| Fig 4C-H | GSEA of core genes | `03_machine_learning` + clusterProfiler |
| Fig 5A-G | scRNA-seq clustering, annotation, cell proportions, DEG | `04_scRNA-seq/01` (QC → Harmony integration → clustering → annotation → cell-proportion comparison → DEG) |
| Fig 6A-L | Hepatocyte/myeloid/T-cell subset trajectories | `04_scRNA-seq/01` (monocle3 / monocle trajectory sections) |
| Fig 7A-F | Molecular docking of PFOA with 6 core targets | `05_molecular_docking` (AutoDock) |
| Fig 8A-N | Molecular dynamics simulations (FABP4-PFOA, SHBG-PFOA) | `05_molecular_docking` (GROMACS) |

> **Note on the scRNA-seq script:** `01_scRNA_seq_full_pipeline.R` is the complete working script covering QC, integration, clustering, annotation, cell-proportion comparison, differential expression, and trajectory analysis (monocle/monocle3). A few segments (e.g., the IGF1-related visualization and the CellChat section) originate from the authors' other projects and are not part of this manuscript; they are retained only for completeness of the working record.
>
> For readability, explanatory comments have been stripped from the scripts; the executable code is unchanged and each script maps to the analyses listed in the table above.

## Data Sources

All raw/processed data analyzed in this study are publicly available:

- **GEO datasets** (NCBI Gene Expression Omnibus): GSE66676, GSE89632, GSE164760 (training set); GSE135251, GSE63067 (validation sets); scRNA-seq samples listed in Table 1 of the manuscript (series GSE136103 and GSE174748).
- **Target databases**: ChEMBL, STITCH, SwissTargetPrediction.
- **Disease gene databases**: GEO, OMIM, GeneCards, TTD.
- **Protein structures**: RCSB Protein Data Bank.
- **In vitro experimental data** (qRT-PCR, Western blot, SPR, CCK-8, Oil Red O): provided in the manuscript and its Supporting Information.

## Software Environment

- R ≥ 4.2; key packages: limma, sva, Seurat (v5), harmony, scDblFinder, clusterProfiler, org.Hs.eg.db, ggplot2, ggvenn, glmnet, randomForestSRC, gbm, caret, mboost, e1071, BART, xgboost, pROC, ComplexHeatmap, monocle3.
- Molecular docking: AutoDock Tools / AutoDock Vina; molecular dynamics: GROMACS (see `05_molecular_docking/README.md`).

## License

This code is released under the [MIT License](LICENSE). It is provided as-is for reproducibility of the published analyses; please cite the manuscript if you use it.

## Contact

Corresponding author: Hongzhen Tang (corresponding author of the manuscript).
