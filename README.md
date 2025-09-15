# ALS RNA-seq Mini-Project - Nilsou Chalats (aka Nilsu Halac)

Data: GSE234297 (whole blood), balanced subset 3 ALS vs 3 CTRL
Method: DESeq2 + apeglm LFC shrinkage
QC: PCA, MA, Volcano; optional GO:BP enrichment

Key results
- Significant genes (FDR < 0.05): 7
- Top upregulated (examples): CXCR2, KCNJ15, PADI2, NIBAN1, BASP1, LITAF
- Top downregulated (examples): RNA5-8SN2
- Top GO:BP: midbrain development (FDR=0.00041); substantia nigra development (FDR=0.0076); neural nucleus development (FDR=0.011)

Files
- Results: results/DESeq2_results_3v3_GSE234297.csv
- GO table: results/go_bp.csv
- Figures: figs/PCA_3v3_GSE234297.pdf, figs/MAplot_3v3_GSE234297.pdf, figs/Volcano_3v3_GSE234297.pdf
- Session: system/sessionInfo.txt

Limitations
- Small n (3x3), limited power
- Tissue: whole blood; CNS generalization limited
- Model: ~condition only; no batch covariates
- One DE gene (RNA5-8SN2) shows tiny effect size

Notes
- IDs are ~98% valid human ENTREZ; symbols provided elsewhere for readability.
- Data source: GEO GSE234297.

Data & Attribution
- Data: GEO GSE234297 (whole blood).
- Notice: GEO data belong to the original authors; not covered by this repo’s license.
- Reuse: Follow GEO terms; cite the study and GSE234297.
- Code/docs: © 2025 Nilsou Chalats, MIT.

