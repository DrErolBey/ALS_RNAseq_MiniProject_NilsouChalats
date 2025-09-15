# ALS RNA-seq Mini-Project — Summary

Dataset: GSE234297 (whole blood). Design: 3x3 ALS vs CTRL.

QC:
- PCA: partial separation
- MA and Volcano generated.

Results:
- Significant genes (FDR<0.05): 7
- Top upregulated (examples): CXCR2, KCNJ15, PADI2, NIBAN1, BASP1, LITAF
- Top downregulated (examples): RNA5-8SN2
- Top GO:BP: midbrain development (FDR=0.00041); substantia nigra development (FDR=0.0076); neural nucleus development (FDR=0.011)

Limitations:
- Small n (3×), limited power
- Tissue: whole blood; CNS generalization limited
- Model: ~condition only; no batch covariates
- One DE gene (RNA5-8SN2) shows tiny effect size

Reproducibility:
- Figures: figs/PCA_3v3_GSE234297.pdf, figs/MAplot_3v3_GSE234297.pdf, figs/Volcano_3v3_GSE234297.pdf
- Full table: results/DESeq2_results_3v3_GSE234297.csv
- GO table: results/go_bp.csv
- Session: system/sessionInfo.txt

