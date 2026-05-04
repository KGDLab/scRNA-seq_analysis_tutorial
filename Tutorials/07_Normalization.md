# Normalization & Feature Selection: SCTransform vs LogNormalize

Normalization transforms raw UMI counts into values that are comparable between cells while preserving biological variation. It accounts for differences in sequencing depth and technical noise. Feature selection identifies highly variable genes that drive downstream analyses like PCA and clustering.

## Traditional normalization (LogNormalize)

A common approach uses **size‑factor normalization**:

1. **Sum counts per cell** to obtain the library size.
2. **Divide each cell’s counts by its total counts** to normalize sequencing depth.
3. **Multiply** by a scale factor (e.g., 10,000) to bring values into a comparable range.
4. **Add a pseudocount** (often 1) and apply log transformation (`log1p`)【77240560611190†L40-L104】.

This method is simple but may not adequately model the mean–variance relationship of UMI counts. Highly expressed genes dominate the variance, limiting the detection of subtle signals. Variants include methods like scran and SCnorm that compute cell‑specific size factors.

## SCTransform

**SCTransform** models counts using a regularized negative binomial regression. It fits a generalized linear model for each gene, using sequencing depth and other covariates (e.g., percent.mt or cell cycle scores) as predictors, and returns Pearson residuals as normalized values【970707360990613†L55-L139】. Key advantages include:

- **Eliminates the need** for pseudocount and log transformation; the Pearson residuals are already variance‑stabilized【970707360990613†L55-L139】.
- **Effective removal** of technical effects like sequencing depth and mitochondrial content, enabling regression of multiple covariates in one step【970707360990613†L55-L139】.
- **Improved detection** of variable genes: SCTransform typically retains more genes and uses a larger number of principal components than log normalization【970707360990613†L55-L139】.
- **Default in Seurat v5**: SCTransform is now Seurat’s recommended normalization method; log normalization is provided mainly for educational purposes【970707360990613†L55-L139】.

SCTransform can be accelerated using the `glmGamPoi` backend; this speeds up negative binomial regression without affecting results【970707360990613†L55-L139】.

### Usage

```r
obj <- SCTransform(
  obj,
  variable.features.n = 3000,
  vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"),
  verbose = FALSE
)
DefaultAssay(obj) <- "SCT"
```

Regressing variables like `percent.mt` or cell‑cycle scores ensures that these confounders do not drive variability in downstream analyses.

## Choosing between methods

- **SCTransform** is recommended for most datasets because it models count distributions and removes technical noise more effectively.

- **LogNormalize** remains useful for small or simple datasets and for teaching purposes. When using LogNormalize, follow up with `FindVariableFeatures()` to identify variable genes, then `ScaleData()` to z‑score genes and regress unwanted covariates.

## Feature selection

Seurat’s `FindVariableFeatures()` (or the variable features returned by SCTransform) identifies genes with high cell‑to‑cell variation. By default, SCTransform selects the top 3,000 genes; you can adjust `variable.features.n` based on dataset size. More features may be beneficial for complex tissues but increase computational cost.

[Back to main tutorial](../README.md)
