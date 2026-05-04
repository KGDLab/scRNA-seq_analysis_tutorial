# Doublet Detection

Droplet‑based single‑cell experiments encapsulate cells into droplets via Poisson loading. Although the aim is one cell per droplet, a fraction of droplets capture two or more cells (doublets or multiplets). Doublets can confound clustering and downstream analyses because they exhibit expression signatures of multiple cell types.

## Why detect doublets?

Doublets often exhibit high UMI counts and increased numbers of expressed genes. Heterotypic doublets (containing cells from different lineages) are more problematic because they can form artificial clusters. Identifying and removing doublets improves clustering accuracy and prevents misinterpretation of marker genes.

## scDblFinder

The **scDblFinder** package is a fast and accurate doublet detection tool within the Bioconductor ecosystem. It synthesizes insights from existing algorithms and introduces improvements such as adaptive estimation of doublet rates and heterotypic/homotypic scoring【550163560382897†L220-L305】. Benchmarking shows that scDblFinder outperforms alternatives and detects most heterotypic doublets while maintaining speed【550163560382897†L220-L305】.

### Usage

1. **Convert Seurat object to `SingleCellExperiment`**:

```r
library(scDblFinder)
library(SingleCellExperiment)
sce <- as.SingleCellExperiment(obj)
```

2. **Run scDblFinder**:

```r
sce <- scDblFinder(sce)
table(sce$scDblFinder.class)
```

The result labels each cell as “singlet” or “doublet” in `sce$scDblFinder.class`.

3. **Transfer labels back to Seurat and filter**:

```r
obj$doublet <- sce$scDblFinder.class
obj <- subset(obj, subset = doublet == "singlet")
```

By default, scDblFinder estimates the expected doublet rate based on the number of cells; you can override this using the `dbr` parameter. For experiments mixing samples or conditions, use `samples=...` to perform per‑sample doublet detection.

## Other doublet detection methods

Alternative approaches include:

- **DoubletFinder** (Seurat workflow) – Generates artificial doublets and uses k‑nearest neighbours to score observed cells; slower than scDblFinder and requires parameter tuning.
- **scds** – Provides multiple scoring methods (`cxds`, `bcds`, `hybrid`) but may be less sensitive for heterotypic doublets.
- **Scrublet** (Python) – Widely used in Scanpy; uses simulated doublets and nearest‑neighbor classification.

## Best practices

1. Run doublet detection after normalization and feature selection but before clustering.
2. Examine `nCount_RNA` and `nFeature_RNA` distributions to manually identify potential doublets.
3. Adjust expected doublet rates based on cell loading concentration and platform (typical rates range from 5–10 %).
4. After filtering doublets, re‑run key steps such as PCA and clustering.

[Back to main tutorial](../README.md)
