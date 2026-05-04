# Detailed Guide to Quality Control & Filtering

Quality control is a critical first step in single‑cell RNA‑seq analysis to ensure that downstream results reflect biological variation rather than technical artifacts. This guide expands upon the brief QC section in the main tutorial by explaining common QC metrics, recommended thresholds, and advanced methods for identifying low‑quality cells.

## Why quality control?

Single‑cell libraries often contain droplets with broken cells, empty droplets, or multiple cells captured together. These artifacts can produce high mitochondrial or ribosomal RNA content, few detected genes, or abnormally high library sizes. Filtering them out preserves biological signal and prevents spurious clustering and differential expression.

## Common QC metrics

1. **Number of detected genes (`nFeature_RNA`)** – Low values suggest empty droplets or dying cells; very high values may indicate doublets. The typical approach is to examine the distribution and define thresholds using the median absolute deviation (MAD) or interquartile range. For example, the PBMC dataset often filters cells with <500 genes or >4,000 genes, but this should be adjusted for your data.

2. **Number of counts (`nCount_RNA`)** – The total number of UMIs per cell. Low counts are associated with poor quality; extremely high counts can indicate doublets. Use plots (e.g., scatter plots of `nCount_RNA` vs `nFeature_RNA`) to identify outliers and set thresholds accordingly.

3. **Percentage of mitochondrial reads (`percent.mt`)** – High mitochondrial content suggests stressed or broken cells. TenX Genomics recommends that healthy human cells typically have <10–20 % mitochondrial reads; thresholds should be adjusted by sample and cell type because tissues like heart may naturally have higher mitochondrial expression【803225974757926†L239-L285】.

4. **Percentage of ribosomal reads (`percent.ribo`)** – Elevated ribosomal RNA can also indicate poor quality, though some cell types (e.g., B‑cells) may have high ribosomal content.

## Basic QC workflow

The tutorial uses median absolute deviation as a rough guide:

```r
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = if (any(grepl("^MT-", rownames(obj)))) "^MT-" else "^mt-")
mad_cut <- function(x, nmads=3){
  med <- median(x); m <- mad(x)
  c(max(0, med - nmads * m), med + nmads * m)
}
lim_feat  <- mad_cut(obj$nFeature_RNA)
lim_count <- mad_cut(obj$nCount_RNA)
mt_upper  <- 20
obj <- subset(
  obj,
  subset =
    nFeature_RNA > lim_feat[1] & nFeature_RNA < lim_feat[2] &
    nCount_RNA > lim_count[1] & nCount_RNA < lim_count[2] &
    percent.mt < mt_upper
)
```

Adjust `nmads` and `mt_upper` based on the distribution and known biology. Visualizations (violin plots, scatter plots) help identify appropriate cutoffs.

## Advanced QC: empty droplets, ambient RNA and doublets

- **Empty droplets** – In droplet‑based platforms, many droplets capture no cell but still contain ambient RNA. Use methods like **EmptyDrops** (`DropletUtils` package) to distinguish empty droplets from true cells before creating the Seurat object.

- **Ambient RNA** – RNA released from lysed cells contaminates droplets, causing spurious gene expression. Tools such as **SoupX** or **CellBender** model and remove ambient RNA; they are recommended for complex tissues or droplet assays with high ambient RNA【803225974757926†L239-L285】.

- **Doublets** – Wells containing two cells produce artificially high gene counts and complex gene expression profiles. Use computational doublet detection methods such as scDblFinder (see the dedicated guide).

## Summary & best practices

1. Always inspect distributions of QC metrics and adjust thresholds to the sample and experimental design.
2. Consider the biology: some cell types naturally have high mitochondrial or ribosomal content; avoid overly strict filtering.
3. Use advanced methods like EmptyDrops and SoupX to remove empty droplets and ambient RNA before clustering.
4. After filtering, re‑calculate QC metrics and confirm that the remaining cells show expected distributions.

[Back to main tutorial](../README.md)
