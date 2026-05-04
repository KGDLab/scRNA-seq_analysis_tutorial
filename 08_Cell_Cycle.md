# Cell Cycle Scoring & Regression

Cell cycle heterogeneity can introduce strong variation in single‑cell RNA‑seq datasets, obscuring biological signals. For example, stem cells and proliferating immune cells express distinct sets of S‑phase and G2/M‑phase genes. If not corrected, cell cycle can drive principal components and lead to clustering by cell cycle rather than cell type.

## Scoring cell cycle phases

Seurat provides a list of canonical S and G2/M markers (`Seurat::cc.genes.updated.2019`). Use `CellCycleScoring()` to compute two scores per cell:

```r
s.genes   <- Seurat::cc.genes.updated.2019$s.genes
g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes
obj <- CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
```

This adds three columns to `meta.data`: `S.Score`, `G2M.Score`, and `Phase` (assigned based on which score is highest). The default `set.ident = TRUE` sets cell identity to the cell cycle phase.

## Regressing out cell cycle effects

To remove cell cycle variation, regress the scores during scaling or SCTransform. There are two approaches【981570600885476†L187-L256】:

1. **Regress both S and G2/M scores**: This completely removes cell cycle signal but may blur distinctions between subpopulations of proliferating cells. Use this when cell cycle is purely nuisance variation.

   ```r
   obj <- SCTransform(obj, vars.to.regress = c("S.Score","G2M.Score"), verbose = FALSE)
   ```

   Alternatively, with log normalization:

   ```r
   obj <- NormalizeData(obj)
   obj <- FindVariableFeatures(obj)
   obj <- ScaleData(obj, vars.to.regress = c("S.Score","G2M.Score"))
   ```

2. **Regress the difference (`S.Score` – `G2M.Score`)**: This retains separation between cycling and non‑cycling cells (e.g., stem vs progenitor) while removing differences among proliferating cells【981570600885476†L187-L256】. Set `vars.to.regress = "PhaseDiff"` where `PhaseDiff = obj$S.Score - obj$G2M.Score`.

   ```r
   obj$PhaseDiff <- obj$S.Score - obj$G2M.Score
   obj <- SCTransform(obj, vars.to.regress = c("percent.mt", "PhaseDiff"), verbose = FALSE)
   ```

## Best practices

1. Use cell cycle scoring to explore whether your dataset is strongly driven by cell cycle.
2. Decide whether to regress cell cycle based on biological questions. If you study proliferation, avoid full regression.
3. When integrating multiple datasets, regress cell cycle either before or during integration to prevent cell cycle from dominating batch correction.

[Back to main tutorial](../README.md)
