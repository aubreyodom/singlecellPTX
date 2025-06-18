---
  title: "liver_scRNAseq"
author: "Xutao Wang"
date: "`r Sys.Date()`"
output:
  html_document:
  toc: true
toc_float: true
code_folding: hide
---
  
  ```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# Import libraries
```{r, message=FALSE, warning=FALSE, results='hide'}
path_to_project <- "/rprojectnb/singlecell/xutao"
condition_control <- "liver_control"
condition_infected <- "liver_infected"
conditions <- c(condition_control, condition_infected)
folder_name <- "CellRangerGex_results/filteredFeatureBarcodeMatrix"

library(dplyr)
library(Seurat)
library(SeuratData)
library(patchwork)
library(ggplot2)

customize_Seurat_FeaturePlot <- function(p, alpha.use = 1, gradient.use = c("yellow", "red"), expression.threshold = 0, is.log1p.transformed = F) {
  
  #### Main function ####
  main_function <- function(p = p, alpha.use = alpha.use, gradient.use = gradient.use, expression.threshold = expression.threshold, is.log1p.transformed = is.log1p.transformed) {
    
    # Order data by gene expresion level
    p$data <- p$data[order(p$data$gene),]
    
    # Define lower limit of gene expression level
    if (isTRUE(is.log1p.transformed)) {
      expression.threshold <- expression.threshold
    } else {
      expression.threshold <- log1p(expression.threshold)
    }
    
    # Compute maximum value in gene expression
    max.exp <- max(p$data$gene)
    
    # Fill points using the gene expression levels
    p$layers[[1]]$mapping$fill <- p$layers[[1]]$mapping$colour
    
    # Define transparency of points
    p$layers[[1]]$mapping$alpha <- alpha.use
    
    # Change fill and colour gradient values
    p <- p + scale_colour_gradientn(colours = gradient.use, guide = F, limits = c(expression.threshold, max.exp), na.value = "grey") +
      scale_fill_gradientn(colours = gradient.use, name = expression(atop(Expression, (log))), limits = c(expression.threshold, max.exp), na.value = "grey") +
      scale_alpha_continuous(range = alpha.use, guide = F)
  }
  
  #### Execution of main function ####
  # Apply main function on all features
  p <- lapply(X = p, alpha.use = alpha.use, gradient.use = gradient.use, 
              expression.threshold = expression.threshold, is.log1p.transformed = is.log1p.transformed,
              FUN = main_function)
  
  # Arrange all plots using cowplot
  # Adapted from Seurat
  # https://github.com/satijalab/seurat/blob/master/R/plotting.R#L1100
  # ncol argument adapted from Josh O'Brien
  # https://stackoverflow.com/questions/10706753/how-do-i-arrange-a-variable-list-of-plots-using-grid-arrange
  cowplot::plot_grid(plotlist = p, ncol = ceiling(sqrt(length(p))))
}

```

# Import data
```{r, eval=FALSE}
liver_list <- lapply(conditions, function(x) {
  path_to_data <- file.path(path_to_project, x, folder_name)
  raw_counts <- Read10X(path_to_data, gene.column = 2, cell.column = 1,
                        unique.features = TRUE, strip.suffix = FALSE)
  meta_data <- data.frame(Status = rep(x, ncol(raw_counts)), 
                          row.names = colnames(raw_counts))
  sobject <- CreateSeuratObject(counts = raw_counts, min.cells = 3, min.features = 200,
                                project = x,
                                meta.data = meta_data)
  sobject <- NormalizeData(sobject)
  sobject <- FindVariableFeatures(sobject, selection.method = "vst", nfeatures = 2000)
  sobject
})

features <- SelectIntegrationFeatures(object.list = liver_list)
```

# Perform integration
```{r, eval=FALSE}
liver_anchors <- FindIntegrationAnchors(object.list = liver_list, 
                                        anchor.features = features)
liver_combined <- IntegrateData(anchorset = liver_anchors)
```

# Perform an integrated analysis
```{r, eval=FALSE}
DefaultAssay(liver_combined) <- "integrated"

liver_combined <- ScaleData(liver_combined, verbose = FALSE)
liver_combined <- RunPCA(liver_combined, npcs = 80, verbose = FALSE)
ElbowPlot(liver_combined, ndims = 50)
```

```{r, eval=FALSE}
set.seed(1)
liver_combined <- RunUMAP(liver_combined, reduction = "pca", dims = 1:30)
liver_combined <- FindNeighbors(liver_combined, reduction = "pca", dims = 1:30)
liver_combined <- FindClusters(liver_combined, resolution = 0.3)
```

```{r, eval=FALSE, fig.width=10}
p1 <- DimPlot(liver_combined, reduction = "umap", group.by = "Status")
p2 <- DimPlot(liver_combined, reduction = "umap", label = "TRUE")
p1 + p2
```

# Save object
```{r}
# saveRDS(liver_combined,
#         file = "/rprojectnb/singlecell/xutao/output/liver_combined.RDS")
```

# Import object
```{r import object, fig.width=8}
liver_combined <- readRDS("/rprojectnb/singlecell/xutao/output/liver_combined.RDS")
# p1 <- DimPlot(liver_combined, reduction = "umap", group.by = "Status")
# p2 <- DimPlot(liver_combined, reduction = "umap", label = "TRUE")
DimPlot(liver_combined, reduction = "umap", label = "TRUE")
```

```{r, fig.width=10, fig.height=5}
DimPlot(liver_combined, reduction = "umap", label = "TRUE", split.by = "Status")
```

# Identify cell type
Prior knowledge was based on [paper1](https://doi.org/10.1016/j.molcel.2019.07.028) and [paper2](https://doi.org/10.1016/j.isci.2021.103233).

## Hepatocytes
Cluster 0, 2, 6, 12
```{r}
FeaturePlot(liver_combined, features = c("Apoc3", "Alb", "Apoa2", "Fabp1"))
```

## Endothelial cells
Cluster 1, 10, 14
```{r, fig.width=8, fig.height=4}
FeaturePlot(liver_combined, features = c("Stab2", "Clec4g"))
```

## Macrophage/Kupffer cells
Cluster 3, 13, 17
```{r, fig.width=8, fig.height=4}
FeaturePlot(liver_combined, features = c("Csf1r", "Clec4f"))
```

## Hepatic stellate cells (HSC)
Cluster 4
```{r, fig.width=8, fig.height=4}
FeaturePlot(liver_combined, features = c("Dcn", "Reln"))
```

## T cells
Cluster 5?
  ```{r}
FeaturePlot(liver_combined, features = c("Cd3g", "Nkg7", "Ccl5"))
```


## Cholangiocytes
Cluster 15
```{r}
FeaturePlot(liver_combined, features = c("Krt7", "Sox9", "Epcam", "Krt8", 
                                         "Krt18", "Aqp1", "Cd24", "Pigr"))
# KRT7, SOX9 EPCAM KRT8 KRT18 AQP1 CD24 PIGR
```

## Dentritic cells (DC)
Cluster 8
```{r}
FeaturePlot(liver_combined, features = c("Irf8"))
```

## B cells
Cluster 11
```{r}
FeaturePlot(liver_combined, features = c("Ebf1", "Cd79a", "Ms4a1"))
```

## Immune cells 
Cd45 was showed as Ptprc
```{r, fig.width=4, fig.height=4}
FeaturePlot(liver_combined, features = c("Ptprc"))
```

# Cell types not detected

## Kupffer cells (KCs)
```{r, eval=FALSE}
FeaturePlot(liver_combined, features = c("Clec4f", "Vsig4", "Lyz1"))
```

## Plasmacytoidp dentritic cells (pDcs)
```{r, eval=FALSE}
FeaturePlot(liver_combined, features = c("Siglech", "Ccr9"))
```

## Neutrophil
```{r, eval=FALSE}
FeaturePlot(liver_combined, features = c("Lyz1", "Retnlg", "Ccl4"))
```

## NK
```{r, eval=FALSE}
FeaturePlot(liver_combined, features = c("Nkg7", "Xcl1", "Ncr1"))
```

## Myofibrolast
```{r, eval=FALSE}
FeaturePlot(liver_combined, features = c("Msh", "Myl7"))
```

## Dividing
```{r, eval=FALSE}
FeaturePlot(liver_combined, features = c("Top2a", "Hist1h2ap", "Hmgb2", "Stmn1"))
```

## Plasma B
```{r, eval=FALSE}
FeaturePlot(liver_combined, features = c("Jchain", "Igha", "Igkc", "Mzb1", "Sec11c"))
```

# Assign cell types
```{r, fig.width=16, fig.height=8}
new.cluster.ids <- c("Hepatocytes", "Endothelial", "Hepatocytes", "Macrophage/Kupffer cells", "Hepatic stellate cells", "T cells", "Hepatocytes", "Hepatocytes", "DC", "moDCs", "Endothelial", "B cells", "Hepatocytes", "Macrophage/Kupffer cells", "Endothelial", "Cholangiocytes", "Unknown", "Macrophage/Kupffer cells")
names(new.cluster.ids) <- levels(liver_combined)
liver_combined <- RenameIdents(liver_combined, new.cluster.ids)
liver_combined$celltype <- Idents(liver_combined)
liver_combined$Status_rename <- gsub("liver_","", liver_combined$Status)
# Remove unknown cells
meta_data_new <- data.frame(cell_name = names(Idents(liver_combined)), 
                            cell_type = Idents(liver_combined)) |> 
  dplyr::inner_join(liver_combined@meta.data |> 
                      dplyr::mutate(cell_name = row.names(liver_combined@meta.data)))
meta_data_new |> 
  dplyr::group_by(Status) |> 
  dplyr::summarise(total = n())
# liver_control 7451 liver_infected 8602
cell_names_sub <- meta_data_new |> 
  dplyr::filter(!cell_type %in% c("Unknown")) |> 
  dplyr::pull(cell_name)
DimPlot(liver_combined[, cell_names_sub], 
        reduction = "umap", label = TRUE, pt.size = 0.5, 
        repel = TRUE, split.by = "Status_rename", label.size = 6)
```

## Figure 1A UMAP of cell type annotation
```{r, fig.width=16, fig.height=8, eval=FALSE}
f1a <- DimPlot(liver_combined[, cell_names_sub], 
               reduction = "umap", label = TRUE, pt.size = 0.5, 
               repel = TRUE, split.by = "Status_rename", label.size = 4) + NoLegend()
ggsave(file.path(path_to_project, "Figures/Figure1/1a.pdf"), f1a, 
       width = 14, height = 7)
```

## Compute percentage
```{r}
# Remove unknown and immune cells
cell_type_prop <- meta_data_new |> 
  dplyr::filter(!cell_type %in% c("Unknown")) |> 
  dplyr::group_by(Status, cell_type) |> 
  dplyr::summarise(n = n()) |> 
  dplyr::mutate(prop = n/sum(n)) |> 
  dplyr::mutate(prop_percent = paste0(round(prop * 100, 2), "%")) |> 
  dplyr::mutate(new_levels = paste0(cell_type, " (", prop_percent,")")) |> 
  dplyr::mutate(Status_rename = gsub(".*\\_", "",Status), 
                Tissue = gsub("\\_.*", "",Status))
```

## Make Stacked bar plot
```{r}
ggplot(cell_type_prop, aes(fill = cell_type, y = prop, x = Status_rename)) + 
  geom_bar(position = "stack", stat = "identity") +
  ggtitle(cell_type_prop$Tissue |> unique()) +
  xlab(NULL) + ylab(NULL) +
  guides(fill = guide_legend(title = "Cell type")) +
  theme_bw() 
```
## Figure 1B Cell proportion
```{r, eval=FALSE}
f1b <- ggplot(cell_type_prop, aes(fill = cell_type, y = prop, x = Status_rename)) + 
  geom_bar(position = "stack", stat = "identity") +
  ggtitle(cell_type_prop$Tissue |> unique()) +
  xlab(NULL) + ylab(NULL) +
  guides(fill = guide_legend(title = "Cell type")) +
  theme_bw() +
  ggtitle("Cell type percentage for liver")
ggsave(file.path(path_to_project, "Figures/Figure1/1b.pdf"), f1b)
write.table(cell_type_prop[,c("Tissue","cell_type", "Status_rename", "prop_percent")],
            file.path(path_to_project, "Figures/Figure1/1b_numbers.txt"), quote = FALSE, row.names = FALSE)
```

## Make label annotation on UMAP
```{r, eval=FALSE}
# Control cells
control_cell_names <- meta_data_new |>
  dplyr::filter(Status == "liver_control", !cell_type %in% c("Unknown")) |> 
  row.names()
liver_control <- liver_combined[, control_cell_names]
new_levels_control <- data.frame(cell_type = levels(Idents(liver_control))) |> 
  dplyr::left_join(meta_data_cell_type |> 
                     dplyr::filter(Status == "liver_control")) |> 
  dplyr::pull(new_levels)
levels(Idents(liver_control)) <- new_levels_control  

# Infected cells
infected_cell_names <- meta_data_new |> 
  dplyr::filter(Status == "liver_infected", !cell_type %in% c("Unknown")) |> 
  row.names()
liver_infected <- liver_combined[, infected_cell_names]
new_levels_infected <- data.frame(cell_type = levels(Idents(liver_infected))) |> 
  dplyr::left_join(meta_data_cell_type |> 
                     dplyr::filter(Status == "liver_infected")) |> 
  dplyr::pull(new_levels)
levels(Idents(liver_infected)) <- new_levels_infected
```

```{r,eval=FALSE, fig.width=20, fig.height=10}
p_control <- DimPlot(liver_control, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 6,
                     repel = TRUE) + 
  NoLegend() + 
  ggtitle("Liver control (46%)")
p_infected <- DimPlot(liver_infected, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 6,
                      repel = TRUE) + 
  NoLegend() + 
  ggtitle("Liver infected (54%)")
p_control + p_infected 
```

# Find genes in the unknown cluster
```{r, eval=FALSE}
gene_avg <- AverageExpression(liver_combined)
# gene_avg_integrated_rank <- gene_avg$integrated |> 
#   apply(1, rank) |> 
#   t()
unknown_genes <- gene_avg$integrated |> 
  as.data.frame() |> 
  dplyr::select(Unknown) |> 
  dplyr::filter(Unknown >= 10) |> 
  row.names()
unknown_genes2 <- gene_avg$RNA |> 
  as.data.frame() |> 
  dplyr::select(Unknown) |> 
  dplyr::filter(Unknown >= 20) |> 
  row.names()
genes_includes <- c("Gpr141", "Gda", "Samhd1", "Trps1", "Pla2g7")
FeaturePlot(liver_combined, 
            features = genes_includes)
# unknown_genes[-which(unknown_genes %in% c("Rbfox1", "Slc8a1", "Gramd1b"))]
```

# Expression of PPAR module (Fgf21 included)
Cannot find UCP2, CD36
```{r, fig.width=8, fig.height=8}
FeaturePlot(liver_combined, features = c("Fgf21", "Lpl", "Angptl4", "Gdf15"), split.by = "Status")
```

## Figure 1C: FGF21 expression among infected
```{r, eval=FALSE}
# Set do.return = TRUE to return ggplot2 object 
p <- FeaturePlot(object = object, features.plot = "Actb", do.return = TRUE) 

# Customize plot by bringing cells with higher expression levels to the front
customize_Seurat_FeaturePlot(p = p, expression.threshold = 0)

f1c <- FeaturePlot(liver_combined[, liver_combined$Status_rename == "infected" & liver_combined$celltype != "Unknown"], features = c("Fgf21"), order = TRUE) +
  scale_colour_gradient(low = "lightgrey", high = "red") +
  ggtitle("FGF21 expression (infected d7)")
ggsave(file.path(path_to_project, "Figures/Figure1/1c_v1.pdf"), f1c, height = 7, width = 7)
```

## Expression of PPAR module among different cells
```{r, fig.width=12, fig.height=4}
VlnPlot(liver_combined[, cell_names_sub], 
        features = c("Fgf21", "Lpl", "Angptl4"), split.by = "Status", 
        group.by = "celltype", pt.size = 0, assay = "RNA") + theme(legend.position = 'right')
```

## Bin genes
```{r, fig.width=8, fig.height=4}
PPAR_module <- list(c("Fgf21", "Lpl", "Angptl4"))
liver_combined_PPAR_module <- AddModuleScore(object = liver_combined, features = PPAR_module, 
                                             name = "PPAR_module_score", assay = "RNA")
FeaturePlot(object = liver_combined_PPAR_module, 
            features = "PPAR_module_score1", split.by = "Status")
```

# Expression of interferon signal markers
TGTP was showed as Tgtp2, Gbp1 was showed as Gbp2b, Ip-10 was showed as Cxcl10
```{r, fig.height=24, fig.width=8}
FeaturePlot(liver_combined, features = c("Iigp1", "Tgtp2", "Gbp2b", "Gbp2", "Gbp3", "Gbp4", "Gbp5", "Cxcl10"), split.by = "Status")
```

## Expression of interferon signal markers among different cells
```{r, fig.width=12, fig.height=12}
VlnPlot(liver_combined[, cell_names_sub], 
        features = c("Iigp1", "Tgtp2", "Gbp2b", "Gbp2", "Gbp3", "Gbp4", "Gbp5", "Cxcl10"), split.by = "Status", 
        group.by = "celltype", pt.size = 0, assay = "RNA", ncol = 3) + theme(legend.position = "right")
```

## Bin genes
```{r, fig.width=8, fig.height=4}
inf_signal_markers_module <- list(c("Iigp1", "Tgtp2", "Gbp2b", "Gbp2", "Gbp3", "Gbp4", "Gbp5", "Cxcl10"))
liver_combined_inf_module <- AddModuleScore(object = liver_combined, features = inf_signal_markers_module, 
                                            name = "inf_signal_markers_module_score", assay = "RNA")
FeaturePlot(object = liver_combined_inf_module, 
            features = "inf_signal_markers_module_score1", split.by = "Status")
```

## Figure 1D: Composite expression of gbp4 and Ip-10
```{r, eval=FALSE}
ifn_gamma <- list(c("Gbp4", "Cxcl10"))
liver_combined_inf_gamma <- AddModuleScore(object = liver_combined, features = ifn_gamma, 
                                           name = "Gbp4_cp_10", assay = "RNA")
f1d <- FeaturePlot(object = liver_combined_inf_gamma[, liver_combined_inf_gamma$Status_rename == "infected" & liver_combined_inf_gamma$celltype != "Unknown"], 
                   features = "Gbp4_cp_101", order = TRUE) +
  scale_colour_gradient(low = "lightgrey", high = "red") +
  ggtitle("Composite GBP and IP-10 expression (infected d7)")

ggsave(file.path(path_to_project, "Figures/Figure1/1d_v1.pdf"), f1d, 
       width = 7, height = 7)
```


# DE test for blue (protruding cluster) vs. red cluster among infected
```{r}
cell_name_want <- meta_data_new |> 
  dplyr::filter(seurat_clusters %in% c(0, 2, 6, 7, 12), 
                Status == "liver_infected") |> 
  dplyr::pull(cell_name)
liver_combine_want <- liver_combined[, cell_name_want]
Idents(liver_combine_want) <- liver_combine_want@meta.data$seurat_clusters
liver_combine_want.response <- FindMarkers(liver_combine_want, 
                                           ident.1 = 7,
                                           ident.2 = c(0, 2, 6, 12), 
                                           verbose = FALSE,
                                           features = c("Fgf21", "Lpl", "Angptl4", "Iigp1", "Tgtp2", "Gbp2b", "Gbp2", "Gbp3", "Gbp4", "Gbp5", "Cxcl10"),
                                           assay = "RNA",
                                           logfc.threshold = 0,
                                           min.cells.feature = 1, min.cells.group = 1, min.pct = 0) |> 
  
  tibble::rownames_to_column("Symbol")
liver_combine_want.response |> 
  DT::datatable()
```

# Expression of IFN-g
```{r, fig.width=8, fig.height=4}
FeaturePlot(liver_combined, features = c("Ifng"), split.by = "Status")
```

## Gene expression among different cell types
```{r}
VlnPlot(object = liver_combined[, cell_names_sub], 
        features = c("Ifng"), split.by = "Status", group.by = "celltype")
```