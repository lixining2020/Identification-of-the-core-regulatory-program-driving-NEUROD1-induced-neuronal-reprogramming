# Invivo Analysis

library(Seurat)
library(presto)
library(patchwork)
library(clustree)
library(limma)
library(ggsignif)
library(ggpubr)
library(ggsci)
library(ggpattern)
library(SeuratWrappers)
library(phylogram)
library(factoextra)


MATRIX <- list(
	E16.5 = "count/Rat-cortex-E16-5_matrix_10X",
	P2 = "count/Rat-cortex-P2_matrix_10X"
)

sc.list <- lapply(names(MATRIX), function(x){

    counts <- Seurat::Read10X(MATRIX[[x]])
    
    # Creat Seurat Object
    exp <- CreateSeuratObject(counts = counts, project = x, min.cells = 10, min.feature = 0)
    exp[["percent.mt"]] <- PercentageFeatureSet(exp, pattern = "^Mt")
    
    exp <- isOutlier(exp, method="mad",degree=3, feature="nCount_RNA")
    exp <- rmDoublet(exp, outDir = QCDIR, name = x)

    # filtering
    exp.filter <-  subset(exp, nFeature_RNA >= 1000 & percent.mt <= 10 & is_outlier == 'no' & is_doublet == 'Singlet')
    return(exp.filter)
})

sc.list <- lapply(sc.list, function(x){
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, mean.function = ExpMean, dispersion.function = LogVMR) 
    x <- RenameCells(x, add.cell.id = unique(x$orig.ident))
})
features <- SelectIntegrationFeatures(object.list = sc.list)
anchors <- FindIntegrationAnchors(sc.list, anchor.features = features)
sc.combined <- IntegrateData(anchorset = anchors)
DefaultAssay(sc.combined) <- 'integrated'
sc.combined <- ScaleData(sc.combined, vars.to.regress = c('nCount_RNA'))
sc.combined <- RunPCA(sc.combined)
sc.combined <- FindNeighbors(sc.combined, reduction="pca", dims=1:30)
sc.combined <- FindClusters(sc.combined, resolution = 0.8)
sc.combined <- RunUMAP(sc.combined, reduction = "pca", dims = 1:30,spread=0.5, min.dist=0.5, return.model = T)

DimPlot(sc.combined, reduction = 'umap', label = T, label.size = 6)

# Mapping in vitro data to in vivo

map.anchors <- FindTransferAnchors(
    reference = sc.combined,
    query = sc.vitro,
    dims = 1:50, reference.reduction = "pca",
    k.filter = 150
)

predictions <- TransferData(
    anchorset = map.anchors,
    refdata = sc.combined$CellType, dims = 1:50
)

sc.vitro$Predictions.Combined <- predictions$predicted.id