### in vitro scRNA Data Analysis

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

#####################################
# Data Integration and Clustering
#####################################

MATRIX = list(
    D1_GFP = "data/count/Rat-Glia-D1-GFP_matrix_10X",
    D1_ND1 = "data/count/Rat-Glia-D1-ND1_matrix_10X",
    D0_Ctrl = "data/count/Rat-Glia-D0_matrix_10X",
    D2_GFP = "data/count/Glia-D2-GFP_matrix_10X",
    D2_ND1 = "data/count/Glia-D2-ND1_matrix_10X",
    D3_GFP = "data/count/Glia-D3-GFP_matrix_10X",
    D3_ND1 = "data/count/Glia-D3-ND1_matrix_10X",
    D5_GFP = "data/count/Glia-D5-GFP_matrix_10X",
    D5_ND1 = "data/count/Glia-D5-ND1_matrix_10X"
)

### Integration
sc.list <- lapply(names(MATRIX), function(x){
	
	# add data info
    Time <- strsplit2(x, "_")[,1]
    Treat <- strsplit2(x, "_")[,2]
    counts <- Seurat::Read10X(MATRIX[[x]])
    
    # Creat Seurat Object
    exp <- CreateSeuratObject(counts = counts, project = x, min.cells = 10, min.feature = 0)
    exp[["percent.mt"]] <- PercentageFeatureSet(exp, pattern = "^Mt")
    exp$Time <- Time
    exp$Treat <- Treat
    
    # QC
    # outlier 
    range <- nCountDensity(exp)
    exp <- isOutlier(exp, method="mad",degree=3, feature="nCount_RNA")
    # remove doublets
    exp <- rmDoublet(exp, outDir = QCDIR, name = x)

    # filtering
    exp.filter <-  subset(exp, nFeature_RNA >= 1000 & percent.mt <= 10 & is_outlier == 'no' & nCount_RNA > range[1] & nCount_RNA < range[2] & is_doublet == 'Singlet')
    return(exp.filter)
})

sc.filter <- merge(sc.list[[1]], sc.list[2:length(sc.list)], add.cell.ids = names(MATRIX))
VlnPlot(sc.filter, features = c("nCount_RNA","nFeature_RNA","percent.mt"),pt.size = 0, ncol = 3,group.by='orig.ident')

# Batch Correction using fastMNN and Clustering
sc.list <- lapply(sc.list, function(x){
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, mean.function = ExpMean, dispersion.function = LogVMR) 
})
sc.vitro <- RunFastMNN(object.list = sc.list)
sc.vitro <- FindNeighbors(sc.vitro, reduction="mnn", dims=1:30)
sc.vitro <- FindClusters(sc.vitro, resolution = 0.5)
sc.vitro <- RunUMAP(sc.vitro, reduction = "mnn", dims = 1:30,spread=0.5, min.dist=0.5, reduction.name = "UMAPmnn")

# DEG
meta = read.csv("inVitro_cell_meta.csv")
sc.vitro <- AddMetaData(object = sc.vitro, metadata = meta)

Idents(sc.vitro) <- "CellType_1"
all.markers <- FindAllMarkers(
    sc.vitro,
    only.pos = T,
    min.pct = 0.25,
    logfc.threshold = 0.25,
    test.use = "MAST",
    assay="RNA" 
)
write.csv(all.markers, "VitroMarkers.csv")

top10.markers <-  all.markers %>% group_by(cluster) %>% top_n(n=10, wt = avg_log2FC)
heatmap_genes <- top10.markers$gene
Doheatmap(sc.vitro, features=heatmap_genes,assay="RNA",group.by="seurat_clusters",size=8) + NoLegend()

cols = c(
        ImA_div='#00008b', ImA = '#8A2BE2', ImA_ND1_hi = '#8B008B',
        Ast_1 = '#6E8B3D', Ast_2 = '#458B00', Ast_3 = '#006400',
        Neu_1 = '#FFC312', Neu_2 =  '#F79F1F', Neu_3 = '#EE5A24', Neu_4 = '#8B0000',
        Others_1 = '#D3D3D3', 'Others_2' = '#BEBEBE')
DimPlot(sc.vitro, group.by = 'CellType_1', label = F, label.size = 5, repel = T, cols = cols, reduction = 'UMAPmnn') 


#####################################
# Trajecctory Analysis
#####################################

### Monocle2 

# use only EGFP cells and Control
sc.filter <- subset(sc.vitro, CellType_1 %!in% c("Others_1", "Others_2")
ND1.cells <- which(GetAssayData(sc.filter,slot="data", assay="RNA")["Neurod1",]>0) %>% names
sc.filter[["is_ND1"]] <- "No"
sc.filter[["is_ND1"]][ND1.cells,] <- "Yes"
sc.filter[["is_ND1"]] <- factor(sc.filter$is_ND1, levels=c("Yes","No"))
GFP.cells <- which(GetAssayData(sc.filter,slot="data", assay="RNA")["EGFP",]>0) %>% names
sc.filter[["is_EGFP"]] <- "No"
sc.filter[["is_EGFP"]][GFP.cells,] <- "Yes"
sc.filter[["is_EGFP"]] <- factor(sc.filter$is_EGFP, levels=c("Yes","No"))
sc.filter <- subset(sc.filter, is_EGFP == 'Yes' | Treat == 'Ctrl')

# define branch signatures to build DDRTree
sc.filter$Branch <- ''
sc.filter@meta.data[which(sc.filter$CellType_2 %in% c('Ast_1', 'Ast_2', 'Ast_3')),]$Branch <- "Ast_potential" 
sc.filter@meta.data[which(sc.filter$CellType_2 %in% c('Neu_1', 'Neu_2', 'Neu_3', 'Neu_4')),]$Branch <- "Neu_potential" 
sc.filter@meta.data[which(sc.filter$CellType_2 %in% c('ImA_div', 'ImA')),]$Branch <- "Initial_state"
sc.filter$BranchA <- ''
sc.filter@meta.data[which(sc.filter$CellType_2 %in% c('Ast_1', 'Ast_2', 'Ast_3', 'Neu_1', 'Neu_2', 'Neu_3', 'Neu_4')),]$BranchA <- "Differentiation" 
sc.filter@meta.data[which(sc.filter$CellType_2 %in% c('ImA_div', 'ImA')),]$BranchA <- "Initial_state"

Idents(sc.filter) <- 'BranchA'
common_deg <- FindMarkers(
    sc.filter, ident.1 = 'Differentiation', ident.2 = 'Initial_state', test.use = 'wilcox'
)
Idents(sc.filter) <- 'Branch'
Ast_potential_deg <- FindMarkers(
    sc.filter, ident.1 = 'Ast_potential', ident.2 = 'Initial_state', tets.use = 'wilcox'
) 
Idents(sc.filter) <- 'Branch'
Neu_potential_deg <- FindMarkers(
    sc.filter, ident.1 = 'Neu_potential', ident.2 = 'Initial_state', tets.use = 'wilcox'
) 
Ast_vs_Neu_deg <- FindMarkers(
    sc.filter, ident.1 = 'Neu_potential', ident.2 = 'Ast_potential', tets.use = 'wilcox'
)

use.genes.for.monocle2 <- unique(
    c(
        rownames(subset(Ast_vs_Neu_deg, p_val_adj <0.05)),
        rownames(subset(Ast_potential_deg, p_val_adj < 0.05)),
        rownames(subset(Neu_potential_deg, p_val_adj < 0.05))
    )
)
common.genes <- rownames(subset(common_deg, p_val_adj < 0.05))
use.genes.for.monocle2.filter <- setdiff(use.genes.for.monocle2, common.genes)

library(monocle)
set.seed(1234)
expr_matrix <- as(as.matrix(sc.filter@assays$RNA@counts), 'sparseMatrix')
p_data <- sc.filter@meta.data
f_data <- data.frame(gene_short_name = row.names(sc.filter), row.names = row.names(sc.filter))
pd <- new('AnnotatedDataFrame', data = p_data)
fd <- new('AnnotatedDataFrame', data = f_data)
cds.sub <- newCellDataSet(
    expr_matrix,
    phenoData = pd, featureData = fd,
    lowerDetectionLimit = 0.25, expressionFamily = negbinomial.size()
)
cds.sub <- estimateSizeFactors(cds.sub)
cds.sub <- estimateDispersions(cds.sub)
cds.sub <- detectGenes(cds.sub, min_expr = 0.1)
cds.sub <- setOrderingFilter(cds.sub, use.genes.for.monocle2.filter)
cds.sub <- reduceDimension(
    cds.sub, 
    max_components = 2, reduction_method = 'DDRTree',
    residualModelFormulaStr = "~orig.ident + Phase"
)
cds.sub <- orderCells(cds.sub)
ggplot(pData(cds.sub), aes(Pseudotime, color = CellType_1, fill=CellType_1)) +
    geom_density(bw=0.5, size=1,alpha=0.5) + 
    theme_classic2()
plot_cell_trajectory(cds.sub, color_by = "Pseudotime",size=1,show_backbone = F, show_branch_points = T) + scale_color_viridis_c()

# Branch Analysis
BEAM_res <- BEAM(cds.sub.filter[use.genes.for.monocle2.filter,], branch_point = 1)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
plot_genes_branched_heatmap(
    cds.sub[rownames(subset(BEAM_res, qval < 1e-4)),],
    branch_point = 1, num_clusters = 6,
    branch_states = c(1,2), branch_labels = c('Ast_potential','Neu_potential'),
    use_gene_short_name = T, show_rownames = T, return_heatmap = T
)


### Monocle3
expr <- GetAssayData(object = sc.filter, slot = slot, assay = assay)
genes <- data.frame(as.character(rownames(expr)))
rownames(genes) <- rownames(expr)
genes <- cbind(genes, genes)
colnames(genes) <- c("GeneSymbol", "gene_short_name")
metadat <- obj@meta.data
cds <- new_cell_data_set(
    expression_data     = expr,
    cell_metadata       = metadat,
    gene_metadata       = genes
)
reducedDimNames(cds) <- "UMAP"
cds <- cluster_cells(cds,verbose = T, cluster_method = "louvain")
cds <- learn_graph(cds)
cds <- order_cells(cds)
cds <- graph_test(cds, neighbor_graph = 'principal_graph', cores = 2, verbose = T)


