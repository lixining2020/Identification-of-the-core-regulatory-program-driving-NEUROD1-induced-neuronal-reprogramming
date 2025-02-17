# ATAC Analysis

library(ArchR)
library(Signac)
library(AnnotationDbi)
library(org.Rn.eg.db)
library(TxDb.Rnorvegicus.UCSC.rn7.refGene)
library(parallel)

load("rn7_liftoff_mm10NcbiRefSeq_ArchR_annotations.rda")

fragments <- 'Neurod1-ATAC/outs/fragments.tsv.gz'

ArrowFiles <- createArrowFiles(
    inputFiles = fragments,
    sampleNames = 'ND1-ATAC',
    filterTSS = 4,
    filterFrags = 1000,
    addTileMat = T,
    addGeneScoreMat = T,
    geneAnnotation = geneAnnotation,
    genomeAnnotation = genomeAnnotation
)

doubScores <- addDoubletScores(
    input = ArrowFiles,
    k = 10,
    knnMethod = "UMAP",
    LSIMethod = 1
)

proj <- ArchRProject(
    ArrowFiles = ArrowFiles,
    outputDirectory = "ND1-ATAC",
    copyArrows = T,
    geneAnnotation = geneAnnotation,
    genomeAnnotation = genomeAnnotation
)

proj <- filterDoublets(proj)

proj <- addIterativeLSI(
    ArchRProj = proj,
    useMatrix = "TileMatrix",
    name = "IterativeLSI",
    iterations = 2,
    clusterParams = list(
        resolution = c(0.2), sampleCells = 7000, n.start = 10
    ),
    varFeatures = 25000, dimsToUse = 1:30
)

proj <- addClusters(
    input = proj,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.8
)

proj <- addUMAP(
    ArchRProj = proj,
    reducedDims = "IterativeLSI",
    name = 'UMAP',
    nNeighbors = 30,
    minDist = 0.6, metric = "cosine",
    force = TRUE
)

proj <- addGeneIntegrationMatrix(
    ArchRProj = proj,
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = sc.vitro,
    addToArrow = FALSE,
    groupRNA = "CellType_1",
    UMAPParams = list(n_neighbors = 30, min_dist = 0.6, metric = "cosine", verbose =FALSE),
    nameCell = "predictedCell",
    nameGroup = "predictedGroup",
    nameScore = "predictedScore",
    threads = 5
)

plotEmbedding(proj, colorBy = "cellColData", name = "predictedGroup")

proj <- addReproduciblePeakSet(
    ArchRProj = proj, 
    groupBy = "predictedGroup", 
    geneAnnotation = geneAnnotation,
    genomeAnnotation = genomeAnnotation,
    genomeSize = 2.75e9
)

proj <- addPeakMatrix(proj)

proj <- addMotifAnnotations(
    ArchRProj = proj, 
    motifSet = "JASPAR2020", 
    name = "Motif",
    force = T,
    cutOff = 0.0001
)

markerPeaks <- getMarkerFeatures(
    proj,
    useMatrix = "PeakMatrix",
    groupBy = "predictedGroup",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

plotBrowserTrack(
    ArchRProj = proj, 
    groupBy = "predictedGroup",
    geneSymbol = "Neurod1",
    features =  getMarkers(markerPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 1", returnGR = TRUE),
    upstream = 50000, downstream = 50000
)



