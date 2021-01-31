# Read Loom files into Seurat, perform LIGER or CCA integration and filtering, export as h5ad for scVelo (in Python)
# Jacqueline Larouche

library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(dplyr)
library(scCATCH)
library(liger)
sessionInfo()

setwd("~/Documents/UMichigan/Aguilar_Lab/ResearchProjects/SingleCell/scVelo")
filtered_old <- readRDS('aged_timecourse_cca_filtered_v2.rds')

# Define Colors
ong <- '#ff8500' #3dpi
blu <- '#0689f3' #0dpi
grn <- '#039934' #7dpi
ppl <- '#9B30FF' #young
teal <- '#00FA9A' #aged
gry <- '#9B9999' #TMS
ltblu <- '#bdf5fb' #KO
dkblu <- '#1224bd' #Rescue
ylw <- '#FCFF33'

d0 <- '#33ceff'
d3 <- '#68ff33'
d7 <- '#ff6433'
IIb <- '#ca33ff'
IIx <- '#9826BF'
nmj_nuclei <- '#651A80'
nmj_pax7 <- '#061A20'

# Read in loom files, convert to seurat objects, and merge
ldat_d0_old <- ReadVelocity(file = "LoomFiles/Sample_D0_Old.loom")
ldat_d3_old <- ReadVelocity(file = "LoomFiles/Sample_D3_Old.loom")
ldat_d7_old <- ReadVelocity(file = "LoomFiles/Sample_D7_Old_MuSC.loom")
ldat_geriatric <- ReadVelocity(file = "LoomFiles/Sample_Geriatric_MuSC_1.loom")
ldat_24m_myonuclei <- ReadVelocity(file = "LoomFiles/24m_myonuclei.loom")
bm_d0_old <- as.Seurat(x = ldat_d0_old)
bm_d3_old <- as.Seurat(x = ldat_d3_old)
bm_d7_old <- as.Seurat(x = ldat_d7_old)
bm_geri <- as.Seurat(x = ldat_geriatric)
bm_24m <- as.Seurat(x = ldat_24m_myonuclei)

bm_24m[["orig.ident"]] <- "24m.myonuclei"
bm_d0_old[["orig.ident"]] <- "d0.old"
bm_d3_old[["orig.ident"]] <- "d3.old"
bm_d7_old[["orig.ident"]] <- "d7.old"
bm_geri[["orig.ident"]] <- 'geriatric'


merged_bm_old <- merge(bm_d0_old, c(bm_d3_old, bm_d7_old, bm_geri, bm_24m))
merged_bm_old[["RNA"]] <- merged_bm_old[["spliced"]]
table(merged_bm_old$orig.ident)
# 24m.myonuclei        d0.old        d3.old        d7.old     geriatric 
# 9979          6646          7595          6734          7554 

saveRDS(merged_bm_old, file = 'age_nmj_nuclei_raw.rds')

# CCA batch correction
merged.list <- SplitObject(merged_bm_old, split.by = "orig.ident")
merged.list <- lapply(X = merged.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 10000)
})
merged.anchors <- FindIntegrationAnchors(object.list = merged.list, dims = 1:10, k.filter = 50, anchor.features = 5000)
seuratObj <- IntegrateData(anchorset = merged.anchors, dims = 1:10)

saveRDS(seuratObj, 'aged_timecourse_integrated.rds')

# Run the standard work flow for visualization and clustering on the integrated data
DefaultAssay(seuratObj) <- "integrated"
seuratObj <- ScaleData(seuratObj, verbose = TRUE)
seuratObj <- RunPCA(seuratObj, verbose = TRUE)
ElbowPlot(seuratObj, ndims = 50)
seuratObj <- RunUMAP(seuratObj, reduction = "pca", dims = 1:30)
seuratObj <- FindNeighbors(seuratObj, reduction = "pca", dims = 1:30)
seuratObj <- FindClusters(seuratObj, resolution = 0.5)
p1 <- DimPlot(seuratObj, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(seuratObj, reduction = "umap", label = TRUE)
p1+p2

FeaturePlot(seuratObj, features = c('Pax7', 'Myf5', 'Myod1', 'Myog', 'Myh4', 'Myh1'), min.cutoff = 'q10', max.cutoff = 'q90')

# Subset out impurities
VlnPlot(seuratObj, features = c('S100b', 'Scx', 'Rgs5', 'Ptprc', 'Pdgfra', 'Cd74', 'Myod1', 'Myog', 'Pax7', 'Myf5', 'Myh4', 'Myh1'), pt.size=0)
FeaturePlot(seuratObj, features = c('S100b', 'Scx', 'Rgs5', 'Ptprc', 'Pdgfra'), min.cutoff = 'q10', max.cutoff = 'q90')

# Separate NMJ cluster, and remove other cell contaminants (FAPS, Pericytes, Tenocytes, Immune Cells)
nmj_cells <- subset(seuratObj, idents = c('8', '17', '19', '20'), invert = FALSE)
muscs_old <- subset(seuratObj, idents = c('8', '17', '19', '20', '9', '18', '10', '16', '14'), invert = TRUE)

#Add metadata for plotting later
Idents(muscs_old) <- 'orig.ident'
current.cluster.ids <- c('d0.old', 'd3.old', 'd7.old', 'geriatric', "24m.myonuclei")
new.cluster.ids <- c('d0', 'd3', 'd7', 'd0', 'myonuclei')
muscs_old$timepoint <- plyr::mapvalues(x = muscs_old$orig.ident, from = current.cluster.ids, to = new.cluster.ids)
type.cluster.ids <- c('MuSC', 'MuSC', 'MuSC', 'MuSC', 'nuclei')
muscs_old$Celltype <- plyr::mapvalues(x = muscs_old$orig.ident, from = current.cluster.ids, to = type.cluster.ids)

Idents(nmj_cells) <- 'orig.ident'
current.cluster.ids <- c('d0.old', 'd3.old', 'd7.old', 'geriatric', "24m.myonuclei")
time.cluster.ids <- c('NMJ_pax7', 'NMJ_pax7', 'NMJ_pax7', 'NMJ_pax7', 'NMJ_nuclei')
type.cluster.ids <- c('NMJ', 'NMJ', 'NMJ', 'NMJ', 'NMJ_nuclei')
nmj_cells$Celltype <- plyr::mapvalues(x = nmj_cells$orig.ident, from = current.cluster.ids, to = type.cluster.ids)
nmj_cells$timepoint <- plyr::mapvalues(x = nmj_cells$orig.ident, from = current.cluster.ids, to = time.cluster.ids)

# Separate out Pax7+ NMJ cells and NMJ nuclei
nmj_pax7 <- subset(nmj_cells, Pax7 > 0)
table(nmj_pax7$orig.ident)
Idents(nmj_cells) <- 'Celltype'
nmj_nuclei <- subset(nmj_cells, idents = 'NMJ_nuclei')

# Merge nmj cells with MuSCs and myonuclei, save
filtered_old <- merge(muscs_old, c(nmj_pax7, nmj_nuclei))
saveRDS(filtered_old, 'aged_timecourse_cca_filtered_v2.rds')

# Repeat CCA
DefaultAssay(filtered_old) <- 'spliced'
filtered.list <- SplitObject(filtered_old, split.by = "orig.ident")
filtered.list <- lapply(X = filtered.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
filtered.anchors <- FindIntegrationAnchors(object.list = filtered.list, dims = 1:20, k.filter = 100, anchor.features = 2000)
filtered_old <- IntegrateData(anchorset = filtered.anchors, dims = 1:20)
saveRDS(filtered_old, 'aged_timecourse_cca_filtered_v2.rds')
filtered_old <- readRDS('aged_timecourse_cca_filtered_v2.rds')

# Dimensional Reduction and Plot
DefaultAssay(filtered_old) <- "integrated"
filtered_old <- ScaleData(filtered_old, verbose = TRUE)
filtered_old <- RunPCA(filtered_old, verbose = TRUE)
ElbowPlot(filtered_old, ndims = 50)
# UMAP and Clustering
filtered_old <- RunUMAP(filtered_old, reduction = "pca", dims = 1:30)
filtered_old <- FindNeighbors(filtered_old, reduction = "pca", dims = 1:30)
filtered_old <- FindClusters(filtered_old, resolution = 0.5)
p1 <- DimPlot(filtered_old, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(filtered_old, reduction = "umap", label = TRUE) + NoLegend()
p3 <- DimPlot(filtered_old, reduction = "umap", group.by = 'timepoint', label = FALSE)
p1+p2+p3

FeaturePlot(filtered_old, features = c('Pax7', 'Myf5', 'Myod1', 'Myog', 'Myh4', 'Myh1', 'S100b'), min.cutoff = 'q10', max.cutoff = 'q90')

marker.genes <- FindAllMarkers(filtered_old, logfc.threshold = 0.25, only.pos = TRUE)
write.table(marker.genes, file = "Figures/Aged/markergenes.csv", sep = ",")

# Cluster 16 represents immune cells
filtered_old <- subset(filtered_old, idents = c('11', '6', '12', '14', '15', '16'), invert = TRUE) #then repeat UMAP and clustering

# Run scCATCH for celltype annotations, then more filtering
Idents(filtered_old) <- 'seurat_clusters'
clu_markers <- findmarkergenes(filtered_old,
                               species = 'Mouse',
                               cluster = 'All',
                               match_CellMatch = TRUE,
                               cancer = NULL,
                               tissue = 'Muscle',
                               cell_min_pct = 0.25,
                               logfc = 0.25,
                               pvalue = 0.05)
clu_ann <- scCATCH(clu_markers$clu_markers,
                   species = 'Mouse',
                   cancer = NULL,
                   tissue = 'Muscle')
write.table(clu_ann, file = "Figures/Aged/annotations_scCATCH.csv", sep=',')
#then repeat UMAP and clustering

# Save Figures
p1 <- DimPlot(filtered_old, group.by = c("seurat_clusters"), label = TRUE) + NoLegend()
p2 <- DimPlot(filtered_old, group.by = c("orig.ident"), label = FALSE)
p3 <- DimPlot(filtered_old, group.by = c("timepoint"), label = FALSE, cols = c(d0, d3, d7, IIb, ylw, nmj_pax7))
p4 <- FeaturePlot(filtered_old, features = c('Myf5', 'Myod1', 'Myog', 'Myh4', 'Myh1', 'S100b'), min.cutoff = 'q10', max.cutoff = 'q90', ncol=3)
Idents(filtered_old) <- "orig.ident"
p5 <- VlnPlot(filtered_old, features = c('nCount_RNA', 'nFeature_RNA'), pt.size = 0, idents = c('d0.old', 'd3.old', 'd7.old'), cols = c(d0, d3, d7), ncol = 1)
p1
p3
p5

tiff("Figures/Aged/umap_dataset_cluster_filtered_v2.tiff", width = 1200, height = 400)
p1+p2+p3
dev.off()

tiff("Figures/Aged/umap_timepoint_filtered_v2.tiff", width = 600, height = 500)
p3
dev.off()

tiff("Figures/Aged/umap_markergenes_v2.tiff", width = 900, height = 600)
p4
dev.off()

tiff("Figures/Seurat/violin_qc_dataset.tiff", width = 350, height = 555)
p5 + theme(text = element_text(family = 'Arial')) + NoLegend()
dev.off()

table(filtered_old$timepoint,filtered_old$orig.ident)
#               24m.myonuclei d0.old d3.old d7.old geriatric
# d0                     0   4646      0      0      4176
# d3                     0      0   5775      0         0
# d7                     0      0      0   5892         0
# myonuclei           8197      0      0      0         0
# NMJ_nuclei            32      0      0      0         0
# NMJ_pax7               0    646      5     18       256

# Output genes up-regulated in the NMJ cells
Idents(filtered_old) <- 'timepoint'
nmj.degs <- FindMarkers(filtered_old, ident.1 = 'NMJ_pax7', only.pos = TRUE, logfc.threshold = 0.25)
DotPlot(filtered_old, assay = 'spliced', features = c('Nrg1', 'Nrg2', 'Nrg3', 'Nrg4', 'Erbb1', 'Erbb2', 'Erbb3', 'Erbb4'))
FeaturePlot(filtered_old, features = c('Erbb3', 'Ngfr'), min.cutoff = 'q10', max.cutoff = 'q90')
write.table(nmj.degs, file = 'nmj_upregulated_genes.csv', sep = ',')

# Export for scVelo package
DefaultAssay(filtered_old) <- 'RNA' # must only do this right before exporting
SaveH5Seurat(filtered_old, filename = "scAging_Old_v2.h5Seurat", overwrite = TRUE)
Convert("scAging_Old_v2.h5Seurat", dest = "h5ad", overwrite = TRUE)
