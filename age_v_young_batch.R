# Read in injured aged timecourse datasets (D0, D3, D7) from 10X v2 and v3
# QC, combine and perform batch correction using LIGER, filter, calculate marker genes, save object, and make plotting dataframe.
# Jacqueline Larouche

library(Seurat) #v3.1.2
library(liger)
library(SeuratWrappers) # v0.1.0
library(dplyr) # v0.8.3
library(cowplot) # v1.0.0
library(ggplot2) # v3.2.1
library(MAST) # v1.12.0
library(tidyverse) # v1.3.0
library(scCATCH) # v2.0

## Load Data
setwd("/nas/homes/jlarouche/Aging/SingleCell/Seurat")
# Aged
MuSC_d0_Aged_1_data <- "10XDatasets/MuSC_d0_Aged_1/mm10" #v2
MuSC_d0_Aged_2_data <- "10XDatasets/MuSC_d0_Aged_2/mm10" #v2
MuSC_d0_Aged_5_data <- "10XDatasets/MuSC_d0_Aged_5/mm10" #v3
MuSC_d0_Aged_6_data <- "10XDatasets/Sample_Geriatric_MuSC_1" #v3

# Young
MuSC_d0_Young_1_data <- "10XDatasets/MuSC_d0_Young_1/mm10" #v2
MuSC_d0_Young_2_data <- "10XDatasets/MuSC_d0_Young_2/mm10" #v2
MuSC_d0_Young_3_data <- "10XDatasets/MuSC_d0_Young_3/mm10" #v3

# TMS
TMS_24mo_58_data <- "10XDatasets/TabulaMurisSenis/24_month/58" #v3
TMS_24mo_59_data <- "10XDatasets/TabulaMurisSenis/24_month/59" #v3
TMS_24mo_60_data <- "10XDatasets/TabulaMurisSenis/24_month/60" #v3
TMS_24mo_61_data <- "10XDatasets/TabulaMurisSenis/24_month/61" #v3

# Create DGE matrices
MuSC_d0_Aged_1_dge <- read10X(sample.dirs = list(MuSC_d0_Aged_1_data), 
                      sample.names = c("MuSCd0Aged1"), min.umis = 300)
MuSC_d0_Aged_2_dge <- read10X(sample.dirs = list(MuSC_d0_Aged_2_data), 
                     sample.names = c("MuSCd0Aged2"), min.umis = 300)
MuSC_d0_Aged_5_dge <- read10X(sample.dirs = list(MuSC_d0_Aged_5_data), 
                     sample.names = c("MuSCd0Aged5"), min.umis = 300)
MuSC_d0_Aged_6_dge <- read10X(sample.dirs = list(MuSC_d0_Aged_6_data), 
                              sample.names = c("MuSCd0Aged6"), min.umis = 300)

MuSC_d0_Young_1_dge <- read10X(sample.dirs = list(MuSC_d0_Young_1_data), 
                               sample.names = c("MuSCd0Young1"), min.umis = 300)
MuSC_d0_Young_2_dge <- read10X(sample.dirs = list(MuSC_d0_Young_2_data), 
                               sample.names = c("MuSCd0Young2"), min.umis = 300)
MuSC_d0_Young_3_dge <- read10X(sample.dirs = list(MuSC_d0_Young_3_data), 
                               sample.names = c("MuSCd0Young3"), min.umis = 300)

TMS_24mo_58_dge <- read10X(sample.dirs = list(TMS_24mo_58_data), 
                     sample.names = c("TMS24mo1"), min.umis = 300)
TMS_24mo_59_dge <- read10X(sample.dirs = list(TMS_24mo_59_data), 
                     sample.names = c("TMS24mo2"), min.umis = 300)
TMS_24mo_60_dge <- read10X(sample.dirs = list(TMS_24mo_60_data), 
                           sample.names = c("TMS24mo3"), min.umis = 300)
TMS_24mo_61_dge <- read10X(sample.dirs = list(TMS_24mo_61_data), 
                           sample.names = c("TMS24mo4"), min.umis = 300)

# Merge into LIGER object
liger10X <- createLiger(list(MuSCd0Aged1 = MuSC_d0_Aged_1_dge, 
                             MuSCd0Aged2 = MuSC_d0_Aged_2_dge,
                             MuSCd0Aged5 = MuSC_d0_Aged_5_dge,
                             MuSCd0Aged6 = MuSC_d0_Aged_6_dge,
                             MuSCd0Young1 = MuSC_d0_Young_1_dge, 
                             MuSCd0Young2 = MuSC_d0_Young_2_dge,
                             MuSCd0Young3 = MuSC_d0_Young_3_dge,
                             TMS24mo1 = TMS_24mo_58_dge,
                             TMS24mo2 = TMS_24mo_59_dge),
                        make.sparse = T, take.gene.union = F, remove.missing = T)
liger10X
# An object of class liger with 12 datasets and 41396 total cells

# continue with other preprocessing steps
ligerex <- liger::normalize(liger10X)
ligerex <- selectGenes(ligerex, var.thresh = 0.1)
ligerex <- scaleNotCenter(ligerex)
saveRDS(ligerex, file = 'aged_v_young_liger')

# Get parameter suggestions
tiff("Figures/Uninjured/suggestK.tiff")
suggestK(ligerex) # try 15
dev.off()

# Perform factorization
ligerex <- optimizeALS(ligerex, k = 15, lambda = 5) 
ligerex <- quantileAlignSNF(ligerex, knn_k = 20, resolution = 1)
# higher resolution for more clusters (default = 1)
# lower knn for more local structure (default = 20)
ligerex <- runUMAP(ligerex)
saveRDS(ligerex, file = 'aged_timecourse_liger_lambda5_k15')
calcAlignment(ligerex) # 0.9645284
calcAgreement(ligerex) # 0.1118749
plots <- plotByDatasetAndCluster(ligerex, return.plots = T)
tiff("Figures/Uninjured/UMAP_dataset_cluster_lambda5_k15.tiff", width = 800, height = 400)
plots[[1]] + plots[[2]]
dev.off()


ligerex <- optimizeNewK(ligerex, k.new = 15, lambda = 3) 
ligerex <- quantileAlignSNF(ligerex, knn_k = 20, resolution = 1)
ligerex <- runUMAP(ligerex)
plots <- plotByDatasetAndCluster(ligerex, return.plots = T)
tiff("Figures/Uninjured/UMAP_dataset_cluster_lambda3_k15.tiff", width = 800, height = 400)
plots[[1]] + plots[[2]]
dev.off()
calcAlignment(ligerex) #0.9645956
calcAgreement(ligerex) #0.1117281

ligerex <- optimizeNewK(ligerex, k.new = 15, lambda = 1) 
ligerex <- quantileAlignSNF(ligerex, knn_k = 20, resolution = 1)
ligerex <- runUMAP(ligerex)
plots <- plotByDatasetAndCluster(ligerex, return.plots = T)
tiff("Figures/Uninjured/UMAP_dataset_cluster_lambda1_k15.tiff", width = 800, height = 400)
plots[[1]] + plots[[2]]
dev.off()
calcAlignment(ligerex) # 0.9028195
calcAgreement(ligerex) # 0.1119011

# Proceed with lambda = 3, k = 15.
ligerex <- optimizeNewK(ligerex, k.new = 15, lambda = 3) 
ligerex <- quantileAlignSNF(ligerex, knn_k = 20, resolution = 1)
ligerex <- runUMAP(ligerex)
saveRDS(ligerex, file = 'aged_v_young_liger')


pdf("Figures/Uninjured/word_clouds.pdf") # important to use PDF to save all pages
plotWordClouds(ligerex)
dev.off()

pdf("Figures/Uninjured/gene_loadings.pdf") # important to use PDF to save all pages
plotGeneLoadings(ligerex)
dev.off()

##### Remove mitochondrial and ribosomal dominated factors
# Factor 3 is largely ribosomal, so remove it from further analysis
factors.use <- c(1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)
ligerex <- quantileAlignSNF(ligerex, dims.use = factors.use)
ligerex <- runUMAP(ligerex, dims.use = factors.use)
plots <- plotByDatasetAndCluster(ligerex, return.plots = T)
tiff("Figures/Uninjured/UMAP_dataset_cluster_liger.tiff", width = 800, height = 400)
plots[[1]] + plots[[2]]
dev.off()
#saveRDS(ligerex, file = 'aged_v_young_liger')

# convert to seurat object for figure plotting
cells.combined <- ligerToSeurat(ligerex, use.liger.genes = T)

##################### Seurat Analysis & Plotting ###########################
ong <- '#ff8500' #3dpi
blu <- '#0689f3' #0dpi
grn <- '#039934' #7dpi
ppl <- '#9B30FF' #young
teal <- '#00FA9A' #aged
gry <- '#9B9999' #TMS

cells.combined <- readRDS('aged_v_young_seurat_v21.rds')
cells.combined@reductions$inmf@assay.used <- 'RNA'
cells.combined@reductions$tsne@assay.used <- 'RNA'
Idents(cells.combined) <- 'orig.ident'
table(cells.combined$orig.ident)
cells.combined.v2 <- subset(cells.combined, idents = c('MuSCd0Aged6','MuSCd0Young4', 'MuSCd0Young5'), invert = TRUE)
table(cells.combined.v2$orig.ident)

# Add metadata for timepoint
Idents(cells.combined) <- 'orig.ident'
current.cluster.ids <- c('MuSCd0Aged1', 'MuSCd0Aged2', 'MuSCd0Aged5', 'MuSCd0Aged6', 'MuSCd0Geriatric', 'MuSCd0Young1', 'MuSCd0Young2', 'MuSCd0Young3', 'MuSCd0Young4', 'MuSCd0Young5','TMS24mo1', 'TMS24mo2')
adj.cluster.ids <- c('Aged1', 'Aged2', 'Aged3', 'Aged3', 'Aged4', 'Young1', 'Young2', 'Young3', 'Young3', 'Young3','TMS24mo1', 'TMS24mo2')
new.cluster.ids <- c('Aged', 'Aged', 'Aged', 'Aged', 'Aged', 'Young', 'Young', 'Young', 'Young', 'Young', 'TMS', 'TMS')
cells.combined$age <- plyr::mapvalues(x = cells.combined$orig.ident, from = current.cluster.ids, to = new.cluster.ids)
cells.combined$adj <- plyr::mapvalues(x = cells.combined$orig.ident, from = current.cluster.ids, to = adj.cluster.ids)
table(cells.combined$age)
# Aged Young   TMS 
# 16227  4557  6466
mean(cells.combined$nCount_RNA) # 5232.825
mean(cells.combined$nFeature_RNA) # 1606.259

factors.use <- c(1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)
cells.combined <- FindNeighbors(cells.combined, reduction = "inmf", dims = factors.use, assay = 'RNA')
cells.combined <- RunUMAP(cells.combined, assay = 'RNA', dims = factors.use, reduction = 'inmf')
cells.combined <- FindClusters(cells.combined, resolution = 1)
saveRDS(cells.combined, 'aged_v_young_seurat_v21.rds')

Idents(cells.combined) <- 'seurat_clusters'
# Visual QC
p1 <- DimPlot(cells.combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(cells.combined, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + NoLegend()
tiff("Figures/Uninjured/umap_seurat_dataset_cluster.tiff", width = 800, height = 400)
p1+p2
dev.off()

p3 <- FeaturePlot(cells.combined, features = c('Pax7', 'Myod1', 'Scx', 'Gpm6b', 'Pecam1', 'Pdgfra', 'Rgs5', 'Ptprc'), min.cutoff = 'q10', max.cutoff = 'q90', ncol = 2)
tiff("Figures/Uninjured/umap_seurat_markergenes_v2.tiff", width = 600, height = 1200)
p3
dev.off()

p4 <- DimPlot(cells.combined, split.by = 'age', group.by = 'orig.ident')#, cols = c(blu, ong, grn))
tiff("Figures/Uninjured/umap_timepoint_split.tiff", width = 900, height = 400)
p4
dev.off()

p5 <- DimPlot(cells.combined, group.by = 'age', cols = c(teal, teal, ppl, gry))
tiff("Figures/Uninjured/umap_uninjured_colored.tiff", width = 400, height = 400)
p5 + theme(text = element_text(family = 'Arial')) + theme_void() + NoLegend()
dev.off()

p6 <- DimPlot(cells.combined, group.by = 'age', split.by = 'age', cols = c(teal, teal, ppl, gry), )
tiff("Figures/Uninjured/umap_uninjured_colored_split.tiff", width = 800, height = 200)
p6 + theme(text = element_text(family = 'Arial')) + theme_void() + NoLegend()
dev.off()

p7 <- DimPlot(cells.combined, group.by = 'Celltype', cols = c(teal, blu, rep(gry,6)))
png("Figures/Uninjured/umap_uninjured_colored_celltype.png", width = 400, height = 400)
p7 + theme(text = element_text(family = 'Arial')) + theme_void() + NoLegend()
dev.off()

Idents(cells.combined) <- 'adj'
my_levels <- c('Young1', 'Young2', 'Young3', 'Aged1', 'Aged2', 'Aged3', 'Aged4', 'TMS24mo1', 'TMS24mo2')
Idents(cells.combined) <- factor(Idents(cells.combined), levels= my_levels)
p7 <- VlnPlot(cells.combined, features = 'nCount_RNA', pt.size = 0, idents = c("Young1", "Young3", "Aged1", "Aged2", "Aged3", "Aged4"))
p7

tiff("Figures/Uninjured/violin_umi_dataset.tiff", width = 600, height = 275)
p7 + theme(text = element_text(family = 'Arial')) + NoLegend()
dev.off()

p8 <- VlnPlot(cells.combined, features = 'nFeature_RNA', pt.size = 0, idents = c("Young1", "Young3", "Aged1", "Aged2", "Aged3", "Aged4"))

tiff("Figures/Uninjured/violin_genes_dataset.tiff", width = 600, height = 275)
p8 + theme(text = element_text(family = 'Arial')) + NoLegend()
dev.off()

# Heatmap of cluster DEGs
Idents(cells.combined) <- 'seurat_clusters'
current.cluster.ids <- as.character(0:11)
new.cluster.ids <- c(rep('MuSC',3),'NMJ', 'MSC', 'SMC', 'Immune', rep('Endo',2), 'MuSC', 'TL', 'MyoN')
cells.combined$Celltype <- plyr::mapvalues(x = cells.combined$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)
Idents(cells.combined) <- 'Celltype'
all.genes <- rownames(cells.combined)
cells.combined <- ScaleData(cells.combined, features = all.genes)
degs.celltype <- FindAllMarkers(cells.combined, only.pos = TRUE, logfc.threshold = 2)
top5 <- degs.celltype %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
head(top5)

png("Figures/Uninjured/heatmap_celltype_degs_top5.png", width = 600, height = 600)
p <- DoHeatmap(cells.combined, features = top5$gene, slot = 'scale.data', group.by = 'Celltype') + scale_fill_gradient2(low="blue", high="red", mid='white')
p
dev.off()


################ Subset MuSCs and UMAP Plot #############################
ligerex_musc <- subsetLiger(ligerex, clusters.use = c(5,10,13, 15))
seurat_musc <- ligerToSeurat(ligerex_musc, use.liger.genes = T)
seurat_musc@reductions$inmf@assay.used <- 'RNA'
seurat_musc@reductions$tsne@assay.used <- 'RNA'
Idents(seurat_musc) <- 'orig.ident'
current.cluster.ids <- c('MuSCd0Aged1', 'MuSCd0Aged2', 'MuSCd0Aged5', 'MuSCd0Aged6', 'MuSCd0Geriatric', 'MuSCd0Young1', 'MuSCd0Young2', 'MuSCd0Young3', 'MuSCd0Young4', 'MuSCd0Young5','TMS24mo1', 'TMS24mo2')
new.cluster.ids <- c('Aged', 'Aged', 'Aged', 'Aged', 'Aged', 'Young', 'Young', 'Young', 'Young', 'Young', 'TMS', 'TMS')
seurat_musc$age <- plyr::mapvalues(x = seurat_musc$orig.ident, from = current.cluster.ids, to = new.cluster.ids)
saveRDS(seurat_musc, "aged_v_young_musc_seurat.rds")

# On Laptop
seurat_musc <- readRDS('aged_v_young_musc_seurat.rds')
seurat_musc <- FindNeighbors(seurat_musc, reduction = "inmf", dims = factors.use, assay = 'RNA')
seurat_musc <- RunUMAP(seurat_musc, assay = 'RNA', dims = factors.use, reduction = 'inmf')
seurat_musc <- FindClusters(seurat_musc, resolution = 1)
saveRDS(seurat_musc, 'aged_v_young_musc_seurat.rds')

Idents(seurat_musc) <- 'orig.ident'
current.cluster.ids <- c('MuSCd0Aged1', 'MuSCd0Aged2', 'MuSCd0Aged5', 'MuSCd0Aged6', 'MuSCd0Geriatric', 'MuSCd0Young1', 'MuSCd0Young2', 'MuSCd0Young3', 'MuSCd0Young4', 'MuSCd0Young5','TMS24mo1', 'TMS24mo2')
adj.cluster.ids <- c('Aged1', 'Aged2', 'Aged3', 'Aged3', 'Aged4', 'Young1', 'Young2', 'Young3', 'Young3', 'Young3','TMS24mo1', 'TMS24mo2')
seurat_musc$adj <- plyr::mapvalues(x = seurat_musc$orig.ident, from = current.cluster.ids, to = adj.cluster.ids)

tiff("Figures/Uninjured/umap_seurat_musc_markergenes.tiff", width = 600, height = 600)
FeaturePlot(seurat_musc, features = c('Pax7', 'Myod1', 'Myf5', 'Myog'), min.cutoff = '0', max.cutoff = 'q90', ncol = 2) + theme(text = element_text(family = 'Arial'))
dev.off()

p6 <- DimPlot(seurat_musc, reduction = 'umap', group.by = 'age', cols = c(teal, teal, ppl, gry))
tiff("Figures/Uninjured/umap_muscs_colored.tiff", width = 400, height = 400)
p6 + theme(text = element_text(family = 'Arial')) + theme_void() + NoLegend()
dev.off()

Idents(seurat_musc) <- 'adj'
my_levels <- c('Young1', 'Young2', 'Young3', 'Young4', 'Young5', 'Aged1', 'Aged2', 'Aged5', 'Aged7', 'TMS24mo1', 'TMS24mo2')
Idents(seurat_musc) <- factor(Idents(seurat_musc), levels= my_levels)
p7 <- VlnPlot(seurat_musc, features = 'Rapsn', idents = c('Aged5', 'Aged7', 'Young3'), pt.size = 0, cols = c(ppl, teal, teal))
png("Figures/Uninjured/violin_muscs_rapsyn.png", width = 500, height = 300)
p7 + theme(text = element_text(family = 'Arial')) + NoLegend()
dev.off()

a <- FindMarkers(seurat_musc, features = c('Rapsn'), ident.1 = 'Aged7', ident.2 = 'Aged5')

# Heatmap of Young vs Aged DEGs
Idents(seurat_musc) <- 'age'
all.genes <- rownames(seurat_musc)
seurat_musc <- ScaleData(seurat_musc, features = all.genes)
degs.musc <- FindMarkers(seurat_musc, ident.1 = 'Young', ident.2 = 'Aged')
top75 <- degs.musc %>% top_n(n = 75, wt = avg_logFC)
top75$gene <- row.names(top75)
head(top75)
bottom75 <- degs.musc %>% top_n(n = -75, wt = avg_logFC)
bottom75$gene <- row.names(bottom75)
head(bottom75)
write.table(degs.musc, file = 'DEGs_MuSC_Aged_Young.csv', sep = ',')

tiff("Figures/Uninjured/heatmap_musc_degs.tif", width = 600, height = 1000)
DoHeatmap(seurat_musc, features = c(top75$gene, bottom75$gene), slot = 'scale.data', group.by = 'age', group.colors = c(ppl, teal, gry)) + scale_fill_gradient2(low="blue", high="red", mid='white')
dev.off()



############### Subset NMJ and Feature Plot #############################
# On Workstation
ligerex_nmj <- subsetLiger(ligerex, clusters.use = c(1,9))
seurat_nmj <- ligerToSeurat(ligerex_nmj, use.liger.genes = T)
seurat_nmj@reductions$inmf@assay.used <- 'RNA'
seurat_nmj@reductions$tsne@assay.used <- 'RNA'
Idents(seurat_nmj) <- 'orig.ident'
current.cluster.ids <- c('MuSCd0Aged1', 'MuSCd0Aged2', 'MuSCd0Aged5', 'MuSCd0Aged6', 'MuSCd0Geriatric', 'MuSCd0Young1', 'MuSCd0Young3', 'MuSCd0Young4', 'MuSCd0Young5','TMS24mo1', 'TMS24mo2')
new.cluster.ids <- c('Aged', 'Aged', 'Aged', 'Aged', 'Aged', 'Young', 'Young', 'Young', 'Young', 'TMS', 'TMS')
seurat_nmj$age <- plyr::mapvalues(x = seurat_nmj$orig.ident, from = current.cluster.ids, to = new.cluster.ids)
saveRDS(seurat_nmj, "aged_v_young_nmj_seurat.rds")

# On Laptop
seurat_nmj <- readRDS('aged_v_young_nmj_seurat.rds')
factors.use <- c(1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)
seurat_nmj <- FindNeighbors(seurat_nmj, reduction = "inmf", dims = factors.use, assay = 'RNA')
seurat_nmj <- RunUMAP(seurat_nmj, assay = 'RNA', dims = factors.use, reduction = 'inmf', n.neighbors = 30L)
FeaturePlot(seurat_nmj, features = c('Pax7'))
seurat_nmj <- FindClusters(seurat_nmj, resolution = 1)
saveRDS(seurat_nmj, 'aged_v_young_nmj_seurat.rds')

p7 <- FeaturePlot(seurat_nmj, features = c('Pax7', 'Myf5', 'Cdh15', 'S100b', 'Cadm1', 'Nrn1'), min.cutoff = '0', max.cutoff = 'q90', ncol = 3)
tiff("Figures/Uninjured/umap_seurat_nmj_markergenes.tiff", width = 900, height = 600)
p7 + theme(text = element_text(family = 'Arial')) & NoAxes()
dev.off()

p8 <- DimPlot(seurat_nmj, reduction = 'umap', group.by = 'age', cols = c(teal, ppl, gry), pt.size = 2) + NoAxes()
png("Figures/Uninjured/umap_seurat_nmj_age.png", res=300, units="in", width=3.5, height=3)
p8
dev.off()

p8 <- DimPlot(seurat_nmj, reduction = "umap", group.by = "age")
p9 <- DimPlot(seurat_nmj, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + NoLegend()
tiff("Figures/Uninjured/umap_seurat_nmj_dataset_cluster.tiff", width = 800, height = 400)
p8+p9
dev.off()

pax7.nmj <- subset(seurat_nmj, subset = Pax7 > 0)

############### Subset Pax7+ cells and Plot MRFs vs celltype ##########################
peach <- '#f7766d'
teal <- '#088992'
nmj_nuclei <- '#651A80'
nmj_pax7 <- '#061A20'
ylw <- '#FCFF33'


png("Figures/Uninjured/umap_seurat_pax7.png", width = 600, height = 500)
FeaturePlot(cells.combined, features = c('Pax7'), min.cutoff = '0', max.cutoff = 'q90', pt.size = 0.05) + theme(text = element_text(family = 'Arial'))
dev.off()

png("Figures/Uninjured/umap_seurat_s100b.png", width = 600, height = 500)
FeaturePlot(cells.combined, features = c('S100b'), min.cutoff = '0', max.cutoff = '2', pt.size = 0.05) + theme(text = element_text(family = 'Arial'))
dev.off()

pax7 <- subset(cells.combined, subset = Pax7 > 0)

Idents(pax7) <- 'Celltype'
pax7 <- subset(pax7, idents = c('MuSC', 'NMJ'))
a <- FindMarkers(pax7, ident.1 = 'NMJ', logfc.threshold = 0.5, only.pos = TRUE)
a$gene_name <- rownames(a)
write.table(a, file = '~/Desktop/nmj_vs_musc_genes.csv', sep = ",")

VlnPlot(pax7, features = c('Agrn', 'Sox2'), pt.size = 0)

pax7.musc <- subset(pax7, idents = c('MuSC'))
pax7.nmj <- subset(pax7, idents = c('NMJ'))

png('Figures/Uninjured/RidgePlot_Pax7_NMJ_MRF.png', width = 700, height = 200)
RidgePlot(pax7.nmj, features = c('Pax7', 'Myf5', 'Myod1'), group.by = 'age', cols = c(teal, ppl)) + NoLegend()
dev.off()

Idents(pax7.nmj) <- 'age'
my_levels <- c('Young', 'Aged')
Idents(pax7.nmj) <- factor(Idents(pax7.nmj), levels= my_levels)

png('Figures/Uninjured/VlnPlot_Pax7_NMJ_MRF.png', width = 1125, height = 250)
VlnPlot(pax7.nmj, features = c('Pax7', 'Myf5', 'Myod1', 'Rapsn', 'S100b'), cols = c(teal), pt.size = 0, idents = 'Aged', ncol = 5) + NoLegend()
dev.off()

Idents(pax7.musc) <- 'age'
my_levels <- c('Young', 'Aged')
Idents(pax7.musc) <- factor(Idents(pax7.musc), levels= my_levels)

png('Figures/Uninjured/VlnPlot_Pax7_MuSC_MRF.png', width = 1125, height = 250)
VlnPlot(pax7.musc, features = c('Pax7', 'Myf5', 'Myod1', 'Rapsn', 'S100b'), cols = c(ppl, teal), pt.size = 0, ncol = 5, idents = c('Young', 'Aged')) + NoLegend()
dev.off()

Idents(cells.combined) <- 'Celltype'
b <- subset(cells.combined, idents = c('NMJ'))
pax7.nmj <- subset(b, subset = Pax7 > 0)
pax7.nmj[['pax7']] <- 'Pax7pos'
other.nmj <- subset(b, subset = Pax7 > 0, invert = TRUE)
other.nmj[['pax7']] <- 'Pax7neg'

pbmc.combined <- merge(pax7.nmj, y = other.nmj)
Idents(pbmc.combined) <- 'pax7'

markers.pax7.nmj <- FindMarkers(pbmc.combined, ident.1 = 'Pax7pos', ident.2 = 'Pax7neg')
write.table(markers.pax7.nmj, file = 'deg_nmj_pax7.csv', sep = ',', row.names = TRUE)

png('Figures/Uninjured/RidgePlot_NMJ_schwann_genes.png', width = 800, height = 5000)
RidgePlot(pbmc.combined, features = c('S100b', 'Rapsn', 'Chodl', 'Nrn1','Chrna1', 'Cadm1', 'Bche', 'Mbp', 'Nrxn1')) + NoLegend()
dev.off()

png('Figures/Uninjured/VlnPlot_musc_schwann_genes.png', width = 600, height = 600)
VlnPlot(pbmc.combined, features = c('Rapsn', 'Chodl', 'Myf5', 'S100b', 'Nrn1', 'Cadm1', 'Mpz', 'Fxyd1', 'Igfbp7'), group.by = 'Celltype', split.by = 'pax7', split.plot=TRUE, pt.size = 0, cols = c(blu, teal))
dev.off()
