# Read in uninjured young, aged, SOD1 ko and SOD1 rescue datasets from 10X v2 and v3 as well as TMS datasets from 24mo mice
# QC, combine and perform batch correction using LIGER, filter, calculate markergenes, save object, and make plotting dataframe.
# Jacqueline Larouche

##################### Perform on Workstation ###########################
library(Seurat) #v3.1.2
library(liger) #v0.5.0
library(SeuratWrappers) # v0.1.0
library(dplyr) # v0.8.3
library(cowplot) # v1.0.0
library(ggplot2) # v3.2.1

## Load Data
setwd("/nas/homes/jlarouche/Aging/SingleCell/Seurat")
# Aged
MuSC_d0_Aged_1_data <- "10XDatasets/MuSC_d0_Aged_1/mm10" #v2
MuSC_d0_Aged_2_data <- "10XDatasets/MuSC_d0_Aged_2/mm10" #v2
MuSC_d0_Aged_5_data <- "10XDatasets/MuSC_d0_Aged_5/mm10" #v3
MuSC_d0_Aged_6_data <- "10XDatasets/Sample_Old_uninjured" #v3
MuSC_d0_Geriatric_data <- "10XDatasets/Sample_Geriatric_MuSC_1" #v3

# Young
MuSC_d0_Young_1_data <- "10XDatasets/MuSC_d0_Young_1/mm10" #v2
MuSC_d0_Young_2_data <- "10XDatasets/MuSC_d0_Young_2/mm10" #v2
MuSC_d0_Young_3_data <- "10XDatasets/MuSC_d0_Young_3/mm10" #v3
MuSC_d0_Young_4_data <- "10XDatasets/D0_Young_MuSC" #v3
MuSC_d0_Young_5_data <- "10XDatasets/Sample_0dYoungMuSCs" #v3

# SOD1 KO
MuSC_SOD1_KO_4_data <- "10XDatasets/MuSC_SOD1_KO_4/mm10" #v3

# SOD1 Rescue
MuSC_SOD1_rescue_3_data <- "10XDatasets/MuSC_SOD1_rescue_3/mm10"
MuSC_SOD1_rescue_7_data <- "10XDatasets/MuSC_SOD1_rescue_7/mm10"

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
MuSC_d0_Geriatric_dge <- read10X(sample.dirs = list(MuSC_d0_Geriatric_data), 
                              sample.names = c("MuSCd0Geriatric"), min.umis = 300)

MuSC_d0_Young_1_dge <- read10X(sample.dirs = list(MuSC_d0_Young_1_data), 
                               sample.names = c("MuSCd0Young1"), min.umis = 300)
MuSC_d0_Young_2_dge <- read10X(sample.dirs = list(MuSC_d0_Young_2_data), 
                               sample.names = c("MuSCd0Young2"), min.umis = 300)
MuSC_d0_Young_3_dge <- read10X(sample.dirs = list(MuSC_d0_Young_3_data), 
                               sample.names = c("MuSCd0Young3"), min.umis = 300)
MuSC_d0_Young_4_dge <- read10X(sample.dirs = list(MuSC_d0_Young_4_data), 
                               sample.names = c("MuSCd0Young4"), min.umis = 300)
MuSC_d0_Young_5_dge <- read10X(sample.dirs = list(MuSC_d0_Young_5_data), 
                               sample.names = c("MuSCd0Young5"), min.umis = 300)

MuSC_SOD1_KO_4_dge <- read10X(sample.dirs = list(MuSC_SOD1_KO_4_data), 
                               sample.names = c("MuSCSOD1KO4"), min.umis = 300)
MuSC_SOD1_rescue_3_dge <- read10X(sample.dirs = list(MuSC_SOD1_rescue_3_data), 
                               sample.names = c("MuSCSOD1rescue3"), min.umis=300)
MuSC_SOD1_rescue_7_dge <- read10X(sample.dirs = list(MuSC_SOD1_rescue_7_data), 
                               sample.names = c("MuSCSOD1rescue7"), min.umis=300)

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
                             MuSCd0Geriatric = MuSC_d0_Geriatric_dge,
                             MuSCd0Young1 = MuSC_d0_Young_1_dge, 
                             MuSCd0Young2 = MuSC_d0_Young_2_dge,
                             MuSCd0Young3 = MuSC_d0_Young_3_dge,
                             MuSCd0Young4 = MuSC_d0_Young_4_dge,
                             MuSCd0Young5 = MuSC_d0_Young_5_dge,
                             MuSCSOD1KO4 = MuSC_SOD1_KO_4_dge,
                             MuSCSOD1rescue3 = MuSC_SOD1_rescue_3_dge,
                             MuSCSOD1rescue7 = MuSC_SOD1_rescue_7_dge,
                             TMS24mo1 = TMS_24mo_58_dge,
                             TMS24mo2 = TMS_24mo_59_dge),
                        make.sparse = T, take.gene.union = F, remove.missing = T)
liger10X
# An object of class liger with 15 datasets and 59795 total cells

# continue with other preprocessing steps
ligerex <- liger::normalize(liger10X)
ligerex <- selectGenes(ligerex, var.thresh = 0.1)
ligerex <- scaleNotCenter(ligerex)
saveRDS(ligerex, file = 'sod1ko_v_rescue_liger')

# Get parameter suggestions
tiff("Figures/SOD1/suggestK.tiff")
suggestK(ligerex) # try 30
dev.off()

# Perform factorization
ligerex <- optimizeALS(ligerex, k = 30, lambda = 5) 
ligerex <- quantileAlignSNF(ligerex, knn_k = 20, resolution = 1)
# higher resolution for more clusters (default = 1)
# lower knn for more local structure (default = 20)
ligerex <- runUMAP(ligerex)
a <- calcAlignment(ligerex)
b <- calcAgreement(ligerex)
plots <- plotByDatasetAndCluster(ligerex, return.plots = T)
tiff("Figures/SOD1/UMAP_dataset_cluster_lambda5_k30.tiff", width = 800, height = 400)
plots[[1]] + plots[[2]]
dev.off()


ligerex <- optimizeNewK(ligerex, k.new = 30, lambda = 3) 
ligerex <- quantileAlignSNF(ligerex, knn_k = 20, resolution = 1)
ligerex <- runUMAP(ligerex)
plots <- plotByDatasetAndCluster(ligerex, return.plots = T)
tiff("Figures/SOD1/UMAP_dataset_cluster_lambda3_k30.tiff", width = 800, height = 400)
plots[[1]] + plots[[2]]
dev.off()
c <- calcAlignment(ligerex)
d <- calcAgreement(ligerex)

ligerex <- optimizeNewK(ligerex, k.new = 30, lambda = 1) 
ligerex <- quantileAlignSNF(ligerex, knn_k = 20, resolution = 1)
ligerex <- runUMAP(ligerex)
plots <- plotByDatasetAndCluster(ligerex, return.plots = T)
tiff("Figures/SOD1/UMAP_dataset_cluster_lambda1_k30.tiff", width = 800, height = 400)
plots[[1]] + plots[[2]]
dev.off()
e <- calcAlignment(ligerex)
f <- calcAgreement(ligerex)

a #0.9756523
b #0.1034662

c #0.9756762
d #0.103594

e #0.9755847
f #0.1035642

# Proceed with lambda = 3, k = 30.
ligerex <- optimizeNewK(ligerex, k.new = 30, lambda = 3) 
ligerex <- quantileAlignSNF(ligerex, knn_k = 20, resolution = 1)
ligerex <- runUMAP(ligerex)
saveRDS(ligerex, file = 'sod1ko_v_rescue_liger')

pdf("Figures/SOD1/gene_loadings.pdf")
plotGeneLoadings(ligerex)
dev.off()

##### Remove mitochondrial and ribosomal dominated factors
# Factor 3 is largely ribosomal, so remove it from further analysis
# Factor 1 is largely mitochondrial
factors.use <- c(2, 4, 5, 6, 8, 10, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30)
ligerex <- quantileAlignSNF(ligerex, dims.use = factors.use)
ligerex <- runUMAP(ligerex, dims.use = factors.use)
plots <- plotByDatasetAndCluster(ligerex, return.plots = T)
tiff("Figures/SOD1/UMAP_dataset_cluster_liger.tiff", width = 800, height = 400)
plots[[1]] + plots[[2]]
dev.off()
saveRDS(ligerex, file = 'sod1ko_v_rescue_liger')

# convert to seurat object for figure plotting
cells.combined <- ligerToSeurat(ligerex, use.liger.genes = T)
cells.combined@reductions$inmf@assay.used <- 'RNA'
cells.combined@reductions$tsne@assay.used <- 'RNA'

cells.combined <- FindNeighbors(cells.combined, reduction = "inmf", dims = factors.use, assay = 'RNA')
cells.combined <- RunUMAP(cells.combined, assay = 'RNA', dims = factors.use, reduction = 'inmf')
cells.combined <- FindClusters(cells.combined, resolution = 1)
saveRDS(cells.combined, "sod1ko_v_rescue_seurat.rds")

################ Figure Generation ####################
ong <- '#ff8500' #3dpi
blu <- '#0689f3' #0dpi
grn <- '#039934' #7dpi
ppl <- '#9B30FF' #young
teal <- '#00FA9A' #aged
gry <- '#9B9999' #TMS
ltblu <- '#bdf5fb' #KO
dkblu <- '#1224bd' #Rescue

setwd("~/Documents/UMichigan/Aguilar_Lab/ResearchProjects/SingleCell/Seurat")
cells.combined <- readRDS('sod1ko_v_rescue_seurat.rds')

# Add metadata for condition
Idents(cells.combined) <- 'orig.ident'
current.cluster.ids <- c('MuSCd0Aged1', 'MuSCd0Aged2', 'MuSCd0Aged5', 'MuSCd0Aged6', 'MuSCd0Geriatric', 'MuSCd0Young1', 'MuSCd0Young2', 'MuSCd0Young3', 'MuSCd0Young4', 'MuSCd0Young5', 'MuSCSOD1KO4', 'MuSCSOD1rescue3', 'MuSCSOD1rescue7', 'TMS24mo1', 'TMS24mo2')
new.cluster.ids <- c('Aged', 'Aged', 'Aged', 'Aged', 'Aged', 'Young', 'Young', 'Young', 'Young', 'Young','SOD1.KO', 'SOD1.Rescue', 'SOD1.Rescue', 'TMS', 'TMS')
cells.combined$condition <- plyr::mapvalues(x = cells.combined$orig.ident, from = current.cluster.ids, to = new.cluster.ids)
table(cells.combined$condition)
# Aged   Geriatric       Young     SOD1.KO SOD1.Rescue         TMS 
# 15319        7554       12057        6609       11790        6466 
mean(cells.combined$nFeature_RNA) # 1922.987
mean(cells.combined$nCount_RNA) # 6977.408
saveRDS(cells.combined, "sod1ko_v_rescue_seurat.rds")

cells.combined <- FindClusters(cells.combined, resolution = 0.2)
Idents(cells.combined) <- 'seurat_clusters'
# Visual QC
p1 <- DimPlot(cells.combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(cells.combined, reduction = "umap", group.by = "seurat_clusters", label = TRUE, label.size = 8) + NoLegend()
tiff("Figures/SOD1/umap_seurat_dataset_cluster.tiff", width = 800, height = 400)
p1+p2
dev.off()

tiff("Figures/SOD1/umap_seurat_cluster.tiff", width = 400, height = 400)
p2
dev.off()

p3 <- FeaturePlot(cells.combined, features = c('Pax7', 'Myod1', 'S100b', 'Rapsn'), min.cutoff = 'q10', max.cutoff = 'q90', ncol = 4)
tiff("Figures/SOD1/umap_seurat_markergenes_v2.tiff", width = 1150, height = 250)
p3
dev.off()

p4 <- DimPlot(cells.combined, group.by = 'condition', split.by = 'condition', order = c('SOD1.Rescue', 'SOD1.KO', 'Young', 'Geriatric', 'Aged', 'TMS'), cols = c(gry, teal, teal, ppl, ltblu, dkblu))
tiff("Figures/SOD1/umap_uninjured_colored_split.tiff", width = 800, height = 200)
p4 + theme(text = element_text(family = 'Arial')) + theme_void() + NoLegend()
dev.off()

cells.combined <- FindClusters(cells.combined, resolution = 0.1)
Idents(cells.combined) <- 'seurat_clusters'
all.genes <- rownames(cells.combined)
cells.combined <- ScaleData(cells.combined, features = all.genes)
sod1.markers <- FindAllMarkers(cells.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
write.table(sod1.markers, file = 'markergenes_seurat_clusters_sod1_res1.csv', sep = ',')
top20 <- sod1.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
top5 <- c('Nppc', 'Rack1', 'Chodl', 'Crlf1', 'Asb5',
          'Eef1b2', 'Lyz2', 'Ptn', 'Gpm6b', 'Dbi',
          'Cnp', 'Neat1', 'Nfat5', 'Malat1', 'Xist',
          'Kcnq1ot1', 'Myl9', 'Tpm2', 'Acta2', 'Rgs5',
          'Pcp4l1', 'Igkc', 'Cd74', 'H2-Aa', 'H2-Ab1',
          'Tpm1', 'Cxcl14', 'Apod', 'Dcn', 'Lum',
          'Tnfaip6', 'Fmod', 'Thbs4', 'Col1a1', 'Mgp',
          'Fabp4', 'Aqp1', 'Ly6c1', 'Cd36', 'Klf2',
          'Rnd1', 'Icam1', 'Cldn5', 'Il6', 'Ly6a',
          'Cxcl10', 'Hbb-bs', 'Hba-a1', 'Hbb-bt', 'Car2',
          'Gpx1', 'Igfbp6', 'Sfrp5', 'Nbl1', 'Slc2a1',
          'S100a6', 'Mpz', 'Mbp', 'Gatm', 'Fxyd3',
          'Cuedc2', 'Ccl21a', 'Mmrn1', 'Lyve1', 'Nts')
p5 <- DoHeatmap(cells.combined, features = top5, slot = 'scale.data') + scale_fill_gradient2(low="blue", high="red", mid = "white")
tiff("Figures/SOD1/heatmap_cluster_markergenes_scale_data.tiff", width = 1000, height = 700)
p5
dev.off()
DotPlot(cells.combined, features = top5$gene)

Idents(cells.combined) <- 'orig.ident'
my_levels <- c('Young1', 'Young2', 'Young3', 'Young4', 'Young5', 'Aged1', 'Aged2', 'Aged5', 'Aged6', 'Aged7', 'TMS24mo1', 'TMS24mo2')
Idents(cells.combined) <- factor(Idents(cells.combined), levels= my_levels)
p7 <- VlnPlot(cells.combined, features = c('nCount_RNA', 'nFeature_RNA'), pt.size = 0, idents = c('MuSCSOD1KO4', 'MuSCSOD1rescue7'), ncol = 1, cols = c(ltblu, dkblu))

tiff("Figures/SOD1/violin_qc_dataset.tiff", width = 350, height = 475)
p7 + theme(text = element_text(family = 'Arial')) + NoLegend()
dev.off()


#### Find and plot DEGs of SOD1 KO vs Rescue MuSCs
Idents(cells.combined) <- 'seurat_clusters'
sod1.muscs <- subset(cells.combined, idents = c('0', '1', '2', '3', '4', '5', '6', '7', '9', '10', '17', '20', '24'))
DimPlot(sod1.muscs)

Idents(sod1.muscs) <- 'condition'
sod1.muscs.subset <- subset(sod1.muscs, idents = c('SOD1.Rescue', 'SOD1.KO'))
all.genes <- rownames(sod1.muscs.subset)
sod1.muscs.subset <- ScaleData(sod1.muscs.subset, features = all.genes)
degs.musc.sod1 <- FindMarkers(sod1.muscs.subset, ident.1 = 'SOD1.KO', ident.2 = 'SOD1.Rescue')
top75 <- degs.musc.sod1 %>% top_n(n = 75, wt = avg_logFC)
top75$gene <- row.names(top75)
bottom75 <- degs.musc.sod1 %>% top_n(n = -75, wt = avg_logFC)
bottom75$gene <- row.names(bottom75)
write.table(degs.musc.sod1, file = 'DEGs_MuSC_SOD1.csv', sep = ',')

eps("Figures/SOD1/heatmap_sod1_musc_degs.eps", width = 600, height = 1000)
DoHeatmap(sod1.muscs.subset, features = c(top75$gene, bottom75$gene), slot = 'scale.data', group.colors = c(ltblu, dkblu)) + scale_fill_gradient2(low="blue", high="red", mid='white')
dev.off()
