# Load Seurat Object containing Pax7+ NMJ cells and Kallisto-aligned bulk RNA-Seq data from Schwann cells. Generate pseudo-bulk dataset from single cell data, using each sample as a replicate and summing reads across cells in each sample. Perform DESeq2 differential expression analysis between schwann cells and MuSCs. Output volcano plot.
# Jacqueline Larouche

# load necessary packages
library(ggplot2)
library(reshape2)
library(data.table)
library(RColorBrewer)
library(DESeq2)
library(zinbwave)
library(EnhancedVolcano)
library(ensembldb)
library(EnsDb.Mmusculus.v79)
library(tximport)
library(Seurat)
library(tibble)
library("AnnotationDbi")

# Read in schwann cell datasets aligned by kallisto
dir <- "~/Documents/UMichigan/Aguilar_Lab/ResearchProjects/SingleCell/Seurat"
setwd(dir)
list.files(dir)
samples <- read.table(file.path(dir, "samples.csv"), header=TRUE)
samples
samples$celltype <- c('Schwann', 'Schwann', 'Schwann', 'Schwann', 'Schwann', 'Schwann', 'Schwann', 'Schwann', 'Schwann', 'Schwann', 'Schwann', 'Schwann')
samples
files <- file.path(dir, "SchwannCellsRNASeq", samples$run, "abundance.h5")
names(files) <- paste0("samples$run", 1:12)
all(file.exists(files)) #should output TRUE here

# Associate transcript IDs with gene names (not IDs because MuSC datasets using gene names)
edb <- EnsDb.Mmusculus.v79
k <- keys(edb, keytype="TXNAME")
tx2gene <- select(edb, k, "GENENAME", "TXNAME")

# Read kallisto .h5 files using tximport
txi.kallisto <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreTxVersion= TRUE, ignoreAfterBar=TRUE)
head(txi.kallisto$counts)
kallisto.list <- matrix(unlist(txi.kallisto), ncol = 12, byrow = TRUE)

# Read in batch-corrected NMJ cluster and generate sudo-bulk datasets with replicates
seurat_nmj <- readRDS('aged_v_young_nmj_seurat.rds')
Idents(seurat_nmj) <- 'seurat_clusters'
filtered_nmj <- subset(seurat_nmj,  ident = '6')
DimPlot(filtered_nmj, group.by = 'orig.ident')
DimPlot(seurat_nmj, reduction = "umap", label = TRUE)

Idents(filtered_nmj) <- 'orig.ident'
bulk.nmj <- NA #initiate dataframe
for (a in levels(filtered_nmj$orig.ident)[1:9]){
  set <- subset(filtered_nmj, ident = a)
  cts <- rowSums(set@assays$RNA@counts)
  df.cts <- as.data.frame(cts)
  bulk.nmj <- as.data.frame(c(bulk.nmj, df.cts))
}
rownames(bulk.nmj) <- rownames(df.cts)
colnames(bulk.nmj) <- levels(filtered_nmj$orig.ident)[1:9]

# Merge pseudo-bulk NMJ dataset with schwann cell dataset
schwann.counts <- txi.kallisto$counts
schwann.cts <- as.data.frame(schwann.counts)
merged.cts <- merge(schwann.cts, bulk.nmj[levels(filtered_nmj$orig.ident)[1:9]],by="row.names",all.x=TRUE)
merged.cts <- column_to_rownames(merged.cts, "Row.names")
merged.cts[is.na(merged.cts)] <- 0
#convert to matrix and round counts up to integer values
merged.cts <- data.matrix(merged.cts)
round.cts <- ceiling(merged.cts)
# save
write.table(round.cts, file = 'merged_round_counts_nmj_schwann.csv', sep = ',')

# Generate dataframe with sample information for DESeq2
samples.merge <- c(colnames(bulk.nmj), levels(samples$run))
celltype <- c(rep('MuSC',9), samples$celltype)
coldata <- data.frame(samples.merge, celltype)
coldata

# Import matrix into DESeq2
dds <- DESeqDataSetFromMatrix(countData = round.cts,
                              colData = coldata,
                              design= ~ celltype)

# Filtering
keep <- rowSums(counts(dds)) >= 5 #may need to reduce this number
dds <- dds[keep,]
dds <- DESeq(dds)
dds
dds$celltype <- relevel(dds$celltype, ref = "Schwann")

# Differential Expression Analysis
res <- results(dds)
res

# Make volcano plot
vol <- EnhancedVolcano(res,
                       lab = rownames(res),
                       x = 'log2FoldChange',
                       y = 'padj',
                       title = 'Schwann vs NMJ MuSC',
                       xlab = bquote(~Log[2]~ ("fold change")),
                       ylab = bquote(~-Log[10]("p-adj")),
                       axisLabSize = 18,
                       titleLabSize = 18,
                       caption = NULL,
                       labvjust = -4,
                       selectLab = c('Pax7', 'Myf5', 'Myod1', 'Dok4', 'Nrxn1', 'Bche', 'Pdgfa', 'Agrn', 'Mbp', 'Mpzl1', 'Mt1', 'Mt2', 'Tnfrsf12a', 'Crlf1', 'Lamtor5'),
                       ylim = c(0,12.5),
                       xlim = c(-7.5, 15),
                       pCutoff = 0.05,
                       pLabellingCutoff = 0,
                       FCcutoff = 0.585,
                       pointSize = 1.0,
                       labSize = 6.0,
                       boxedLabels = TRUE,
                       colAlpha = 0.5,
                       legendPosition = 'right',
                       legendLabSize = 12,
                       legendVisible = FALSE,
                       legendIconSize = 4.0,
                       drawConnectors = TRUE,
                       widthConnectors = 0.2,
                       colConnectors = 'black',
                       col = c('grey', 'grey', 'grey', 'red'))
vol

tiff("Figures/volcano_nmj_v_schwann_padjusted_v3.tiff", width = 1250, height = 750)
vol
dev.off()
