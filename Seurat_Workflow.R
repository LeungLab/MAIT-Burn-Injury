library(Seurat)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(limma)
library(viridis)
library(MAST)

#Importing and setting up matrix for analysis
#The first file must clear the first 7 lines of csv and add sample tag column to end of the matrix

Abseq_1 <- read.csv(file = '~/Desktop/Leung_Lab/Seurat/Skin_BD/Skin_Seurat/DL2_DBEC_MPC_Abseq_SampleTag.csv', sep = ',', header = TRUE, row.names = 1, check.names = FALSE)
Abse_1 <- Abseq_1[!grepl("Multiplet", Abseq_1$Sample_Name),]
Abse_1 <- Abseq_1[!grepl("Undetermined", Abseq_1$Sample_Name),]
Abse_1 <- Abse_1[!grepl("Multiplet", Abse_1$Sample_Name),]
Abseq_1 <- Abse_1
Abseq_1P <- t(Abseq_1[, str_detect(string = colnames(Abseq_1), pattern = 'pAbO')])
Abseq_1RNA <- (Abseq_1[, !str_detect(string = colnames(Abseq_1), pattern = 'pAbO')])
Abseq_1RNA <- t(Abseq_1RNA[, !str_detect(string = colnames(Abseq_1RNA), pattern = 'Sample')])
dim(Abseq_1P)
dim(Abseq_1RNA)
Sample_Tag1 <- t(Abseq_1[, str_detect(string = colnames(Abseq_1), pattern = 'Sample_Name')])
dim(Sample_Tag1)
P_names <- sapply(X = str_split(string = rownames(Abseq_1P), pattern = '\\|'),
                  FUN = function(x) ifelse(x[1] == 'CD197',
                                           paste(paste(x[1], x[2], sep = '|'), x[3], sep = '_'),
                                           paste(x[1], x[2], sep = '|')))
rownames(Abseq_1P) <- paste0(P_names)
dim(Abseq_1RNA)
dim(Abseq_1P)
SkinBD <- CreateSeuratObject(counts = Abseq_1RNA, min.cells = 1)
SkinBD@meta.data[, "Condition"] <- as.vector(Sample_Tag1)
table(SkinBD@meta.data$Condition)
AG1  AG2  AG5   B1   B2   B5
2190 3267 3404 2468 3000 1091
SkinBD@meta.data$State[SkinBD@meta.data$Condition == "AG1"] <- "NonBurn_1"
SkinBD@meta.data$State[SkinBD@meta.data$Condition == "AG2"] <- "NonBurn_2"
SkinBD@meta.data$State[SkinBD@meta.data$Condition == "AG5"] <- "NonBurn_3"
SkinBD@meta.data$State[SkinBD@meta.data$Condition == "B1"] <- "Burn_1"
SkinBD@meta.data$State[SkinBD@meta.data$Condition == "B2"] <- "Burn_2"
SkinBD@meta.data$State[SkinBD@meta.data$Condition == "B5"] <- "Burn_3"
table(SkinBD@meta.data$State)
SkinBD@meta.data$Cond[SkinBD@meta.data$Condition == "AG1"] <- "Non_Burn"
SkinBD@meta.data$Cond[SkinBD@meta.data$Condition == "AG2"] <- "Non_Burn"
SkinBD@meta.data$Cond[SkinBD@meta.data$Condition == "AG5"] <- "Non_Burn"
SkinBD@meta.data$Cond[SkinBD@meta.data$Condition == "B1"] <- "Burn"
SkinBD@meta.data$Cond[SkinBD@meta.data$Condition == "B2"] <- "Burn"
SkinBD@meta.data$Cond[SkinBD@meta.data$Condition == "B5"] <- "Burn"
table(SkinBD@meta.data$Cond)
Burn_1    Burn_2    Burn_3 NonBurn_1 NonBurn_2 NonBurn_3
2468      3000      1091      2190      3267      3404
Assays(SkinBD)
[1] "RNA"
Abseq_assay <- CreateAssayObject(counts = Abseq_1P)
SkinBD[["AbSeq"]] <- Abseq_assay
Assays(SkinBD)
[1] "RNA"   "AbSeq"
DefaultAssay(SkinBD)
[1] "RNA"
DefaultAssay(SkinBD) <- "AbSeq" # If want to switch to AbSeq
DefaultAssay(SkinBD)
[1] "AbSeq"

#Analysis of CD3 RNA
VlnPlot(SkinBD, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
plot2 <- FeatureScatter(SkinBD, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2
SkinBD <- subset(SkinBD, subset = nFeature_RNA > 20 & nFeature_RNA < 2500)
SkinBD <- NormalizeData(SkinBD, normalization.method = "LogNormalize", scale.factor = 10000)
SkinBD <- FindVariableFeatures(SkinBD, selection.method = "vst", nfeatures = 300)
top10 <- head(VariableFeatures(SkinBD), 10)
plot1 <- VariableFeaturePlot(SkinBD)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
all.genes <- rownames(SkinBD)
SkinBD <- ScaleData(SkinBD, features = all.genes)
SkinBD <- RunPCA(SkinBD, features = VariableFeatures(object = SkinBD))
print(SkinBD[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(SkinBD, dims = 1:2, reduction = "pca")
DimPlot(SkinBD, reduction = "pca")
DimHeatmap(SkinBD, dims = 1, cells = 500, balanced = TRUE)
SkinBD <- JackStraw(SkinBD, num.replicate = 100, prop.freq = 0.05)
SkinBD <- ScoreJackStraw(SkinBD, dims = 1:20)
JackStrawPlot(SkinBD, dims = 1:15)
SkinBD <- FindNeighbors(SkinBD, dims = 1:10)
SkinBD <- FindClusters(SkinBD, resolution = 0.8, verbose = FALSE)
SkinBD <- RunUMAP(SkinBD, dims = 1:15)
DimPlot(SkinBD, label = TRUE, label.size = 15, group.by = "seurat_clusters", pt.size = 2)
DimPlot(SkinBD, label = FALSE, group.by = "Cond", cols = c("red", "black"), pt.size = 2)
Idents(object = SkinBD) <- "State"
NonBurn.markers <- FindMarkers(SkinBD, ident.1 = c("NonBurn_1", "NonBurn_2", "NonBurn_3"), min.pct = 0.25)
Burn.markers <- FindMarkers(SkinBD, ident.1 = c("Burn_1", "Burn_2", "Burn_3"), min.pct = 0.25)
head(NonBurn.markers, n = 10)
VlnPlot(SkinBD, features = c("ITGAE","CD69","S1PR1","ITGAE","CCR7", "SELL"))
FeaturePlot(SkinBD, features = c("CD3","ITGAE","CD69","S1PR1","ITGAE","CCR7", "SELL"), blend = TRUE)
FeaturePlot(SkinBD, features = c("CD3","ITGAE","CD69","S1PR1","ITGAE","CCR7", "SELL"), cols = plasma)
FeaturePlot(SkinBD, features = c("TBX21","RORC","FOXP3","BCL6","EOMES","GATA3","RORA","RUNX3", "ZBTB16"), cols = plasma)
SkinBD.marker.Ab <- FindAllMarkers(SkinBD, only.pos = TRUE, min.pct = 0.0001, logfc.threshold = 0.0001)
SkinBD.marker.cond %>%
    group_by("Cond") %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(SkinBD, features = top50$gene, draw.lines = FALSE, group.colors = c("red","black")) + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))
DoHeatmap(SkinBD, features = SkinBD.marker.Ab$gene, draw.lines = FALSE, group.colors = c("red","black")) + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))
FeaturePlot(SkinBD, features = c("IFNG","IL32","TNF","IL2","NAMPT","CSF3","TNFSF10","FASLG", "CSF2")) & scale_colour_gradientn(colours = rev(RColorBrewer::brewer.pal(n =11, name = "RdBu")))
plasma <- viridis(9, direction = 1, option = "C")
FeaturePlot(SkinBD, features = c("IFNG","IL32","TNF","IL2","NAMPT","CSF3","TNFSF10","FASLG", "CSF2"), cols = plasma)
DoHeatmap(SkinBD, features = top10$gene, draw.lines = FALSE) + scale_fill_viridis(6, direction = 1, option = "H")
VlnPlot(SkinBD.CD4, features = c("ITGAE","CD69","S1PR1","CCR7","SELL"), stack = TRUE, group.by = "Cond", fill.by ="ident", cols = c("red", "black"))
VlnPlot(SkinBD, features = c("ITGAE","CD69","S1PR1","CCR7","SELL"), cols = c("red", "black"), group.by = "Cond", flip = TRUE)
VlnPlot(SkinBD, features = c("IFNG","IL32","TNF","IL2","NAMPT","CSF3","TNFSF10","FASLG", "CSF2"), stack = TRUE, group.by = "Cond", fill.by ="ident", cols = c("red", "black"))
DotPlot(SkinBD, features = c("ITGAE","CD69","S1PR1","CCR7","SELL","SELPG"), group.by = "Cond") + RotatedAxis()+ scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))
FeaturePlot(SkinBD_final, features = c("ITGAE","CD69","S1PR1","CCR7","SELL"), cols = plasma)


#CD4 RNA analysis
SkinBD.CD4 <- subset(x = SkinBD, subset = CD4 > 2)
VizDimLoadings(SkinBD.CD4, dims = 1:2, reduction = "pca")
DimPlot(SkinBD.CD4, reduction = "pca")
DimHeatmap(SkinBD.CD4, dims = 1, cells = 500, balanced = TRUE)
SkinBD.CD4 <- JackStraw(SkinBD.CD4, num.replicate = 100, prop.freq = 0.05)
SkinBD.CD4 <- ScoreJackStraw(SkinBD.CD4, dims = 1:20)
JackStrawPlot(SkinBD.CD4, dims = 1:15)
SkinBD.CD4 <- FindNeighbors(SkinBD.CD4, dims = 1:10)
SkinBD.CD4 <- FindClusters(SkinBD.CD4, resolution = 0.8, verbose = FALSE)
SkinBD.CD4 <- RunUMAP(SkinBD.CD4, dims = 1:10)
DimPlot(SkinBD.CD4, label = TRUE, label.size = 15, group.by = "seurat_clusters", pt.size = 2.5)
DimPlot(SkinBD.CD4, label = FALSE, group.by = "Cond", cols = c("red", "black"), pt.size = 2.5)
levels(SkinBD.CD4)
levels(SkinBD.CD4)<- c("Burn_1", "Burn_2", "Burn_3", "NonBurn_1", "NonBurn_2", "NonBurn_3")
Idents(object = SkinBD.CD4) <- "State"
CD4.State.markers <- FindAllMarkers(SkinBD.CD4, group.by = "State" ,only.pos = TRUE, min.pct = 0.25, test.use ="MAST")
CD4.State.markers %>%
    group_by(cluster) %>%
    top_n(n = 8, wt = avg_log2FC) -> top8_CD4.State
DoHeatmap(SkinBD.CD4, features = top8O_CD4.State$gene, draw.lines = FALSE) + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))
VlnPlot(SkinBD.CD4, features = c("CD4","S1PR1","CD69","LGALS1","CCR7", "SELL","ITGA4", "SELPLG", "CXCR4"), stack = TRUE, group.by = "Cond", fill.by ="ident", cols = c("red", "black"))
VlnPlot(SkinBD.CD4, features = c("TBX21","RORA","RORC","GATA3","FOXP3","EOMES", "BCL6", "BCL11B", "ZBTB16"), stack = TRUE, group.by = "Cond", fill.by ="ident", cols = c("red", "black"))
FeaturePlot(SkinBD.CD4, features = c("TBX21","RORA","RORC","GATA3","FOXP3","EOMES", "BCL6", "BCL11B", "ZBTB16"), cols = plasma)
VlnPlot(SkinBD.CD4, features = c("IFNG","IL32","TNF","IL2","NAMPT","CSF3","TNFSF10","FASLG", "CSF2"), stack = TRUE, group.by = "Cond", fill.by ="ident", cols = c("red", "black"))
VlnPlot(SkinBD.CD4, features = c("PIK3IP1","CCR7"), stack = TRUE, group.by = "Cond", fill.by ="ident", cols = c("red", "black"))
DotPlot(SkinBD.CD4, features = c("CD4","S1PR1","CD69","LGALS1","CCR7", "SELL","ITGA4", "SELPLG", "CXCL4"), group.by = "Cond") + RotatedAxis()+ scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))
#CD4 AbSeq analysis
SkinBD.CD4 <- NormalizeData(object = SkinBD.CD4, assay = "AbSeq", normalization.method = 'CLR')
SkinBD.CD4 <- ScaleData(object = SkinBD.CD4, assay = 'AbSeq')
VlnPlot(SkinBD.CD4, features = c("CD152-CTLA4","CD279:EH12-1-PDCD1","CD25:2A3-IL2RA","CD161:DX12-KLRB1","CD69-CD69","CD38:HIT2-CD38"), assay = "AbSeq", stack = TRUE, group.by = "Cond", fill.by ="ident", cols = c("red", "black"))

#CD8 RNA analysis
SkinBD.CD8 <- subset(x = SkinBD, subset = CD8A > 2, idents = "CD8 T cells")
VizDimLoadings(SkinBD.CD8, dims = 1:2, reduction = "pca")
DimPlot(SkinBD.CD8, reduction = "pca")
DimHeatmap(SkinBD.CD8, dims = 1, cells = 500, balanced = TRUE)
SkinBD.CD8 <- JackStraw(SkinBD.CD8, num.replicate = 100, prop.freq = 0.05)
SkinBD.CD8 <- ScoreJackStraw(SkinBD.CD8, dims = 1:20)
JackStrawPlot(SkinBD.CD8, dims = 1:15)
SkinBD.CD8 <- FindNeighbors(SkinBD.CD8, dims = 1:6)
SkinBD.CD8 <- FindClusters(SkinBD.CD8, resolution = 0.8, verbose = FALSE)
SkinBD.CD8 <- RunUMAP(SkinBD.CD8, dims = 1:15)
DimPlot(SkinBD.CD8, label = TRUE, label.size = 15, group.by = "seurat_clusters", pt.size = 2.5)
DimPlot(SkinBD.CD8, label = FALSE, group.by = "Cond", cols = c("red", "black"), pt.size = 2.5)
levels(SkinBD.CD8)
levels(SkinBD.CD8)<- c("Burn_1", "Burn_2", "Burn_3", "NonBurn_1", "NonBurn_2", "NonBurn_3")
Idents(object = SkinBD.CD8) <- "State"
CD8.State.markers <- FindAllMarkers(SkinBD.CD8, group.by = "State" ,only.pos = TRUE, min.pct = 0.25, test.use ="MAST")
CD8.State.markers %>%
    group_by(cluster) %>%
    top_n(n = 8, wt = avg_log2FC) -> top8_CD8.State
DoHeatmap(SkinBD.CD8, features = top5_CD8.State$gene, draw.lines = FALSE) + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))
VlnPlot(SkinBD.CD8, features = c("CD8A","S1PR1","IL7R","CD69","LGALS1","CCR7", "SELL", "IL4R", "ITGA4"), stack = TRUE, group.by = "Cond", fill.by ="ident", cols = c("red", "black"))
VlnPlot(SkinBD.CD8, features = c("TBX21","RORA","RORC","GATA3","FOXP3","EOMES", "BCL6", "BCL11B", "ZBTB16"), stack = TRUE, group.by = "Cond", fill.by ="ident", cols = c("red", "black"))
FeaturePlot(SkinBD.CD8, features = c("TBX21","RORA","RORC","GATA3","FOXP3","EOMES", "BCL6", "BCL11B", "ZBTB16"), cols = plasma)
VlnPlot(SkinBD.CD8, features = c("IFNG","IL32","TNF","IL2","NAMPT","CSF3","TNFSF10","FASLG", "CSF2"), stack = TRUE, group.by = "Cond", fill.by ="ident", cols = c("red", "black"))
DotPlot(SkinBD.CD8, features = c("CD4","S1PR1","IL7R","CD69","LGALS1","CCR7", "SELL", "KLRB1", "IL4R", "ITGA4"), group.by = "Cond") + RotatedAxis()+ scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))

#CD8 Abseq analysis
SkinBD.CD8 <- NormalizeData(object = SkinBD.CD8, assay = "AbSeq", normalization.method = 'CLR')
SkinBD.CD8 <- ScaleData(object = SkinBD.CD8, assay = 'AbSeq')
VlnPlot(SkinBD.CD8, features = c("CD152-CTLA4","CD279:EH12-1-PDCD1","CD25:2A3-IL2RA","CD161:DX12-KLRB1","CD69-CD69","CD38:HIT2-CD38"), assay = "AbSeq", stack = TRUE, group.by = "Cond", fill.by ="ident", cols = c("red", "black"))


#TCRgd RNA analysis
VlnPlot(SkinBD, features = c("TRDC", "TARP-refseq", "CD300A"), ncol = 3, assay = "RNA")
SkinBD.gdT <- subset(x = SkinBD, subset = TRDC > 2 & `TARP-refseq` > 2)
VizDimLoadings(SkinBD.gdT, dims = 1:2, reduction = "pca")
DimPlot(SkinBD.gdT, reduction = "pca")
DimHeatmap(SkinBD.gdT, dims = 1, cells = 500, balanced = TRUE)
SkinBD.gdT <- JackStraw(SkinBD.gdT, num.replicate = 100, prop.freq = 0.05)
SkinBD.gdT <- ScoreJackStraw(SkinBD.gdT, dims = 1:20)
JackStrawPlot(SkinBD.gdT, dims = 1:15)
SkinBD.gdT <- FindNeighbors(SkinBD.gdT, dims = 1:8)
SkinBD.gdT <- FindClusters(SkinBD.gdT, resolution = 0.8, verbose = FALSE)
SkinBD.gdT <- RunUMAP(SkinBD.gdT, dims = 1:15)
DimPlot(SkinBD.gdT, label = FALSE, pt.size = 1)
DimPlot(SkinBD.gdT, label = FALSE, group.by = "Cond", cols = c("red", "black"), pt.size = 3)
gdT.markers <- FindAllMarkers(SkinBD.gdT, only.pos = TRUE, min.pct = 0.25)
gdT.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10_gdT
DoHeatmap(SkinBD.gdT, features = top10_gdT$gene, draw.lines = FALSE) + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))
VlnPlot(SkinBD.gdT, features = c("TRDC","S1PR1","CD69","LGALS1","CCR7", "SELL"), stack = TRUE, group.by = "Cond", fill.by ="ident", cols = c("red", "black"))
VlnPlot(SkinBD.gdT, features = c("TBX21","RORA","RORC","GATA3","FOXP3","EOMES", "BCL6", "BCL11B", "ZBTB16"), stack = TRUE, group.by = "Cond", fill.by ="ident", cols = c("red", "black"))
FeaturePlot(SkinBD.gdT, features = c("TBX21","RORA","RORC","GATA3","FOXP3","EOMES", "BCL6", "BCL11B", "ZBTB16"), cols = plasma)
VlnPlot(SkinBD.CD8, features = c("IFNG","IL32","TNF","IL2","NAMPT","CSF3","TNFSF10","FASLG", "CSF2"), stack = TRUE, group.by = "Cond", fill.by ="ident", cols = c("red", "black"))
DotPlot(SkinBD.CD8, features = c("CD4","S1PR1","IL7R","CD69","LGALS1","CCR7", "SELL", "KLRB1", "IL4R", "ITGA4"), group.by = "Cond") + RotatedAxis()+ scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))

#TCRgd AbSeq Analysis
SkinBD.gdT <- NormalizeData(object = SkinBD.gdT, assay = "AbSeq", normalization.method = 'CLR')
SkinBD.gdT <- ScaleData(object = SkinBD.gdT, assay = 'AbSeq')
VlnPlot(SkinBD.gdT, features = c("CD152-CTLA4","CD279:EH12-1-PDCD1","CD25:2A3-IL2RA","CD161:DX12-KLRB1","CD69-CD69","CD38:HIT2-CD38"), assay = "AbSeq", stack = TRUE, group.by = "Cond", fill.by ="ident", cols = c("red", "black"))


#CD3 AbSeq analysis
DefaultAssay(object = SkinBD) <- 'AbSeq'
SkinBD.AB <- SkinBD
SkinBD.AB <- NormalizeData(object = SkinBD.AB, assay = "AbSeq", normalization.method = 'CLR')
SkinBD.AB <- FindVariableFeatures(object = SkinDB.AB, verbose = FALSE)
SkinBD.AB@assays$AbSeq@var.features
SkinBD.AB <- ScaleData(object = SkinBD.AB, assay = 'AbSeq')
SkinBD.AB <- RunPCA(SkinBD.AB, features = VariableFeatures(object = SkinBD.AB))
ElbowPlot(object = SkinBD.AB, ndims = 5)
SkinBD.AB <- FindNeighbors(object = SkinBD.AB, reduction = 'pca', dims = 1:5)
SkinBD.AB <- FindClusters(object = SkinBD.AB, dims.use = 1:5, verbose = TRUE)
SkinBD.AB <- RunUMAP(SkinBD.AB, dims = 1:5)
DimPlot(SkinBD.AB, label = FALSE, pt.size = 1)
DimPlot(SkinBD.AB, label = FALSE, group.by = "Cond", cols = c("red", "black"), pt.size = 1)
VlnPlot(SkinBD.AB, features = c("CD152-CTLA4","CD279:EH12-1-PDCD1","CD25:2A3-IL2RA","CD161:DX12-KLRB1","CD69-CD69","CD366-Havcr2","CD38:HIT2-CD38","TCR-Valpha7-TRAV7"), stack = TRUE, group.by = "Cond", fill.by ="ident", cols = c("red", "black"))

#MAIT analysis
SkinBD <- NormalizeData(object = SkinBD, assay = "AbSeq", normalization.method = 'CLR')
SkinBD <- ScaleData(object = SkinBD, assay = 'AbSeq')
VlnPlot(SkinBD, features = c("CD161:DX12-KLRB1", "TCR-Valpha7-TRAV7"), ncol = 2, assay = "AbSeq")
SkinBD.MAIT <- subset(x = SkinBD, subset = `TCR-Valpha7-TRAV7` > 1.5 & `CD161:DX12-KLRB1` > 1)
DefaultAssay(SkinBD.MAIT) <- "RNA"
VizDimLoadings(SkinBD.MAIT, dims = 1:2, reduction = "pca")
DimPlot(SkinBD.MAIT, reduction = "pca")
DimHeatmap(SkinBD.MAIT, dims = 1, cells = 500, balanced = TRUE)
SkinBD.MAIT <- JackStraw(SkinBD.MAIT, num.replicate = 100, prop.freq = 0.08)
SkinBD.MAIT <- ScoreJackStraw(SkinBD.MAIT, dims = 1:15)
JackStrawPlot(SkinBD.MAIT, dims = 1:10)
SkinBD.MAIT <- FindNeighbors(SkinBD.MAIT, dims = 1:6)
SkinBD.MAIT <- FindClusters(SkinBD.MAIT, resolution = 0.5, verbose = FALSE)
SkinBD.MAIT <- RunUMAP(SkinBD.MAIT, dims = 1:15)
DimPlot(SkinBD.MAIT, label = FALSE, pt.size = 5)
DimPlot(SkinBD.MAIT, label = FALSE, group.by = "Cond", cols = c("red", "black"), pt.size = 3)
MAIT.markers <- FindAllMarkers(SkinBD.MAIT, only.pos = TRUE, min.pct = 0.25)
MAIT.markers %>%
    group_by(cluster) %>%
    top_n(n = 5, wt = avg_log2FC) -> top5_MAIT
DoHeatmap(SkinBD.MAIT, features = top5_MAIT$gene, draw.lines = FALSE) + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))
VlnPlot(SkinBD.MAIT, features = c("KLRB1","S1PR1","CD69","LGALS1","CCR7", "SELL"), stack = TRUE, group.by = "Cond", fill.by ="ident", cols = c("red", "black"))
VlnPlot(SkinBD.MAIT, features = c("TBX21","RORA","RORC","GATA3","FOXP3","EOMES", "BCL6", "BCL11B", "ZBTB16"), stack = TRUE, group.by = "Cond", fill.by ="ident", cols = c("red", "black"))
FeaturePlot(SkinBD.MAIT, features = c("TBX21","RORA","RORC","GATA3","FOXP3","EOMES", "BCL6", "BCL11B", "ZBTB16"), cols = plasma)
VlnPlot(SkinBD.MAIT, features = c("PRF1","IL32","CXCL8","ICOS","NAMPT","GZMB","GZMK","FASLG"), stack = TRUE, group.by = "Cond", fill.by ="ident", cols = c("red", "black"))
DotPlot(SkinBD.MAIT, features = c("CD4","S1PR1","IL7R","CD69","LGALS1","CCR7", "SELL", "KLRB1", "IL4R", "ITGA4"), group.by = "Cond") + RotatedAxis()+ scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))

#MAIT AbSeq Analysis
SkinBD.MAIT <- NormalizeData(object = SkinBD.MAIT, assay = "AbSeq", normalization.method = 'CLR')
SkinBD.MAIT <- ScaleData(object = SkinBD.MAIT, assay = 'AbSeq')
VlnPlot(SkinBD.MAIT, features = c("CD152-CTLA4","CD279:EH12-1-PDCD1","CD25:2A3-IL2RA","CD161:DX12-KLRB1","CD69-CD69","CD38:HIT2-CD38","TCR-Valpha7-TRAV7"), assay = "AbSeq", stack = TRUE, group.by = "Cond", fill.by ="ident", cols = c("red", "black"))

saveRDS(SkinBD, file = "~/Desktop/Leung_Lab/Seurat/Skin_BD/Skin_Seurat/SkinBD_all")

#DE analysis between Burn and Non-Burn mRNA

Idents(object = SkinBD.CD4) <- "Cond"
CD4.DE <- FindMarkers(SkinBD.CD4, ident.1 = "Burn", ident.2 = "Non_Burn", min.pct = 0.15)

Idents(object = SkinBD.CD8) <- "Cond"
CD8.DE <- FindMarkers(SkinBD.CD8, ident.1 = "Burn", ident.2 = "Non_Burn", min.pct = 0.15)

Idents(object = SkinBD.gdT) <- "Cond"
gdT.DE <- FindMarkers(SkinBD.gdT, ident.1 = "Burn", ident.2 = "Non_Burn", min.pct = 0.2)

Idents(object = SkinBD.MAIT) <- "Cond"
MAIT.DE <- FindMarkers(SkinBD.MAIT, ident.1 = "Burn", ident.2 = "Non_Burn", min.pct = 0.05)

#DE analysis between Burn and Non-Burn AbSeq

DefaultAssay(SkinBD.CD4) <- "AbSeq"
Idents(object = SkinBD.CD4) <- "Cond"
CD4Ab.DE <- FindMarkers(SkinBD.CD4, ident.1 = "Burn", ident.2 = "Non_Burn", min.pct = 0.05)

DefaultAssay(SkinBD.CD8) <- "AbSeq"
Idents(object = SkinBD.CD8) <- "Cond"
CD8Ab.DE <- FindMarkers(SkinBD.CD8, ident.1 = "Burn", ident.2 = "Non_Burn", min.pct = 0.05)

DefaultAssay(SkinBD.gdT) <- "AbSeq"
Idents(object = SkinBD.gdT) <- "Cond"
gdTAb.DE <- FindMarkers(SkinBD.gdT, ident.1 = "Burn", ident.2 = "Non_Burn", min.pct = 0.05)

DefaultAssay(SkinBD.MAIT) <- "AbSeq"
Idents(object = SkinBD.MAIT) <- "Cond"
MAITAb.DE <- FindMarkers(SkinBD.MAIT, ident.1 = "Burn", ident.2 = "Non_Burn", min.pct = 0.05)

DefaultAssay(SkinBD) <- "AbSeq"
Idents(object = SkinBD) <- "Cond"
SkinBDAb.DE <- FindMarkers(SkinBD, ident.1 = "Burn", ident.2 = "Non_Burn", min.pct = 0.05)

#GSEA analysis

library(devtools)
install_github('immunogenomics/presto')
library(presto)
CD4.genes <- wilcoxauc(SkinBD.CD4, 'Cond')
dplyr::count(CD4.genes, group)
library(msigdbr)
library(fgsea)
library(dplyr)
library(ggplot2)
msigdbr_show_species()
CD4.genes %>%
     dplyr::filter(group == "Burn") %>%
     arrange(desc(logFC), desc(auc)) %>%
     head(n = 10)

CD4.Burn.genes<- CD4.genes %>%
    dplyr::filter(group == "Burn") %>%
    arrange(desc(auc)) %>%
    dplyr::select(feature, auc)

CD4.Burn.ranks <- deframe(CD4.Burn.genes)
head(CD4.Burn.ranks)
CD4.Burn.fgseaRes <- fgsea(fgsea_sets, stats = CD4.Burn.ranks, nperm = 1000)
CD4.Burn.fgseaResTidy <- CD4.Burn.fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))
CD4.Burn.fgseaResTidy %>%
    dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>%
    arrange(pval) %>%
    head()
ggplot(CD4.Burn.fgseaResTidy %>% filter(pval < 0.008) %>% head(n= 40), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 7.5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
        title="Hallmark pathways NES from GSEA") +
  theme_minimal()
plotEnrichment(fgsea_sets[["GSE11057_NAIVE_VS_EFF_MEMORY_CD4_TCELL_DN"]],
  CD4.Burn.ranks) + labs(title="GSE11057_NAIVE_VS_EFF_MEMORY_CD4_TCELL_DN")


CD4.genes %>%
     dplyr::filter(group == "Non_Burn") %>%
     arrange(desc(logFC), desc(auc)) %>%
     head(n = 10)

CD4.NBurn.genes<- CD4.genes %>%
    dplyr::filter(group == "Non_Burn") %>%
    arrange(desc(auc)) %>%
    dplyr::select(feature, auc)

CD4.NBurn.ranks <- deframe(CD4.NBurn.genes)
head(CD4.NBurn.ranks)
CD4.NBurn.fgseaRes <- fgsea(fgsea_sets, stats = CD4.NBurn.ranks, nperm = 1000)
CD4.NBurn.fgseaResTidy <- CD4.NBurn.fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))
CD4.NBurn.fgseaResTidy %>%
    dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>%
    arrange(pval) %>%
    head()
ggplot(CD4.NBurn.fgseaResTidy %>% filter(pval < 0.01) %>% head(n= 40), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 7.5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
        title="Hallmark pathways NES from GSEA") +
  theme_minimal()
plotEnrichment(fgsea_sets[["GSE11057_NAIVE_VS_EFF_MEMORY_CD4_TCELL_DN"]],
                 CD4.NBurn.ranks) + labs(title="GSE11057_NAIVE_VS_EFF_MEMORY_CD4_TCELL_DN")
