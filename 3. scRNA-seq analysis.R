library(Seurat)
# BiocManager::install("UCell")
library(UCell)
library(dplyr)
library(ggplot2)
library(stringr)
library(viridisLite)
load("bfreaname_seurat_obj.Rdata")
up_genes<-read.table("up_regulated_genes.txt",sep = "\t",header = T)[,1]
down_genes<-read.table("down_regulated_genes.txt",sep = "\t",header = T)[,1]
up_down_genes<-read.table("updown_regulated_genes.txt",sep = "\t",header = T)[,1]
down_up_genes<-read.table("downup_regulated_genes.txt",sep = "\t",header = T)[,1]
markers <- list()
markers$Consis_up <- up_genes
markers$Consis_down <- down_genes
markers$Up_then_down <- up_down_genes
markers$Down_then_up <- down_up_genes

BPPARAM <- BiocParallel::MulticoreParam(workers=14)


marker_score <- AddModuleScore_UCell(bfreaname_seurat_obj,
                                     features=markers)
a <- colnames(marker_score@meta.data) %>% str_subset("_UCell")
model_gene_score <- FeaturePlot(marker_score,features = a,order = T, ncol = 2, cols = viridis(256))
ggsave(filename = "model_gene_score.pdf",plot = model_gene_score,width = 10,height = 8)
dim_result<-DimPlot(bfreaname_seurat_obj, reduction = "umap", label = TRUE, pt.size = 0.5, cols = viridis(17),label.color = "white",label.box = T)
ggsave(filename = "celltype_DimResults.pdf",plot = dim_result,width = 10,height = 8)

mesenchymal_stem_cell <-c("CD3D","CD3E","LEF1","TCF7","IL7R","CD40LG",####CD4+ T细胞---0
                          "CD160","CD247","GNLY","NKG7",###NK细胞---1
                          "CD14","S100A12","CCL3",####单核细胞(CCL3)---2
                          "CD8A","CD8B","GZMK",####CD8+ T细胞---3
                          "TCF7L2","CEBPB","HES4",###MONO---4
                          "PF4","PPBP",####巨核祖细胞（巨核细胞）---5
                          "CD79A","CD79B","MS4A1",'CD19',####B细胞---6
                          "CD1C","FCER1A","ENHO",####树突状细胞---7
                          "MZB1","FKBP11","PRDM1","IGHG1",####B细胞-MZB1---8
                          "LILRA4","CLEC4C","IL3RA",####树突状细胞---9
                          "IFI6","IFIT3",####单核细胞(IFI6-IFIT3)---10
                          "PF4","PPBP",####巨核祖细胞（巨核细胞）---11
                          "EAF1","SLC16A3",####单核细胞(EAF1-SLC16A3)---12
                          "SPINK2","PRSS57","CD34","ALDH1A1",###造血干细胞---13
                          "ALAS2","HBB","AHSP","HBD","HBM",####红细胞---14
                          "CLEC9A","BATF3","WDFY4",####树突状细胞(CLEC9A)---15
                          "CD79A","CD79B","MS4A1",'CD19')####B细胞---16
marker_cell <- DotPlot(marker_score, features = unique(mesenchymal_stem_cell),group.by = "anno",cols = viridis(2),col.min = 0,col.max = 2)+RotatedAxis()
ggsave(filename = "marker_cell.pdf",plot = marker_cell ,width = 14,height = 8)
bfreaname_seurat_obj[["percent.mt"]] <- PercentageFeatureSet(bfreaname_seurat_obj, pattern = "^MT-")
countvol_plot <- VlnPlot(bfreaname_seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(filename = "QC.pdf",plot = countvol_plot,width = 12,height = 8)
