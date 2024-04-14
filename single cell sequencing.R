
library(magrittr)
library(data.table)
library(Seurat)
library(SCP)
library(tidyverse)

    a <- fread("GSE138852_counts.csv.gz") 
    a <- a[!duplicated(a$V1)]
    a <-  tibble::column_to_rownames(a,var = "V1")
b=fread("GSE138852_covariates.csv.gz",header = T)
b <-  tibble::column_to_rownames(b,var = "V1")
  sce=CreateSeuratObject(a,meta.data = b)
 
qs::qsave(sce,file=paste0("sce.qs"))




sce=qs::qread(file=paste0("sce.qs"))



# merge -------------------------------------------------------------------

DIR="GSE138852"
dir.create(DIR)
qs::qsave(sce,file=paste0(DIR,"/GSE138852.qs"))

panc8_sub=qs::qread(file=paste0(DIR,"/GSE138852.qs"))

# LIGER
# panc8_sub=pancreas_sub
panc8_sub <- Integration_SCP(srtMerge = panc8_sub, batch = "oupSample.batchCond", integration_method = "LIGER")
# qs::qsave(panc8_sub,file=paste0(DIR,"/LIGER.qs"))

# Seurat
# panc8_sub=pancreas_sub
panc8_sub <- Integration_SCP(srtMerge = panc8_sub, batch = "oupSample.batchCond", integration_method = "Seurat")
# qs::qsave(panc8_sub,file=paste0(DIR,"/Seurat.qs"))

# Harmony
panc8_sub <- Integration_SCP(srtMerge = panc8_sub, batch = "oupSample.batchCond", integration_method = "Harmony",cluster_resolution = seq(0.1,1,0.1))
qs::qsave(panc8_sub,file=paste0(DIR,"/Harmony.qs"))


panc8_sub$LIGER_SNN_res.0.6
panc8_sub$Seurat_SNN_res.0.6
panc8_sub$Harmony_SNN_res.0.6

CellDimPlot(
  srt = panc8_sub, group.by = c("Harmony_SNN_res.0.6","Seurat_SNN_res.0.6","LIGER_SNN_res.0.6","oupSample.cellType"),
  reduction = "HarmonyUMAP2D",
  title = "Harmony", theme_use = "theme_blank"
)
ggsave(paste0(DIR,"/3.CellDimPlot_batch_HarmonyUMAP2D.pdf"),width = 20,height = 20)


CellDimPlot(
  srt = panc8_sub, group.by = c("Harmony_SNN_res.0.6","Seurat_SNN_res.0.6","LIGER_SNN_res.0.6","oupSample.cellType"),
  reduction = "LIGERUMAP2D",
  title = "LIGER", theme_use = "theme_blank"
)
ggsave(paste0(DIR,"/3.CellDimPlot_batch_LIGERUMAP2D.pdf"),width = 20,height = 20)

CellDimPlot(
  srt = panc8_sub, group.by = c("Harmony_SNN_res.0.6","Seurat_SNN_res.0.6","LIGER_SNN_res.0.6","oupSample.cellType"),
  reduction = "SeuratUMAP2D",
  title = "Seurat", theme_use = "theme_blank"
)
ggsave(paste0(DIR,"/3.CellDimPlot_batch_SeuratUMAP2D.pdf"),width = 20,height = 20)



hub=readxl::read_xlsx("hub基因.xlsx",col_names = F)

FeaturePlot(panc8_sub,features =hub$...1,reduction = "HarmonyUMAP2D" ,cols = c("white","red"),ncol = 5)
  ggsave(paste0("FeaturePlot_hub.pdf"),width = 20,height = 15)


CellDimPlot(
  srt = panc8_sub, group.by = c("oupSample.cellType"),split.by = "oupSample.subclustCond",
  reduction = "HarmonyUMAP2D",
  title = "Harmony", theme_use = "theme_blank"
)
ggsave(paste0("CellDimPlot_group_HarmonyUMAP2D.pdf"),width = 20,height = 10)





# sub ---------------------------------------------------------------------



panc8_sub=subset(panc8_sub,oupSample.subclustCond%in%c("AD","ct"))


CellDimPlot(
  srt = panc8_sub, group.by = c("Harmony_SNN_res.0.6","Seurat_SNN_res.0.6","LIGER_SNN_res.0.6","oupSample.cellType"),
  reduction = "HarmonyUMAP2D",
  title = "Harmony", theme_use = "theme_blank"
)
ggsave(paste0(DIR,"/3.CellDimPlot_batch_HarmonyUMAP2D_sub.pdf"),width = 20,height = 20)


CellDimPlot(
  srt = panc8_sub, group.by = c("Harmony_SNN_res.0.6","Seurat_SNN_res.0.6","LIGER_SNN_res.0.6","oupSample.cellType"),
  reduction = "LIGERUMAP2D",
  title = "LIGER", theme_use = "theme_blank"
)
ggsave(paste0(DIR,"/3.CellDimPlot_batch_LIGERUMAP2D_sub.pdf"),width = 20,height = 20)

CellDimPlot(
  srt = panc8_sub, group.by = c("Harmony_SNN_res.0.6","Seurat_SNN_res.0.6","LIGER_SNN_res.0.6","oupSample.cellType"),
  reduction = "SeuratUMAP2D",
  title = "Seurat", theme_use = "theme_blank"
)
ggsave(paste0(DIR,"/3.CellDimPlot_batch_SeuratUMAP2D_sub.pdf"),width = 20,height = 20)



hub=readxl::read_xlsx("hub基因.xlsx",col_names = F)

panc8_sub_ad=subset(panc8_sub,oupSample.subclustCond=="AD")
FeaturePlot(panc8_sub_ad,features =hub$...1,reduction = "HarmonyUMAP2D" ,cols = c("white","red"),ncol = 5)
ggsave(paste0("FeaturePlot_hub_AD.pdf"),width = 20,height = 15)


panc8_sub_ct=subset(panc8_sub,oupSample.subclustCond=="ct")
FeaturePlot(panc8_sub_ct,features =hub$...1,reduction = "HarmonyUMAP2D" ,cols = c("white","red"),ncol = 5)
ggsave(paste0("FeaturePlot_hub_ct.pdf"),width = 20,height = 15)

CellDimPlot(
  srt = panc8_sub, group.by = c("oupSample.cellType"),split.by = "oupSample.subclustCond",
  reduction = "HarmonyUMAP2D",
  title = "Harmony", theme_use = "theme_blank"
)
ggsave(paste0("CellDimPlot_group_HarmonyUMAP2D_sub.pdf"),width = 20,height = 10)

Idents(panc8_sub)=panc8_sub$oupSample.cellType
CellChat::StackedVlnPlot(panc8_sub,features = intersect(hub$...1,rownames(panc8_sub)),split.by = "oupSample.subclustCond")

ggsave(paste0("StackedVlnPlot.pdf"),width = 10,height = 10)
