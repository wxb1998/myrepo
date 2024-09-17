library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)
library(DoubletFinder)
library(future)
library(RColorBrewer)
library(SeuratData)
library(patchwork)
library(Rserve)
set.seed(2)
###读取数据
brain <- readRDS("/data/home/quj_lab/wangxuebao/01_results/03_mouse_sport/05_umap/brain/brain_umap-final.rds")
DefaultAssay(brain) <- "RNA"
brain <- NormalizeData(object = brain, normalization.method = "LogNormalize")
plan("multicore", workers = 28)
options(future.globals.maxSize = 100000 * 1024^2)#100000MB~=100G
dir.create('/data/home/quj_lab/wangxuebao/01_results/03_mouse_sport/06_cellmarker/brain')
setwd('/data/home/quj_lab/wangxuebao/01_results/03_mouse_sport/06_cellmarker/brain')







pdf("brain_universal.pdf", width=16,height=12)
################################################################################
# Pericyte
FeaturePlot(brain, raster=FALSE,features = c('Pdgfrb','Mcam','Rgs5'))
VlnPlot(brain, features = c('Pdgfrb','Mcam','Rgs5'), pt.size = 0)
# EC
FeaturePlot(brain, raster=FALSE,features = c("Plvap", "Egfl7", "Adgrl4", "Flt1", "Cdh5", "Tm4sf1", "Pecam1"))
VlnPlot(brain, features = c("Plvap", "Egfl7", "Adgrl4", "Flt1", "Cdh5", "Tm4sf1", "Pecam1"), pt.size = 0)
dev.off()
pdf("brain_2.pdf", width=16,height=12)
################################################################################
#OPC 
FeaturePlot(brain, raster=FALSE,features = c('Ppfibp1',"Opalin", "Qdpr","Pcdh15"))
VlnPlot(brain, features = c('Ppfibp1',"Opalin", "Qdpr","Pcdh15"), pt.size = 0)
# Macrophage  Microglia
FeaturePlot(brain, raster=FALSE,features = c('C1qa',' Mrc1','Apbb1ip',"Ctss","Itgam","Ptprc"))
VlnPlot(brain, features = c('C1qa',' Mrc1','Apbb1ip',"Ctss","Itgam","Ptprc"), pt.size = 0)
#Ependymal 
FeaturePlot(brain, raster=FALSE,features = c('Foxj1','Cfap54','Dnah6'))
VlnPlot(brain, features = c('Foxj1','Cfap54','Dnah6'), pt.size = 0)
#Oligodendrocyte
FeaturePlot(brain, raster=FALSE,features = c("Mobp", "Mog", "Mbp", "Plp1"))
VlnPlot(brain, features = c("Mobp", "Mog", "Mbp", "Plp1"), pt.size = 0)
###Meningeal
FeaturePlot(brain, raster=FALSE,features = c('Pdgfra',"Dcn" , "Col3a1" ,"Igf2"))
VlnPlot(brain, features = c('Pdgfra',"Dcn" , "Col3a1" ,"Igf2"), pt.size = 0)
# Astrocytes
FeaturePlot(brain, raster=FALSE,features = c('Sparc', 'Id1', 'Id3', 'Hes5', 'Aldh1l1', 'Aldoc', 'Hopx', 'Timp4', 'Ndrg2', 'Gdf10', 'Sox9', 'Slc1a3', 'Sparcl1'))
VlnPlot(brain, features = c('Sparc', 'Id1', 'Id3', 'Hes5', 'Aldh1l1', 'Aldoc', 'Hopx', 'Timp4', 'Ndrg2', 'Gdf10', 'Sox9', 'Slc1a3', 'Sparcl1'), pt.size = 0)
dev.off()

pdf("brain_5.pdf", width=16,height=12)
################################################################################
# Excitatory neurons
FeaturePlot(brain, raster=FALSE,features = c('Slc17a7', 'Nrgn', 'Camk2a', 'Satb2', 'Tbr1','Thy1'))
VlnPlot(brain, features = c('Slc17a7', 'Nrgn', 'Camk2a', 'Satb2', 'Tbr1','Thy1'), pt.size = 0)
# inhibitory neurons 
FeaturePlot(brain, raster=FALSE,features = c('Gad1', 'Gad2', 'Slc32a1'))
VlnPlot(brain, features = c('Gad1', 'Gad2', 'Slc32a1'), pt.size = 0)
# Interneuron
FeaturePlot(brain, raster=FALSE,features = c("Lgr5","Cntn4","Rbfox1",'Syt1','Snap25','Grin1','Dlx1','Dlx2','Lhx6'))
VlnPlot(brain, features = c("Lgr5","Cntn4","Rbfox1",'Syt1','Snap25','Grin1','Dlx1','Dlx2','Lhx6'), pt.size = 0)
# NPC 
FeaturePlot(brain, raster=FALSE,features = c('Hmgb2', 'Hmgn2', 'H2afz', 'Ccnd2', 'Pax6', 'Sfrp1', 'Vim'))
VlnPlot(brain, features = c('Hmgb2', 'Hmgn2', 'H2afz', 'Ccnd2', 'Pax6', 'Sfrp1', 'Vim'), pt.size = 0)
#VLMC 
FeaturePlot(brain, raster=FALSE,features = c('Slc6a13'))
VlnPlot(brain, features = c('Slc6a13'), pt.size = 0)
#NSC
FeaturePlot(brain, raster=FALSE,features = c('Gfap','Sox2','Prox1'))
VlnPlot(brain, features = c('Gfap','Sox2','Prox1'), pt.size = 0)
dev.off()
