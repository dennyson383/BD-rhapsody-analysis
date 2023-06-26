#### sever script
library(EnsDb.Mmusculus.v79)
library(Seurat)
library(dplyr)
library(tidyverse)
library(SeuratWrappers)
library(Matrix)
library(collapse)
library(data.table)
library(cowplot)
library(multtest)
library(metap)

#Aux function
tomara_q_de = function(txi3.cartridge.1, project, nCount = 30, cart1 = F) {
  #Cartridge1 was already exported as $counts, therefore this line is not needed
  if(!cart1) {txi3.cartridge.1 <- txi3.cartridge.1$counts}
  
  sumsofcols = colSums(txi3.cartridge.1)
  sumofrows = rowSums(txi3.cartridge.1)
  
  txi3.cartridge.1 = txi3.cartridge.1[sumofrows > 0 , sumsofcols > nCount]
  transcript = as.data.frame(txi3.cartridge.1@Dimnames[[1]])
  colnames(transcript)[1] = "tx_id"
  
  txmap = left_join(transcript, tx2gene, copy = T)
  
  tx2 = as.matrix(txmap)
  rownames(tx2) = txmap$tx_id
  #tail(tx2)
  
  txi3.cartridge.1 <- txi3.cartridge.1[!is.na(txmap$gene_name),]
  
  breaks_txi_cartridge.1 = ceiling(nrow(txi3.cartridge.1)/8000)
  
  list_txi_cartridge.1 = lapply(1:breaks_txi_cartridge.1, function(x){
    
    begin = (x-1)*8000 + 1
    end = min(x*8000, nrow(txi3.cartridge.1))
    print(begin)
    
    breaks_txi = txi3.cartridge.1[begin:end, ]
    txi_merged = merge(breaks_txi, tx2, by = "row.names")
    
    lastcol = ncol(txi_merged)
    txi_merged = txi_merged[-(lastcol-1)]
    txi_merged.1 = fsum(get_vars(txi_merged, is.numeric), g = txi_merged[lastcol-1])
    
    return(list(test.1 = data.frame(txi_merged.1)))
    
  })
  
  sparse_catridge.1 = do.call(rbind, lapply(list_txi_cartridge.1, "[[", 1))
  
  n0 <- which(rownames(sparse_catridge.1) == "")
  n1 <- which(rownames(sparse_catridge.1) == "1")
  
  sparse_catridge.1 = sparse_catridge.1[-c(n0, n1),]
  
  sparse_catridge.1.1 = fsum(get_vars(sparse_catridge.1, is.numeric), g = rownames(sparse_catridge.1))
  
  sparse_catridge.1.1 = as.matrix(sparse_catridge.1.1)
  
  txi_sparse_catridge.1 = as(sparse_catridge.1.1, "dgCMatrix")
  
  cartridge.data = CreateSeuratObject(txi_sparse_catridge.1, min.cells = 0, min.features = 0, project = project)
  
  return(cartridge.data)  
}
stag_filter = function(sObj, sampletag).                                  {
  
  if (sampletag == 1){
    sample_obj <- subset(sObj, subset = 
                           sampletag1 > 0 & 
                           sampletag2 < 1 & 
                           sampletag3 < 1 & 
                           sampletag4 < 1 &
                           sampletag5 < 1 &
                           sampletag6 < 1 & 
                           sampletag7 < 1 & 
                           sampletag8 < 1 & 
                           sampletag9 < 1 &
                           sampletag10 < 1 &
                           sampletag11 < 1 &
                           sampletag12 < 1)
    
    return(sample_obj) 
  }
  
  if (sampletag == 2){
    sample_obj <- subset(sObj, subset = 
                           sampletag1 < 1 & 
                           sampletag2 > 0 & 
                           sampletag3 < 1 & 
                           sampletag4 < 1 &
                           sampletag5 < 1 &
                           sampletag6 < 1 & 
                           sampletag7 < 1 & 
                           sampletag8 < 1 & 
                           sampletag9 < 1 &
                           sampletag10 < 1 &
                           sampletag11 < 1 &
                           sampletag12 < 1)
    
    return(sample_obj)
  }
  
  if (sampletag == 3){
    sample_obj <- subset(sObj, subset = 
                           sampletag1 < 1 & 
                           sampletag2 < 1 & 
                           sampletag3 > 0 & 
                           sampletag4 < 1 &
                           sampletag5 < 1 &
                           sampletag6 < 1 & 
                           sampletag7 < 1 & 
                           sampletag8 < 1 & 
                           sampletag9 < 1 &
                           sampletag10 < 1 &
                           sampletag11 < 1 &
                           sampletag12 < 1)
    return(sample_obj)
  }
  
  if (sampletag == 4){
    sample_obj <- subset(sObj, subset = 
                           sampletag1 < 1 & 
                           sampletag2 < 1 & 
                           sampletag3 < 1 & 
                           sampletag4 > 0 &
                           sampletag5 < 1 &
                           sampletag6 < 1 & 
                           sampletag7 < 1 & 
                           sampletag8 < 1 & 
                           sampletag9 < 1 &
                           sampletag10 < 1 &
                           sampletag11 < 1 &
                           sampletag12 < 1)
    return(sample_obj)
  }
  
  if (sampletag == 5){
    sample_obj <- subset(sObj, subset = 
                           sampletag1 < 1 & 
                           sampletag2 < 1 & 
                           sampletag3 < 1 & 
                           sampletag4 < 1 &
                           sampletag5 > 0 &
                           sampletag6 < 1 & 
                           sampletag7 < 1 & 
                           sampletag8 < 1 & 
                           sampletag9 < 1 &
                           sampletag10 < 1 &
                           sampletag11 < 1 &
                           sampletag12 < 1)
    return(sample_obj)
  }
  
  if (sampletag == 6){
    sample_obj <- subset(sObj, subset = 
                           sampletag1 < 1 & 
                           sampletag2 < 1 & 
                           sampletag3 < 1 & 
                           sampletag4 < 1 &
                           sampletag5 < 1 &
                           sampletag6 > 0 & 
                           sampletag7 < 1 & 
                           sampletag8 < 1 & 
                           sampletag9 < 1 &
                           sampletag10 < 1 &
                           sampletag11 < 1 &
                           sampletag12 < 1)
    return(sample_obj)
  }
  
  if (sampletag == 7){
    sample_obj <- subset(sObj, subset = 
                           sampletag1 < 1 & 
                           sampletag2 < 1 & 
                           sampletag3 < 1 & 
                           sampletag4 < 1 &
                           sampletag5 < 1 &
                           sampletag6 < 1 & 
                           sampletag7 > 0 & 
                           sampletag8 < 1 & 
                           sampletag9 < 1 &
                           sampletag10 < 1 &
                           sampletag11 < 1 &
                           sampletag12 < 1)
    return(sample_obj)
  }
  
  if (sampletag == 8){
    sample_obj <- subset(sObj, subset = 
                           sampletag1 < 1 & 
                           sampletag2 < 1 & 
                           sampletag3 < 1 & 
                           sampletag4 < 1 &
                           sampletag5 < 1 &
                           sampletag6 < 1 & 
                           sampletag7 < 1 & 
                           sampletag8 > 0 & 
                           sampletag9 < 1 &
                           sampletag10 < 1 &
                           sampletag11 < 1 &
                           sampletag12 < 1)
    return(sample_obj)
  }
  
  if (sampletag == 9){
    sample_obj <- subset(sObj, subset = 
                           sampletag1 < 1 & 
                           sampletag2 < 1 & 
                           sampletag3 < 1 & 
                           sampletag4 < 1 &
                           sampletag5 < 1 &
                           sampletag6 < 1 & 
                           sampletag7 < 1 & 
                           sampletag8 < 1 & 
                           sampletag9 > 0 &
                           sampletag10 < 1 &
                           sampletag11 < 1 &
                           sampletag12 < 1)
    return(sample_obj)
  }
  
  if (sampletag == 10){
    sample_obj <- subset(sObj, subset = 
                           sampletag1 < 1 & 
                           sampletag2 < 1 & 
                           sampletag3 < 1 & 
                           sampletag4 < 1 &
                           sampletag5 < 1 &
                           sampletag6 < 1 & 
                           sampletag7 < 1 & 
                           sampletag8 < 1 & 
                           sampletag9 < 1 &
                           sampletag10 > 0 &
                           sampletag11 < 1 &
                           sampletag12 < 1)
    return(sample_obj)
  }
  
  if (sampletag == 11){
    sample_obj <- subset(sObj, subset = 
                           sampletag1 < 1 & 
                           sampletag2 < 1 & 
                           sampletag3 < 1 & 
                           sampletag4 < 1 &
                           sampletag5 < 1 &
                           sampletag6 < 1 & 
                           sampletag7 < 1 & 
                           sampletag8 < 1 & 
                           sampletag9 < 1 &
                           sampletag10 < 1 &
                           sampletag11 > 0 &
                           sampletag12 < 1)
    return(sample_obj)
  }
  
  if (sampletag == 12){
    sample_obj <- subset(sObj, subset = 
                           sampletag1 < 1 & 
                           sampletag2 < 1 & 
                           sampletag3 < 1 & 
                           sampletag4 < 1 &
                           sampletag5 < 1 &
                           sampletag6 < 1 & 
                           sampletag7 < 1 & 
                           sampletag8 < 1 & 
                           sampletag9 < 1 &
                           sampletag10 < 1 &
                           sampletag11 < 1 &
                           sampletag12 > 0)
    return(sample_obj)
  }
}
saveScatter = function(sObj,
                       path,
                       mt         = 20  ,
                       nfeat.min  = 100 ,
                       nfeat.max  = 1000,
                       ncount.min = 100 ,
                       ncount.max = 2000,
                       cols       = NULL)                                 {
  
  VlnPlot(sObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  
  plot1 = FeatureScatter(sObj, feature1 = "nCount_RNA", feature2 = "percent.mt", cols = cols)   +
    geom_hline(yintercept = mt, linewidth = 0.3) +
    geom_vline(xintercept = ncount.min, linewidth = 0.3)
  
  plot2 = FeatureScatter(sObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = cols) +
    geom_hline(yintercept = c(nfeat.min, nfeat.max), linewidth = 0.3) +
    geom_vline(xintercept = c(ncount.min, ncount.max), linewidth = 0.3)
  
  plot1 + plot2
  
  ggsave(path, width = 30, height = 20, dpi = 400, units = "cm")

}
SingelCellQC = function(sObj, name, quality_before = T, quality_after = T){
  
  sObj[["percent.mt"]] = PercentageFeatureSet(sObj, pattern = "^mt-")
  print(dim(sObj))
  
  dir = paste0("Output/", name, "_before.png")
  quality_before = saveScatter(sObj, path = dir)
  
  subset_sObj = subset(sObj, subset =
                      nFeature_RNA > 100  & 
                      nFeature_RNA < 1000 &
                      percent.mt   < 20   &
                      nCount_RNA   > 100  &
                      nCount_RNA   < 2000 )
  print(dim(subset_sObj))
  dir.a = paste0("Output/", name, "_after.png")
  quality_after = saveScatter(subset_sObj, dir.a)
  
  return(subset_sObj)
}
NormalizeSeuratObj = function(sObj)                                       {
  
  vec = c("voce esta fazendo besteira", "tem certeza que voce está fazendo certo ?", "se eu fosse voce, parava!",
          "roda novamente isso, bicho", "tu é fresco é ?", "PAYSANDUUUU", "REMOOOO", "o paraná não é no sul",
          "o denny é demais", "igor bundão", "desde quando voce é inteligente?", "Senta que lavem história")
  
  i = sample(vec, 1)
  
  print(vec[round(runif(1, 1, length(vec)))])
  
  counts = GetAssayData(sObj, assay = "RNA")
  counts = counts[-(which(rownames(counts)  %in%  c('sampletag1','sampletag2','sampletag3','sampletag4',
                                                    'sampletag5','sampletag6','sampletag7','sampletag8',
                                                    'sampletag9','sampletag10','sampletag11','sampletag12'))),]
  sObj = subset(sObj, features = rownames(counts))
  
  return(NormalizeData(sObj, normalization.method = "LogNormalize", scale.factor = 1000))
}
######
#Preparing tx2gn.                                                          ----
txdb = EnsDb.Mmusculus.v79
tx2gene = transcripts(txdb, columns = c(listColumns(txdb, "tx_id"), "tx_id", "gene_name"), return.type = "DataFrame")

sample_tag = data.frame(tx_id     = c("sampletag1",
                                      "sampletag2",
                                      "sampletag3",
                                      "sampletag4",
                                      "sampletag5",
                                      "sampletag6",
                                      "sampletag7",
                                      "sampletag8",
                                      "sampletag9",
                                      "sampletag10",
                                      "sampletag11",
                                      "sampletag12"),
                        gene_name = c("sampletag1",
                                      "sampletag2",
                                      "sampletag3",
                                      "sampletag4",
                                      "sampletag5",
                                      "sampletag6",
                                      "sampletag7",
                                      "sampletag8",
                                      "sampletag9",
                                      "sampletag10",
                                      "sampletag11",
                                      "sampletag12"))

tx2gene = rbind(tx2gene, sample_tag)


#Reading txi.                                                              ----
obj1 = readRDS("~/Desktop/cart1.RDS")
obj2 = readRDS("~/Desktop/cart2.RDS")
obj3 = readRDS("~/Desktop/cart3.RDS")
obj4 = readRDS("~/Desktop/cart4.RDS")

#Quality control naive                                                     ----
# naive to multiple cartridge.
naive_s1  = stag_filter(obj1, 1)
naive_s6  = stag_filter(obj2, 6)
naive_s11 = stag_filter(obj3, 11)
naive_s4  = stag_filter(obj4, 4)

#naive s1
naive_s1  = SingelCellQC(naive_s1, name = "/naive/naive_s1")

#naive s6
naive_s6  = SingelCellQC(naive_s6, name = "/naive/naive_s6")

#naive s11
naive_s11 = SingelCellQC(naive_s11, name = "/naive/naive_s11_before")

#naive s4
naive_s4  = SingelCellQC(naive_s4, name = "/naive/naive_s4_before")

#Quality control pbs                                                       ----
#pbs to multiple cartridge
pbs_s2   = stag_filter(obj1, 2)
pbs_s7   = stag_filter(obj2, 7)
pbs_s12  = stag_filter(obj3, 12)
pbs_s5   = stag_filter(obj4, 5)

#pbs s2
pbs_s2  = SingelCellQC(pbs_s2, name = "/pbs/pbs_s2")

#pbs s7
pbs_s7  = SingelCellQC(pbs_s7, name = "/pbs/pbs_s7")

#pbs s12
pbs_s12  = SingelCellQC(pbs_s12, name = "/pbs/pbs_s12")

#pbs s5
pbs_s5  = SingelCellQC(pbs_s5, name = "/pbs/pbs_s5")

#Quality contron
#Quality control igg.                                                      ----
#igg to multiple cartridge
igg_s3  = stag_filter(obj1, 3)
igg_s8  = stag_filter(obj2, 8)
igg_s1  = stag_filter(obj3, 1)
igg_s6  = stag_filter(obj4, 6)

#igg s3
igg_s3  = SingelCellQC(igg_s3, name = "/igg/igg_s3")

#igg s8
igg_s8  = SingelCellQC(igg_s8, name = "/igg/igg_s8")

#igg s1
igg_s1  = SingelCellQC(igg_s1, name = "/igg/igg_s1")

#igg s6
igg_s6  = SingelCellQC(igg_s6, name = "/igg/igg_s6")

#Quality control antipd1.                                                  ----
#antipd1. to multiple cartridge
antipd1_s4 = stag_filter(obj1, 4)
antipd1_s9 = stag_filter(obj2, 9)
antipd1_s2 = stag_filter(obj3, 2)
antipd1_s7 = stag_filter(obj4, 7)

#antipd1 s4
antipd1_s4  = SingelCellQC(antipd1_s4, name = "/antipd1/antipd1_s4")

#antipd1 s9
antipd1_s9  = SingelCellQC(antipd1_s9, name = "/antipd1/antipd1_s9")

#antipd1 s2
antipd1_s2  = SingelCellQC(antipd1_s2, name = "/antipd1/antipd1_s2")

#antipd1 s7
antipd1_s7  = SingelCellQC(antipd1_s7, name = "/antipd1/antipd1_s7")

#Quality control antictla4                                                 ----
#antictla4 to multiple cartridge
antictla4_s5 = stag_filter(obj1, 5)
antictla4_s10 = stag_filter(obj2, 10)
antictla4_s3 = stag_filter(obj3, 3)
antictla4_s8 = stag_filter(obj4, 8)

#antictla4 s5
antictla4_s5  = SingelCellQC(antictla4_s5, name = "/antictla4/antictla4_s5")

#antictla4 s10
antictla4_s10  = SingelCellQC(antictla4_s10, name = "/antictla4/antictla4_s10")

#antictla4 s3
antictla4_s3  = SingelCellQC(antictla4_s3, name = "/antictla4/antictla4_s3")

#antictla4 s8
antictla4_s8  = SingelCellQC(antictla4_s8, name = "/antictla4/antictla4_s8")


#Normalize                                                                 ----
naive_s1_norm  = NormalizeSeuratObj(naive_s1)
naive_s6_norm  = NormalizeSeuratObj(naive_s6)
naive_s11_norm = NormalizeSeuratObj(naive_s11)
naive_s4_norm  = NormalizeSeuratObj(naive_s4)

pbs_s2_norm    = NormalizeSeuratObj(pbs_s2)
pbs_s7_norm    = NormalizeSeuratObj(pbs_s7)
pbs_s12_norm   = NormalizeSeuratObj(pbs_s12)
pbs_s5_norm    = NormalizeSeuratObj(pbs_s5)

igg_s3_norm    = NormalizeSeuratObj(igg_s3)
igg_s8_norm    = NormalizeSeuratObj(igg_s8)
igg_s1_norm    = NormalizeSeuratObj(igg_s1)
igg_s6_norm    = NormalizeSeuratObj(igg_s6)

antipd1_s4_norm  = NormalizeSeuratObj(antipd1_s4)
antipd1_s9_norm  = NormalizeSeuratObj(antipd1_s9)
antipd1_s2_norm  = NormalizeSeuratObj(antipd1_s2)
antipd1_s7_norm  = NormalizeSeuratObj(antipd1_s7)

antictla4_s5_norm   = NormalizeSeuratObj(antictla4_s5)
antictla4_s10_norm  = NormalizeSeuratObj(antictla4_s10)
antictla4_s3_norm   = NormalizeSeuratObj(antictla4_s3)
antictla4_s8_norm   = NormalizeSeuratObj(antictla4_s8)

#Merge group                                                               ----
naive_norm     = merge(naive_s1_norm,     y = c(naive_s6_norm, naive_s11_norm, naive_s4_norm), add.cell.ids = c("naive1", "naive6", "naive11", "naive4"))
pbs_norm       = merge(pbs_s2_norm,       y = c(pbs_s7_norm, pbs_s12_norm, pbs_s5_norm), add.cell.ids = c("pbs2", "pbs7", "pbs12", "pbs5"))
igg_norm       = merge(igg_s3_norm,       y = c(igg_s8_norm, igg_s1_norm, igg_s6_norm), add.cell.ids = c("igg3", "igg8", "igg1", "igg6"))
antipd1_norm   = merge(antipd1_s4_norm,   y = c(antipd1_s9_norm, antipd1_s2_norm, antipd1_s7_norm), add.cell.ids = c("antipd1.4", "antipd1.9", "antipd1.2", "antipd1.7"))
antictla4_norm = merge(antictla4_s5_norm, y = c(antictla4_s10_norm, antictla4_s3_norm, antictla4_s8_norm), add.cell.ids = c("antictla4.5", "antictla4.10", "antictla4.3", "antictla4.8"))

naive_norm$groups     = "naive"
pbs_norm$groups       = "pbs"
igg_norm$groups       = "igg"
antipd1_norm$groups   = "antipd1"
antictla4_norm$groups = "antictla4"

head(colnames(naive_norm))
head(colnames(pbs_norm))
head(colnames(igg_norm))
head(colnames(antipd1_norm))
head(colnames(antictla4_norm))

tail(colnames(naive_norm))
tail(colnames(pbs_norm))
tail(colnames(igg_norm))
tail(colnames(antipd1_norm))
tail(colnames(antictla4_norm))

#Create metadata                                                           ----
# A = naive cartridge 1
# B = naive cartridge 2
# C = naive cartridge 3
# D = naive cartridge 4
meta_naive        = LETTERS[Idents(object = naive_norm)]
names(meta_naive) = colnames(x = naive_norm)
naive_norm        = AddMetaData(object = naive_norm, metadata = meta_naive, col.name = 'letter.idents')
head(x = naive[[]]) #unique(sapply(X = strsplit(colnames(naive), split = "_"), FUN = "[", 1))
table(naive$orig.ident)
table(naive$letter.idents)


######
#Integrate variables                                                       ----
naive_pbs_norm       = FindIntegrationAnchors(object.list = list(naive_norm, pbs_norm, antipd1_norm), dims = 1:20)
naive_antipd1_norm   = FindIntegrationAnchors(object.list = list(naive_norm, antipd1_norm), dims = 1:20)
naive_antictla4_norm = FindIntegrationAnchors(object.list = list(naive_norm, antictla4_norm), dims = 1:20)

naive_pbs_norm     = IntegrateData(anchorset = naive_pbs_norm, dims = 1:20)
naive_antipd1_norm = IntegrateData(anchorset = naive_antipd1_norm, dims = 1:20)

DefaultAssay(naive_pbs_norm) = "integrated"

#naive vs pbs
# Run the standard workflow for visualization and clustering
naive_pbs_norm = ScaleData(naive_pbs_norm, verbose = FALSE)
naive_pbs_norm = RunPCA(naive_pbs_norm, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
naive_pbs_norm = RunUMAP(naive_pbs_norm, reduction = "pca", dims = 1:20)
naive_pbs_norm = FindNeighbors(naive_pbs_norm, reduction = "pca", dims = 1:20)
naive_pbs_norm = FindClusters(naive_pbs_norm, resolution = 0.5)

# Visualization
p1 = DimPlot(naive_pbs_norm, reduction = "umap", group.by = "groups")
p2 = DimPlot(naive_pbs_norm, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

DefaultAssay(naive_pbs_norm) = "RNA"
naive_pbs_norm.markers = FindConservedMarkers(naive_pbs_norm, ident.1 = 7, grouping.var = "groups", verbose = FALSE)
head(naive_pbs_norm.markers)
FeaturePlot(naive_pbs_norm, features = c("Cd4", "Cd8", "Cd25"), min.cutoff = "q9", split.by = "groups")

#naive vs antipd1
# Run the standard workflow for visualization and clustering
naive_antipd1_norm = ScaleData(naive_antipd1_norm, verbose = FALSE)
naive_antipd1_norm = RunPCA(naive_antipd1_norm, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
naive_antipd1_norm = RunUMAP(naive_antipd1_norm, reduction = "pca", dims = 1:20)
naive_antipd1_norm = FindNeighbors(naive_antipd1_norm, reduction = "pca", dims = 1:20)
naive_antipd1_norm = FindClusters(naive_antipd1_norm, resolution = 0.5)

# Visualization
p1 = DimPlot(naive_antipd1_norm, reduction = "umap", group.by = "groups")
p2 = DimPlot(naive_antipd1_norm, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

DefaultAssay(naive_antipd1_norm) = "RNA"
naive_antipd1_norm.markers = FindConservedMarkers(naive_antipd1_norm, ident.1 = 7, grouping.var = "groups", verbose = FALSE)
head(naive_antipd1_norm.markers)
FeaturePlot(naive_antipd1_norm, features = c("Cd3d", "Sell", "Crem", "Cd8a",  "Gnly", "Cd79", "Fcgr3a","Ccl2", "Ppbp",
                                             "foxp3"), min.cutoff = "q5")

##    

######
#PCA to every conditions                                                   ----
## naive
all.genes.naive = rownames(naive_norm)
naive_norm      = FindVariableFeatures(naive_norm, selection.method = "vst", nfeatures = 2000)
naive_norm      = ScaleData(naive_norm, features = all.genes.naive)
naive_norm      = RunPCA(naive_norm, features = VariableFeatures(naive_norm))

## pbs
all.genes.pbs   = rownames(pbs_norm)
pbs_norm        = FindVariableFeatures(pbs_norm, selection.method = "vst", nfeatures = 2000)
pbs_norm        = ScaleData(pbs_norm, features = all.genes.pbs)
pbs_norm        = RunPCA(pbs_norm, features = VariableFeatures(pbs_norm))

##igg
all.genes.igg   = rownames(igg_norm)
igg_norm        = FindVariableFeatures(igg_norm, selection.method = "vst", nfeatures = 2000)
igg_norm        = ScaleData(igg_norm, features = all.genes.igg)
igg_norm        = RunPCA(igg_norm, features = VariableFeatures(igg_norm))

##antipd1
all.genes.antipd1    = rownames(antipd1_norm)
antipd1_norm         = FindVariableFeatures(antipd1_norm, selection.method = "vst", nfeatures = 2000)
antipd1_norm         = ScaleData(antipd1_norm, features = all.genes.antipd1)
antipd1_norm         = RunPCA(antipd1_norm, features = VariableFeatures(antipd1_norm))

##antipd1
all.genes.antictla4   = rownames(antictla4_norm)
antictla4_norm        = FindVariableFeatures(antictla4_norm, selection.method = "vst", nfeatures = 2000)
antictla4_norm        = ScaleData(antictla4_norm, features = all.genes.antictla4)
antictla4_norm        = RunPCA(antictla4_norm, features = VariableFeatures(antictla4_norm))

######
#Dimplots plots.                                                           -----
head(naive_norm@meta.data)
DimPlot(naive_norm, reduction = "pca", split.by = rownames(naive_norm@meta.data))
DimPlot(pbs_norm, reduction = "pca", split.by = "orig.ident")
DimPlot(igg_norm, reduction = "pca", split.by = "orig.ident")
DimPlot(antipd1_norm, reduction = "pca", split.by = "orig.ident")
DimPlot(antictla4_norm, reduction = "pca", split.by = "orig.ident")
# run the above code to every obj
# FeaturePlot(naive_norm, features = c("nCount_RNA"))
# DimHeatmap(naive_norm, dims = 1, cells = 500, balanced = T)
# DimHeatmap(naive_norm, dims = 1:5, cells = 500, balanced = T)

######
#Determine the ‘dimensionality’ of the dataset                             ----
naive_norm = JackStraw(naive_norm, num.replicate = 100)
naive_norm = ScoreJackStraw(naive_norm, dims = 1:20)
JackStrawPlot(naive_norm, dims = 1:15)

# Cells cluster
######
#Investigating resolution and how each cluster splits                      ----
naive_norm.cluster = FindNeighbors(naive_norm, dims = 1:10)
naive_norm.cluster = FindClusters(naive_norm.cluster, resolution = seq(0.2, 1, by = 0.5))
naive_norm.cluster = RunUMAP(naive_norm.cluster, dims = 1:10)
DimPlot(naive_norm.cluster, reduction = "umap", split.by = "orig.ident")
clustree(naive_norm.cluster, prefix = "RNA_snn_res.")

pbs_norm.cluster = FindNeighbors(pbs_norm, dims = 1:10)
pbs_norm.cluster = FindClusters(pbs_norm.cluster, resolution = seq(0.2, 1, by = 0.5))
pbs_norm.cluster = RunUMAP(pbs_norm.cluster, dims = 1:10)
DimPlot(pbs_norm.cluster, reduction = "umap", split.by = "orig.ident")
clustree(pbs_norm.cluster, prefix = "RNA_snn_res.")

igg_norm.cluster = FindNeighbors(igg_norm, dims = 1:10)
igg_norm.cluster = FindClusters(igg_norm.cluster, resolution = seq(0.2, 1, by = 0.5))
igg_norm.cluster = RunUMAP(igg_norm.cluster, dims = 1:10)
DimPlot(igg_norm.cluster, reduction = "umap", split.by = "orig.ident")
clustree(igg_norm.cluster, prefix = "RNA_snn_res.")

antipd1_norm.cluster = FindNeighbors(antipd1_norm, dims = 1:10)
antipd1_norm.cluster = FindClusters(antipd1_norm.cluster, resolution = seq(0.2, 1, by = 0.5))
antipd1_norm.cluster = RunUMAP(antipd1_norm.cluster, dims = 1:10)
DimPlot(antipd1_norm.cluster, reduction = "umap", split.by = "orig.ident")
clustree(antipd1_norm.cluster, prefix = "RNA_snn_res.")

antictla4_norm.cluster = FindNeighbors(antictla4_norm, dims = 1:10)
antictla4_norm.cluster = FindClusters(antictla4_norm.cluster, resolution = seq(0.2, 1, by = 0.5))
antictla4_norm.cluster = RunUMAP(antictla4_norm.cluster, dims = 1:10)
DimPlot(antictla4_norm.cluster, reduction = "umap", split.by = "orig.ident")
clustree(antictla4_norm.cluster, prefix = "RNA_snn_res.")

# find all markers of cluster 2 differential expression genes.             ----
##naive
# cluster2.markers_naive_norm = FindMarkers(naive_norm.cluster, ident.1 = 2, min.pct = 0.25)
# head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
# cluster5.markers_naive_norm = FindMarkers(naive_norm.cluster, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
# head(cluster5.markers_naive_norm, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
naive_norm.cluster.new = FindAllMarkers(naive_norm.cluster, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbs_norm.cluster.new = FindAllMarkers(pbs_norm.cluster, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
igg_norm.cluster.new = FindAllMarkers(igg_norm.cluster, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
antipd1_norm.cluster.new = FindAllMarkers(antipd1_norm.cluster, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
antictla4_norm.cluster.new = FindAllMarkers(antictla4_norm.cluster, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# integrate
naive_pbs_norm.cluster = FindAllMarkers(naive_pbs_norm, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)




naive_norm.cluster.new |> 
     group_by(cluster) |>
     slice_max(n = 2, order_by = avg_log2FC)

pbs_norm.cluster.new |> 
   group_by(cluster) |>
   slice_max(n = 2, order_by = avg_log2FC)

igg_norm.cluster.new |> 
   group_by(cluster) |>
   slice_max(n = 2, order_by = avg_log2FC)

antipd1_norm.cluster.new |> 
       group_by(cluster) |>
       slice_max(n = 2, order_by = avg_log2FC)

antictla4_norm.cluster.new |> 
         group_by(cluster) |>
         slice_max(n = 2, order_by = avg_log2FC)

top10 = naive_pbs_norm.cluster |>
  group_by(cluster) |>
  dplyr::top_n(n = 10, wt = avg_log2FC)

subset(top10, cluster == 5)

top10[which(!top10$gene %in% top10[duplicated(top10$gene),"gene"]),] #why is this returning everyone??
#Plot key genes in the ridge plot
RidgePlot(Her2p, assay = "SCT", features = c("IFI27","IFI6"), ncol = 2, group.by = "SCT_snn_res.0.4", cols = color)


VlnPlot(naive_pbs_norm, features = c("Scgb1a1", "cd4", "Foxp3"))
VlnPlot(antictla4_norm.cluster, features = c("Cd3e", "Cd14", "Foxp3"))

FeaturePlot(naive_pbs_norm, features =c("Scgb1a1", "Cd4", "Foxp3"))
FeaturePlot(antictla4_norm.cluster, features =c("Cd3e", "Cd14", "Foxp3"))

naive_norm.cluster.new |>
  group_by(cluster) |>
  top_n(n = 10, wt = avg_log2FC) = top10
head(top10)

DoHeatmap(naive_norm.cluster, features = top10$gene) + NoLegend()

new.cluster.ids.naive = c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Treg", "Macrophage", "Treg foxp3+", "Treg Tbet", "pDC", "Mono/Mk Doublets")
names(new.cluster.ids.naive) = levels(naive_norm.cluster)
naive_norm.cluster = RenameIdents(naive_norm.cluster, new.cluster.ids.naive)
pd = DimPlot(naive_norm.cluster, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
plot_grid(p1, pd)

head(naive_norm.cluster)
Idents(object = naive_norm)








