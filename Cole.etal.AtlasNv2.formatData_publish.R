# Setup ----
# load libraries 
#* you have to install these all the first time and POSSIBLY their dependencies until you get no errors.
library(easypackages)
libraries("Seurat", "Matrix", "readxl","RColorBrewer",'Rmagic',
          'patchwork','dplyr','viridis','ggplot2','pals','SeuratWrappers')

## load workspace with relevant objects ---- 
load (file = 'Genes.Nv2.RData') #gene annotations and colour palettes.

# Generate datasets ----

#* all rawdata available from *GEO* under accession numbers *GSE200198* and *GSE154105*

Developmental = T #load and process samples from developmental stages
Tissues = T # load and process samples from adult tissue dissections
Alldata = T # merge full dataset
if (Developmental == F &&
    Tissues == F &&
    Alldata == F)
{
  if (!exists("AllData.Nv2"))
    load(file='Alldata.Nv2.separate.Robj')
  
  if (!exists("Alldata.dev.Nv2"))
    load (file = 'Nv2.Alldata.dev.RObj')
  
  if (!exists("Alldata.tissue.Nv2"))
    load (file = 'Nv2.Alldata.tissues.RObj')

}

if (Developmental)  
  ## Developmental ----
  {
    ### load data ----
    raw.data1 <- Read10X(data.dir = "Z:/sequencing/Alison/Nematostella/Nv2_mapping/wildtype/18hr.gast/filtered_feature_bc_matrix")
    raw.data2 <- Read10X(data.dir = "Z:/sequencing/Alison/Nematostella/Nv2_mapping/wildtype/24h.gast(2)/filtered_feature_bc_matrix")
    raw.data3 <- Read10X(data.dir = "Z:/sequencing/Alison/Nematostella/Nv2_mapping/wildtype/24h.gast(3)/filtered_feature_bc_matrix")
    raw.data4 <- Read10X(data.dir = "Z:/sequencing/Alison/Nematostella/Nv2_mapping/wildtype/24h.gast(SN)/filtered_feature_bc_matrix")
    raw.data5 <- Read10X(data.dir = "Z:/sequencing/Alison/Nematostella/Nv2_mapping/wildtype/48h.pla(1)/filtered_feature_bc_matrix") #need to re-run this.
    raw.data6 <- Read10X(data.dir = "Z:/sequencing/Alison/Nematostella/Nv2_mapping/wildtype/48h.pla(2)/filtered_feature_bc_matrix")
    raw.data7 <- Read10X(data.dir = "Z:/sequencing/Alison/Nematostella/Nv2_mapping/wildtype/72h.pla(1)/filtered_feature_bc_matrix")
    raw.data8 <- Read10X(data.dir = "Z:/sequencing/Alison/Nematostella/Nv2_mapping/wildtype/4d.pla(1)/filtered_feature_bc_matrix")
    raw.data9 <- Read10X(data.dir = "Z:/sequencing/Alison/Nematostella/Nv2_mapping/wildtype/4d.pla(2cryo)/filtered_feature_bc_matrix")
    raw.data10 <- Read10X(data.dir = "Z:/sequencing/Alison/Nematostella/Nv2_mapping/wildtype/5d.pla(1)/filtered_feature_bc_matrix")
    raw.data11 <- Read10X(data.dir = "Z:/sequencing/Alison/Nematostella/Nv2_mapping/wildtype/8d.ppolyp(1)/filtered_feature_bc_matrix")
    raw.data12 <- Read10X(data.dir = "Z:/sequencing/Alison/Nematostella/Nv2_mapping/wildtype/8d.ppolyp(2)/filtered_feature_bc_matrix")
    raw.data13 <- Read10X(data.dir = "Z:/sequencing/Alison/Nematostella/Nv2_mapping/wildtype/16d.ppolyp(1)/filtered_feature_bc_matrix")
    # 
    # # set the gene names to the annotations
    rownames(raw.data1) <- genes$gene_short_name
    rownames(raw.data2) <- genes$gene_short_name
    rownames(raw.data3) <- genes$gene_short_name
    rownames(raw.data4) <- genes$gene_short_name
    rownames(raw.data5) <- genes$gene_short_name
    rownames(raw.data6) <- genes$gene_short_name
    rownames(raw.data7) <- genes$gene_short_name
    rownames(raw.data8) <- genes$gene_short_name
    rownames(raw.data9) <- genes$gene_short_name
    rownames(raw.data10) <- genes$gene_short_name
    rownames(raw.data11) <- genes$gene_short_name
    rownames(raw.data12) <- genes$gene_short_name
    rownames(raw.data13) <- genes$gene_short_name
    
    #### generate Seurat objects ----   
    
    h18  <- CreateSeuratObject(counts = raw.data1, project = "D.18h gastrula")
    h24.2  <- CreateSeuratObject(counts = raw.data2, project = "D.24h2 gastrula")
    h24.3  <- CreateSeuratObject(counts = raw.data3, project = "D.24h3 gastrula")
    h24.4  <- CreateSeuratObject(counts = raw.data4, project = "D.24h4 gastrula")
    d2  <- CreateSeuratObject(counts = raw.data5, project = "D.2d planula")
    d2.2  <- CreateSeuratObject(counts = raw.data6, project = "D.2d2 planula")
    d3  <- CreateSeuratObject(counts = raw.data7, project = "D.3d planula")
    d4  <- CreateSeuratObject(counts = raw.data8, project = "D.4d planula")
    d4c  <- CreateSeuratObject(counts = raw.data9, project = "D.4dc planula")
    d5  <- CreateSeuratObject(counts = raw.data10, project = "D.5d planula")
    d8.1  <- CreateSeuratObject(counts = raw.data11, project = "D.8d1 polyp")
    d8.2  <- CreateSeuratObject(counts = raw.data12, project = "D.8d2 polyp")
    d16  <- CreateSeuratObject(counts = raw.data13, project = "D.16d polyp")
    
    #re-name the samples to make barcodes unique and identifiable
    h18  <- RenameCells(h18, add.cell.id ='h18')
    h24.2  <- RenameCells(h24.2, add.cell.id ='h24.2')
    h24.3  <- RenameCells(h24.3, add.cell.id ='h24.3')
    h24.4  <- RenameCells(h24.4, add.cell.id ='h24.4')
    d2  <- RenameCells(d2, add.cell.id ='d2')
    d2.2  <- RenameCells(d2.2, add.cell.id ='d2.2')
    d3  <- RenameCells(d3, add.cell.id ='d3')
    d4  <- RenameCells(d4, add.cell.id ='d4')
    d4c  <- RenameCells(d4c, add.cell.id ='d4c')
    d5  <- RenameCells(d5, add.cell.id ='d5')
    d8.1  <- RenameCells(d8.1, add.cell.id ='d8.1')
    d8.2  <-RenameCells(d8.2, add.cell.id ='d8.2')
    d16  <- RenameCells(d16, add.cell.id ='d16')
    
    #### filter the libraries ----
    VlnPlot(object = h18, features = c("nFeature_RNA", "nCount_RNA"),
            group.by = 'orig.ident',cols = LibCP)
    h18 <- subset(x = h18, subset = nFeature_RNA > 350 & nCount_RNA < 100000)
    
    VlnPlot(object = h24.2, features = c("nFeature_RNA", "nCount_RNA"),
            group.by = 'orig.ident',cols = LibCP)
    h24.2 <- subset(x = h24.2, subset = nFeature_RNA > 250 & nCount_RNA < 4000)
    
    VlnPlot(object = h24.3, features = c("nFeature_RNA", "nCount_RNA"),
            group.by = 'orig.ident',cols = LibCP)
    h24.3 <- subset(x = h24.3, subset = nFeature_RNA > 250 & nCount_RNA < 15000)
    
    VlnPlot(object = h24.4, features = c("nFeature_RNA", "nCount_RNA"),
            group.by = 'orig.ident',cols = LibCP)
    h24.4 <- subset(x = h24.4, subset = nFeature_RNA > 350 & nCount_RNA < 20000)
    
    VlnPlot(object = d2.2, features = c("nFeature_RNA", "nCount_RNA"),
            group.by = 'orig.ident',cols = LibCP)
    d2.2 <- subset(x = d2.2, subset = nFeature_RNA > 350 & nCount_RNA < 35000)
    
    VlnPlot(object = d2, features = c("nFeature_RNA", "nCount_RNA"),
            group.by = 'orig.ident',cols = LibCP)
    d2 <- subset(x = d2, subset = nFeature_RNA > 250 & nCount_RNA < 50000)
    
    VlnPlot(object = d3, features = c("nFeature_RNA", "nCount_RNA"),
            group.by = 'orig.ident',cols = LibCP)
    d3 <- subset(x = d3, subset = nFeature_RNA > 250 & nCount_RNA < 15000)
    
    VlnPlot(object = d4, features = c("nFeature_RNA", "nCount_RNA"),
            group.by = 'orig.ident',cols = LibCP)
    d4 <- subset(x = d4, subset = nFeature_RNA > 250 & nCount_RNA < 25000)
    
    VlnPlot(object = d4c, features = c("nFeature_RNA", "nCount_RNA"),
            group.by = 'orig.ident',cols = LibCP)
    d4c <- subset(x = d4c, subset = nFeature_RNA > 250 & nCount_RNA < 25000)
    
    VlnPlot(object = d5, features = c("nFeature_RNA", "nCount_RNA"),
            group.by = 'orig.ident',cols = LibCP)
    d5 <- subset(x = d5, subset = nFeature_RNA > 250 & nCount_RNA < 10000)
    
    VlnPlot(object = d8.1, features = c("nFeature_RNA", "nCount_RNA"),
            group.by = 'orig.ident',cols = LibCP)
    d8.1 <- subset(x = d8.1, subset = nFeature_RNA > 250 & nCount_RNA < 10000)
    
    VlnPlot(object = d8.2, features = c("nFeature_RNA", "nCount_RNA"),
            group.by = 'orig.ident',cols = LibCP)
    d8.2 <- subset(x = d8.2, subset = nFeature_RNA > 250 & nCount_RNA < 10000)
    
    VlnPlot(object = d16, features = c("nFeature_RNA", "nCount_RNA"),
            group.by = 'orig.ident',cols = LibCP)
    d16 <- subset(x = d16, subset = nFeature_RNA > 250 & nCount_RNA < 5000)
    
    ### generate the merged object ----
    libraries = c(h18,h24.2,h24.3,h24.4,d2,d2.2,d3,d4,d4c,d5,d8.1,d8.2,d16)
    
    Alldata.dev.Nv2=merge(h18,c(h24.2,h24.3,h24.4,d2,d2.2,d3,d4,d4c,d5,d8.1,d8.2,d16))
    Alldata.dev.Nv2@meta.data$orig.ident<-as.factor(Alldata.dev.Nv2@meta.data$orig.ident)
    
    libs=c("D.18h gastrula","D.24h2 gastrula", "D.24h3 gastrula", "D.24h4 gastrula",
           "D.2d planula","D.2d2 planula", "D.3d planula",  "D.4d planula",  "D.4dc planula",
           "D.5d planula","D.8d1 polyp","D.8d2 polyp", 'D.16d polyp')
    order.ind = match(libs,levels(Alldata.dev.Nv2@meta.data$orig.ident))
    Alldata.dev.Nv2@meta.data$orig.ident = factor(Alldata.dev.Nv2@meta.data$orig.ident,
                                                  levels(Alldata.dev.Nv2@meta.data$orig.ident)[order.ind])#
    Alldata.dev.Nv2$lifehistory='DevSubset'
    #### clean up the workspace ----  
    rm(h18,h24.2,h24.3,h24.4,d2,d2.2,d3,d4,d4c,d5,d8.1,d8.2,d16)
    rm (raw.data1, raw.data2, raw.data3,raw.data4,raw.data5,raw.data6,raw.data7,
        raw.data8,raw.data9,raw.data10,raw.data11,raw.data12,raw.data13)
    
    save(Alldata.dev.Nv2,file='Nv2.Alldata.dev.Robj')
  }


if (Tissues)
  
  ## Tissues ----
{
  ### load raw data ----
{
  raw.data1 <- Read10X(data.dir = "Z:/sequencing/Alison/Nematostella/Nv2_mapping/wildtype/Tissues4.tent/filtered_feature_bc_matrix")
  raw.data2 <- Read10X(data.dir = "Z:/sequencing/Alison/Nematostella/Nv2_mapping/wildtype/Tissues4.pha/filtered_feature_bc_matrix")
  raw.data3 <- Read10X(data.dir = "Z:/sequencing/Alison/Nematostella/Nv2_mapping/wildtype/Tissues4.mes/filtered_feature_bc_matrix")
  raw.data4 <- Read10X(data.dir = "Z:/sequencing/Alison/Nematostella/Nv2_mapping/wildtype/Tissues4.bw/filtered_feature_bc_matrix")
  raw.data5 <- Read10X(data.dir = "Z:/sequencing/Alison/Nematostella/Nv2_mapping/wildtype/Tissue3.intPbw/filtered_feature_bc_matrix")
  raw.data6 <- Read10X(data.dir = "Z:/sequencing/Alison/Nematostella/Nv2_mapping/wildtype/Tissue2.phbw/filtered_feature_bc_matrix")
  raw.data7 <- Read10X(data.dir = "Z:/sequencing/Alison/Nematostella/Nv2_mapping/wildtype/Tissue1.mes.M/filtered_feature_bc_matrix")
  raw.data8 <- Read10X(data.dir = "Z:/sequencing/Alison/Nematostella/Nv2_mapping/wildtype/Tissue1.mes.F/filtered_feature_bc_matrix")
  
  # # set the gene names to the annotations | needs updating to Nv2
  rownames(raw.data1) <- genes$gene_short_name
  rownames(raw.data2) <- genes$gene_short_name
  rownames(raw.data3) <- genes$gene_short_name
  rownames(raw.data4) <- genes$gene_short_name
  rownames(raw.data5) <- genes$gene_short_name
  rownames(raw.data6) <- genes$gene_short_name
  rownames(raw.data7) <- genes$gene_short_name
  rownames(raw.data8) <- genes$gene_short_name
  
  ### generate Seurat object ----   
  
  t4.tent  <- CreateSeuratObject(counts = raw.data1, project = "T.tentacle1")
  t4.pha  <- CreateSeuratObject(counts = raw.data2, project = "T.pharynx1")
  t4.mes  <- CreateSeuratObject(counts = raw.data3, project = "T.mesentery1")
  t4.bw  <- CreateSeuratObject(counts = raw.data4, project = "T.bodywall1")
  t3.intPbw  <- CreateSeuratObject(counts = raw.data5, project = "T.intP.bw")
  t2.phbw  <- CreateSeuratObject(counts = raw.data6, project = "T.phbw")
  t1.mesM  <- CreateSeuratObject(counts = raw.data7, project = "T.mes.M")
  t1.mesF  <- CreateSeuratObject(counts = raw.data8, project = "T.mes.F")
  
  #re-name the samples to make barcodes unique and identifiable
  
  t4.tent  <- RenameCells(t4.tent, add.cell.id ='tent')
  t4.pha  <- RenameCells(t4.pha, add.cell.id ='pha')
  t4.mes  <- RenameCells(t4.mes, add.cell.id ='mes.j')
  t4.bw  <- RenameCells(t4.bw, add.cell.id ='bw')
  t3.intPbw  <- RenameCells(t3.intPbw, add.cell.id ='intp.bw')
  t2.phbw  <- RenameCells(t2.phbw, add.cell.id ='phbw')
  t1.mesM  <- RenameCells(t1.mesM, add.cell.id ='mesM')
  t1.mesF  <- RenameCells(t1.mesF, add.cell.id ='mesF')
  
  ####filter the libraries: ----
  VlnPlot(object = t4.tent, features = c("nFeature_RNA", "nCount_RNA"),
          group.by = 'orig.ident',cols = LibCP)
  t4.tent <- subset(x = t4.tent, subset = nFeature_RNA > 250 & nCount_RNA < 10000)
  
  VlnPlot(object = t4.pha, features = c("nFeature_RNA", "nCount_RNA"),
          group.by = 'orig.ident',cols = LibCP)
  t4.pha <- subset(x = t4.pha, subset = nFeature_RNA > 250 & nCount_RNA < 10000)
  
  VlnPlot(object = t4.mes, features = c("nFeature_RNA", "nCount_RNA"),
          group.by = 'orig.ident',cols = LibCP)
  t4.mes <- subset(x = t4.mes, subset = nFeature_RNA > 250 & nCount_RNA < 10000)
  
  VlnPlot(object = t4.bw, features = c("nFeature_RNA", "nCount_RNA"),
          group.by = 'orig.ident',cols = LibCP)
  t4.bw <- subset(x = t4.bw, subset = nFeature_RNA > 250 & nCount_RNA < 6000)
  
  
  VlnPlot(object = t3.intPbw, features = c("nFeature_RNA", "nCount_RNA"),
          group.by = 'orig.ident',cols = LibCP)
  t3.intPbw <- subset(x = t3.intPbw, subset = nFeature_RNA > 250 & nCount_RNA < 5000)
  
  VlnPlot(object = t2.phbw, features = c("nFeature_RNA", "nCount_RNA"),
          group.by = 'orig.ident',cols = LibCP)
  t2.phbw <- subset(x = t2.phbw, subset = nFeature_RNA > 250 & nCount_RNA < 15000)
  
  VlnPlot(object = t1.mesM, features = c("nFeature_RNA", "nCount_RNA"),
          group.by = 'orig.ident',cols = LibCP)
  t1.mesM <- subset(x = t1.mesM, subset = nFeature_RNA > 350 & nCount_RNA < 10000)
  
  VlnPlot(object = t1.mesF, features = c("nFeature_RNA", "nCount_RNA"),
          group.by = 'orig.ident',cols = LibCP)
  t1.mesF<- subset(x = t1.mesF, subset = nFeature_RNA > 350 & nCount_RNA < 10000)
  

}
  

#generate the merged object
{
Alldata.tissue.Nv2=merge(t4.tent ,c(t4.pha,t4.mes,t4.bw  ,t3.intPbw,t2.phbw,t1.mesM,t1.mesF))
Alldata.tissue.Nv2@meta.data$orig.ident<-as.factor(Alldata.tissue.Nv2@meta.data$orig.ident)

libs=c("T.tentacle1","T.pharynx1","T.mesentery1","T.bodywall1","T.intP.bw","T.phbw","T.mes.M","T.mes.F")
order.ind = match(libs,levels(Alldata.tissue.Nv2@meta.data$orig.ident))
Alldata.tissue.Nv2@meta.data$orig.ident = factor(Alldata.tissue.Nv2@meta.data$orig.ident,
                                                 levels(Alldata.tissue.Nv2@meta.data$orig.ident)[order.ind])#
#clean up the workspace  
rm(t4.tent ,t4.pha,t4.mes,t4.bw  ,t3.intPbw,t2.phbw,t1.mesM,t1.mesF)
rm (raw.data1, raw.data2, raw.data3,raw.data4,raw.data5,raw.data6,raw.data7,
    raw.data8,raw.data9,raw.data10,raw.data11,raw.data12,raw.data13)

# save (Alldata.tissue.Nv2,file='Alldata.tissue.Nv2.unanalyzed.RObj')
data1=Alldata.tissue.Nv2
data1@meta.data$orig.ident<-as.factor(data1@meta.data$orig.ident)

#print the levels:
levels(data1$orig.ident)

#order by eye (or write out them out and run a match index script)
order.ind=c(4:6,2,1,3,7:8)

#re-order the factors as 
data1@meta.data$orig.ident = factor(data1@meta.data$orig.ident,
                                    levels(data1@meta.data$orig.ident)[order.ind])#
}
### Run Seurat pipeline ----
{
  #### normalize the dataset (merged)  ----
  data1 <- NormalizeData(data1, scale.factor = 5000) 
  #this can be discussed; most people will use 10,000 (Seurat default I think) but this might be over inflating.
  
  #### calculate variable genes ----
  #select variable features from each library and use this set.
  {
    list=  NULL
    vargenelist <- SplitObject(data1, split.by = "orig.ident")
    for (i in 1:length(vargenelist)) {
      vargenelist[[i]] <- FindVariableFeatures(vargenelist[[i]],nfeatures = 500, verbose = FALSE)
    }
    for (i in 1:length(vargenelist)) {
      x <- vargenelist[[i]]@assays$RNA@var.features
      list=c(list,x)}
    list=unique(list)
    length(list)
    
    data1@assays$RNA@var.features = list
    
  }
  data1 <- ScaleData(data1,split.by = 'orig.ident')#
  
  #### Dimensional reductions ----    
  data1 <- RunPCA(data1, pcs.compute = 50)
  PCAPlot(data1,group.by='orig.ident')
  ElbowPlot(object = data1, ndims = 50)
  
  # set dimensions
  d=as.integer(which(data1@reductions$pca@stdev>1.5))
  d=1:30
  #UMAP
  data1 <- RunUMAP(data1, dims = d,
                   reduction = 'pca',
                   reduction.name ='umap',reduction.key ='umap', #this is default; can change name to save different maps
                   #the following parameters control your output.
                   n.neighbors = 20L, #how many similar cells do you expect
                   spread =1, #how big is the axis space: 0.2 for lineages; 1 for separation
                   min.dist = 0.6, #how close to plot each cell higher=spread
                   local.connectivity = 1)#overall how connected is the graph (AllData = 100)
}
#### Image cellplots ----
  DimPlot(data1,cols=alpha(LibCP[c(1,3,5,7,9,11,13,15,17)],0.8),group.by = 'orig.ident',order=rev(levels(data1$orig.ident)))&NoAxes()

## save Tissues ----
Alldata.tissue.Nv2 = data1
# add life cycle information:
Alldata.tissue.Nv2$lifehistory = 'NA'
Alldata.tissue.Nv2$lifehistory[colnames(Alldata.tissue.Nv2)]='AdultSubset' 

 
save(Alldata.tissue.Nv2,file='Nv2.Alldata.tissues.Robj')
}

if (Alldata)
  # Alldata ----
{
  if(!exists("Alldata.dev.Nv2"))
    load(file='Nv2.Alldata.dev.Robj')
  if(!exists("Alldata.tissue.Nv2"))
     load(file='Nv2.Alldata.tissues.Robj')
  data1=merge(Alldata.dev.Nv2,Alldata.tissue.Nv2)
  data1@meta.data$orig.ident<-as.factor(data1@meta.data$orig.ident)
  
  #print the levels:
  levels(data1$orig.ident)
  
  #order by eye (or write out them out and run a match index script)
  order.ind=c(2:13,1,21,14:20)
  
  #re-order the factors as 
  data1@meta.data$orig.ident = factor(data1@meta.data$orig.ident,
                                      levels(data1@meta.data$orig.ident)[order.ind])#
}

Alldata.Nv2=data1

# add life cycle information

Alldata.Nv2$lifehistory = 'NA'
Alldata.Nv2$lifehistory[colnames(Alldata.dev.Nv2)]='DevSubset'
Alldata.Nv2$lifehistory[colnames(Alldata.tissue.Nv2)]='AdultSubset'

#save workspace ----
rm(vargenelist,d,i,libs,list,order.ind,libraries,data1,AllData.Nv2,x,Tissues,Alldata,Developmental,Alldata.tissue.Nv2, Alldata.dev.Nv2)
save.image(file="Cole.Atlas.Nv2.preprocess.RData")
