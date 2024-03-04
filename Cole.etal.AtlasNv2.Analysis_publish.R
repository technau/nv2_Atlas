# Setup ----
# load libraries 
library(easypackages)
libraries("Seurat", "Matrix", "readxl","RColorBrewer",'Rmagic',
          'patchwork','dplyr','viridis','ggplot2','pals','SeuratWrappers')

## load workspace with relevant objects ---- 
load (file = 'Cole.Atlas.Nv2.preprocess.RData') 

## process Alldata ----
data1=Alldata.Nv2

  ###standard Seurat pipleline ----
  {
    #normalize the dataset (merged)  
    data1 <- NormalizeData(data1, scale.factor = 5000) 

    #### calculate variable genes ----
    #select variable features from each library and use this set.
    {
      list=  NULL
      vargenelist <- SplitObject(data1, split.by = "orig.ident")
      for (i in 1:length(vargenelist)) {
        vargenelist[[i]] <- FindVariableFeatures(vargenelist[[i]],nfeatures = 1000, verbose = FALSE)
      }
      for (i in 1:length(vargenelist)) {
        x <- vargenelist[[i]]@assays$RNA@var.features
        list=c(list,x)}
      list=unique(list)
      length(list)
      
      data1@assays$RNA@var.features = list
      
    }
    data1 <- ScaleData(data1,split.by = 'orig.ident')#
    
    ### Dimensional reductions ----    
    data1 <- RunPCA(data1, pcs.compute = 50)
    PCAPlot(data1,group.by='orig.ident')
    ElbowPlot(object = data1, ndims = 50)
    
    # set dimensions
    d=1:30
    #UMAP
    data1 <- RunUMAP(data1, dims = d,
                     reduction = 'pca',
                     reduction.name ='umap',reduction.key ='umap',n.neighbors = 10L, 
                     spread =0.6, 
                     min.dist = 0.3, 
                     local.connectivity = 10)
    
    ### Find Clusters ----
    data1 <- FindNeighbors(object = data1,reduction ="pca",dims = d,
                           annoy.metric = 'cosine',
                           k.param = 20)
    
    data1 <- FindClusters(object = data1,resolution = 0.2)
    
    #calculate relationship between clusters
    data1 <- BuildClusterTree(object = data1, reorder = TRUE,
                              dims = d,
                              reorder.numeric = T    )

    #image cell plot and clusters    
    DimPlot(data1,reduction='umap', label = T,label.size = 5, repel = T,cols = clust.cp.separate)&
      NoLegend()&NoAxes()&
    DimPlot(data1,cols=alpha(LibCP,0.6),group.by = 'orig.ident',order=rev(levels(data1$orig.ident)))&NoAxes()
    
  }

  ## assign cluster names ----
  #run a semi-automated script for updating cluster names based on marker gene expression
  {
    clusterNames<- read_excel("SupplementalData2.xlsx", sheet='Tissues.Nv2')
    goi = clusterNames$Marker
    
    data1 <- SetIdent(data1, value = 'seurat_clusters') #if cluster names are not integers
    data1 <- BuildClusterTree(data1,reorder = T,reorder.numeric = T) #this re-sets everything to numeric; required by script below
    DotPlot(data1,features= goi)+RotatedAxis()+NoLegend()& #checks that all your genes are available and actually express specifically; useful for building the gene list
      DimPlot(data1,label=T,cols = clust.cp.separate)+NoAxes()
    #assign cluster ID to the individual libraries
    data1<-ScaleData(data1,features = goi, split.by = 'orig.ident') #works with scaled gene values, so must scale all of these first (in case they are not within the VarFeatures set!)
    
    cl <-length(levels(data1@active.ident)) #how many clusters
    C.suffix <-seq(1:cl) #the numbers
    
    g=length(goi)
    clName = vector()
    m=matrix(0L,g,cl) #generate a matrix of values of each cluster for each gene:
    for (j in 1:cl) #for each cluster set
    {
      for (i in 1:g) #for each gene
        m[i,j]=mean(data1@assays$RNA@scale.data[goi[i],WhichCells(data1,idents = C.suffix[j])]) #average scaled value for cluster/gene
      clName[j]=as.integer(which.max(m[,j])) #choose the highest value
    }
    sort(clName) #this prints which IDs were selected; sanity check that it did what you wanted... should have no unwanted replicated IDs
    clusterNames$ID[clName] #you can also look at the names to see where the problem might have been
    
    #use the desired order from the spreadsheet to re-order the clusters:
    #first order the identities..
    data1@active.ident = factor(data1@active.ident,
                                levels(data1@active.ident)[order(clName)])
    #set the names to your IDs; once the gene list works, it is really easy to then update your cluster names with your excel spreadsheet.
    levels(data1@active.ident) = clusterNames$ID[clName][order(clName)]
    #save the IDs in metadata:
    data1@meta.data$IDs = data1@active.ident
  }
  ## save  ----
  Alldata.Nv2=data1
  DimPlot(Alldata.Nv2,group.by = 'IDs',
          cols=clust.cp.separate)&NoAxes()
  Alldata.Nv2<-SetIdent(Alldata.Nv2,value = 'IDs')
  save (Alldata.Nv2,file='Alldata.Nv2.publish.Robj')
  
## Generate DEG lists ----
data1=SetIdent(Alldata.Nv2,value='IDs')

{
  data1@active.assay='RNA'
  all.markers <- FindAllMarkers(data1,
                                logfc.threshold = 1,
                                return.thresh = 0.001,
                                min.pct = 0.2,
                                only.pos = TRUE)

    # add annotations associated with this list:
    l=c(names(all.markers),colnames(genes)[c(1,4:6,3)])
    all.markers[,8:12]<-'NA'
    all.markers=setNames(all.markers,l) 
    # 
    for (i in 1:length(levels(data1@active.ident))) # 
    {
      x=all.markers[as.numeric(all.markers$cluster)==i,][1:length(which(as.numeric(all.markers$cluster)==i)),7]
      anInd = match(genes[match(x,genes$gene_short_name),1],genes$geneID)
      all.markers[as.numeric(all.markers$cluster)==i,][1:length(which(as.numeric(all.markers$cluster)==i)),8:12]<-genes[anInd,c(1,4:6,3)]
    }  

  
  #also for transcription factors only: increase the number of available genes
  data1=FindVariableFeatures(data1,nfeatures = 5000)
  
  all.markers_TF <- FindAllMarkers(data1,
                                   logfc.threshold = 0.4,
                                   features = intersect(TF_list,data1@assays$RNA@var.features),
                                   min.pct = 0.05,
                                   only.pos = TRUE,
                                   return.thresh = 0.001)
  
  
    l=c(names(all.markers_TF),colnames(genes)[c(1,4:6,3)])
    all.markers_TF[,8:12]<-'NA'
    all.markers_TF=setNames(all.markers_TF,l) 
    # 
    for (i in 1:length(levels(data1@active.ident))) # 
    {
      x=all.markers_TF[as.numeric(all.markers_TF$cluster)==i,][1:length(which(as.numeric(all.markers_TF$cluster)==i)),7]
      anInd = match(genes[match(x,genes$gene_short_name),1],genes$geneID)
      all.markers_TF[as.numeric(all.markers_TF$cluster)==i,][1:length(which(as.numeric(all.markers_TF$cluster)==i)),8:12]<-genes[anInd,c(1,4:6,3)]
    }  
}
  
alldata.markers=all.markers
alldata.markers.TF=all.markers_TF

save(alldata.markers,alldata.markers.TF,file='AlldataDEGS.updated.RData')

#process subsets ----
  subsets.single = F
  if (subsets.single)
  {
    #First, merge the neurons and neurogland, and the cnidocytes:
    levels(Alldata.Nv2$IDs)=c("pSC", "PGCs", "neurogland.all", "gland.mucous", "unchar.immune","neurogland.all", "cnidocyte","cnidocyte", "pharyngeal.ect", "epithelia.ect","ectoderm.embryonic", "retractor muscle", "gastrodermis", "mesendoderm.embryonic")

    Alldata.Nv2=SetIdent(Alldata.Nv2,value='IDs')
    Alldata.Nv2=droplevels(Alldata.Nv2)
    
    
    data1.subsets <- SplitObject(Alldata.Nv2)
    all.analyses = T
    order = match(levels(Alldata.Nv2),names(data1.subsets))
    data1.subsets=data1.subsets[order]
    
    ## subset analyses ----
    if (all.analyses)
    {
      ### analysis ----
      
      for (i in 1:length(names (data1.subsets)))
      {
        data1=data1.subsets[[i]]
        data1@active.assay = 'RNA'
        data1 <- FindVariableFeatures(data1,nfeatures = 2000)#,
        
        coi = as.vector(data1@assays$RNA@counts@Dimnames[[2]])
        t=ScaleData(Alldata.Nv2,split.by = 'orig.ident', features = data1@assays$RNA@var.features)
        data1@assays$RNA@scale.data = t@assays$RNA@scale.data[,coi]
       
        # run different reduction algorythms:
        # #PCA
        data1 <- RunPCA(data1, pcs.compute = 50)
        ElbowPlot(object = data1, ndims = 50)
        
        # set dimensions
        d=as.integer(which(data1@reductions$pca@stdev>2))
        if (length(d) <10 || length(d)>=30)
          d=1:20
        
        data1 <- RunUMAP(data1, reduction ="pca", 
                         n.neighbors = 10L,spread =1,
                         dims = d,reduction.name ='umap',
                         reduction.key ='umap',min.dist = 0.3,
                         local.connectivity = 10)#
        data1 <- FindNeighbors(object = data1,
                               reduction ="pca",dims = d,
                               nn.method = 'annoy',  
                               annoy.metric = 'cosine',
                               k.param = 10)
        data1 <- FindClusters(object = data1,resolution = 1,random.seed = 0)
        # #look at relationship between clusters
        data1 <- BuildClusterTree(object = data1, reorder = TRUE,
                                  reorder.numeric = TRUE)#, dims = c(1:d))
        data1$IDs.fine = data1@active.ident
        data1 <- FindClusters(object = data1,resolution = 0.2,random.seed = 0)
        # #look at relationship between clusters
        data1 <- BuildClusterTree(object = data1, reorder = TRUE,
                                  reorder.numeric = TRUE)#, dims = c(1:d))
        data1$IDs.coarse = data1@active.ident
        
        data1.subsets[[i]] = data1
      }#this sets RNA clustering at res = 0.2 as default
      save(data1.subsets,file='data1.subsets.scalealldata.Robj')
 ## CytoTRACE ----
#This wasn't saved... need to re-run
      for (i in 1:length(names (data1.subsets)))
      {
          library(CytoTRACE)
          data1 = data1.subsets[[i]]
          data1@active.assay='RNA'
          cyto<-CytoTRACE(as.matrix(data1@assays$RNA@counts))#, batch = data1@meta.data$orig.ident)
          data1@meta.data$cytoTRACE = cyto$CytoTRACE
          data1@meta.data$cytoTRACEorder = cyto$CytoTRACErank
          print(FeaturePlot(data1,'cytoTRACE', pt.size = 2, reduction = 'umap', cols = rev(brewer.pal(11 , "Spectral" ))))
          data1.subsets[[i]]=data1
        }

    }  
  }

#Check and annotate clusters from subsets ----
  ## cluster Ident pSC ----
  data1=data1.subsets$pSC
  DimPlot(data1,cols=clust.cp.separate)
  {
    
    DimPlot(data1,label = T,cols=clust.cp.separate,group.by = 'IDs.fine')&NoAxes()&
      DimPlot(data1,label = T,cols=clust.cp.separate,group.by = 'IDs.coarse')&NoAxes()
    data1<- SetIdent(data1,value = 'IDs.coarse')
    
    #run a semi-automated script for updating cluster names based on marker gene expression
    # try generating this by 
    generate.lists=F
    if (generate.lists){
      #endodermal:
      x=DotPlot(data1,'RNA','Six1-2')
      cl.Six12=which(x$data$pct.exp>=30)
      (cl.Six12)
      g=NULL
      for (i in 1:length(cl.Six12)){
        x=all.markers$gene[all.markers$cluster==cl.Six12[i]][2]
        g=c(g,x)
      }  
      g
      #used that to populate the marker list
      
      #mitotic:
      x=DotPlot(data1,'RNA','PCNA')
      x
      cl.tll=which(x$data$pct.exp>=60)
      (cl.tll)
      g=NULL
      for (i in 1:length(cl.tll)){
        x=all.markers$gene[all.markers$cluster==cl.tll[i]][1]
        g=c(g,x)
      }  
      g
  #then hand picked the first DEG from each remaining cluster... 
    }
    
    {
      data1<-SetIdent(data1,value = 'IDs.coarse')
      clusterNames<- read_excel("SupplementalData2.xlsx",sheet = 'Alldata.pSC.Nv2')
      goi = clusterNames$Marker
      DimPlot(data1,label=T,cols = clust.cp.separate)+NoAxes()
      #assign cluster ID to the individual libraries
      data1<-ScaleData(data1,features = goi, split.by = 'orig.ident') #works with scaled gene values, so must scale all of these first (in case they are not within the VarFeatures set!)
      
      #set empty name variable:
      clName = vector()
      
      #generate a matrix of values of each cluster for each gene:
      m=AverageExpression(data1,features=goi,slot='scale.data') #Seurat has a function for this
      
      for (j in 1:length(levels(data1@active.ident))) #for each cluster set
      {
        clName[j]=as.integer(which.max(m$RNA[,j])) #choose the highest value
      }
      sort(clName) #this prints which IDs were selected; sanity check that it did what you wanted... should have no unwanted replicated IDs
      clusterNames$ID[clName] #you can also look at the names to see where the problem might have been
      
      #first order the identities..
      data1@active.ident = factor(data1@active.ident,
                                  levels(data1@active.ident)[order(clName)])
      #set the names to your IDs; once the gene list works, it is really easy to then update your cluster names with your excel spreadsheet.
      levels(data1@active.ident) = clusterNames$ID[clName][order(clName)]
      #save the IDs in metadata:
      data1@meta.data$IDs = data1@active.ident

    }
  }
  DimPlot(data1,cols=clust.cp.separate)
  data1.subsets$pSC=data1
  
  DotPlot(data1.subsets$pSC,'RNA',unique(c('SoxC','AshC','INSM1-like-1','AshA',
                                           'Six1-2','NKx2.2D','ZNF845','Nem64')),
                                           # goi[sort(clName)])),
          cols=gene.cp[c(2,11)],scale.by = 'size',col.min = 0
  )&RotatedAxis()
  
  
  
  ## cluster Ident PGC ----
  data1=data1.subsets$PGCs
  data1<-SetIdent(data1,value = 'IDs.coarse')
  DimPlot(data1,cols=clust.cp.separate)
  
  #run a semi-automated script for updating cluster names based on marker gene expression
  {
    clusterNames<- read_excel("SupplementalData2.xlsx",sheet = 'Alldata.PGCs.Nv2')
    goi = clusterNames$Marker
    
    data1 <- SetIdent(data1, value = 'seurat_clusters') #if cluster names are not integers
    data1 <- BuildClusterTree(data1,reorder = T,reorder.numeric = T) #this re-sets everything to numeric; required by script below
    # DotPlot(data1,features= goi)+RotatedAxis()+NoLegend() #checks that all your genes are available and actually express specifically; useful for building the gene list
    DimPlot(data1,label=T,cols = clust.cp.separate)+NoAxes()
    #assign cluster ID to the individual libraries
    data1<-ScaleData(data1,features = goi, split.by = 'orig.ident') #works with scaled gene values, so must scale all of these first (in case they are not within the VarFeatures set!)
    
    #set empty name variable:
    clName = vector()
    
    #generate a matrix of values of each cluster for each gene:
    m=AverageExpression(data1,features=goi,slot='scale.data') #Seurat has a function for this
    
    for (j in 1:length(levels(data1@active.ident))) #for each cluster set
    {
      clName[j]=as.integer(which.max(m$RNA[,j])) #choose the highest value
    }
    sort(clName) #this prints which IDs were selected; sanity check that it did what you wanted... should have no unwanted replicated IDs
    clusterNames$ID[clName] #you can also look at the names to see where the problem might have been
    
    #first order the identities..
    data1@active.ident = factor(data1@active.ident,
                                levels(data1@active.ident)[order(clName)])
    #set the names to your IDs; once the gene list works, it is really easy to then update your cluster names with your excel spreadsheet.
    levels(data1@active.ident) = clusterNames$ID[clName][order(clName)]
    #save the IDs in metadata:
    data1@meta.data$IDs = data1@active.ident
    data1.subsets$PGCs=data1
    
    DotPlot(data1.subsets$PGCs,'RNA',unique(c('PIWL1-like-5','PIWL2-like-1','BOULE-like-1',
                                              'Nanos2',
                                              'Nanos1','Sox3','SoxC','FoxO3',
                                              'FOXO-like-3','Myc2','myc1','Myc3',
                                              goi[sort(clName)])),
            cols=gene.cp[c(2,11)],scale.by = 'size',col.min = 0
    )&RotatedAxis()
    
    #Then update markers with 
    j=2
  }
  DimPlot(data1,cols=clust.cp.separate)

  ## cluster Ident RM ----
  data1=data1.subsets$`retractor muscle`
  data1<-SetIdent(data1,value = 'IDs.coarse')
  
  DimPlot(data1,cols=clust.cp.separate)
  
  #run a semi-automated script for updating cluster names based on marker gene expression
  
  {
    
    # try generating this by 
    generate.lists=F
    if (generate.lists){
      # tentacle muscle
      x=DotPlot(data1,'RNA','Nem64')
      x
      cl.Six12=which(x$data$pct.exp>=50)
      (cl.Six12)
      g=NULL
      for (i in 1:length(cl.Six12)){
        x=all.markers$gene[all.markers$cluster==cl.Six12[i]][1]
        g=c(g,x)
      }  
      g
      
      #mesentery muscle
      x=DotPlot(data1,'RNA','Nem24')
      x
      cl.tll=which(x$data$pct.exp>=2)
      (cl.tll)
      g=NULL
      for (i in 1:length(cl.tll)){
        x=all.markers$gene[all.markers$cluster==cl.tll[i]][1]
        g=c(g,x)
      }  
      g
      
    }
    
    
    {clusterNames<- read_excel("SupplementalData2.xlsx",sheet = 'Alldata.RM.Nv2')
      goi = clusterNames$Marker
      
       DimPlot(data1,label=T,cols = clust.cp.separate)+NoAxes()
      #assign cluster ID to the individual libraries
      data1<-ScaleData(data1,features = goi, split.by = 'orig.ident') #works with scaled gene values, so must scale all of these first (in case they are not within the VarFeatures set!)
      
      #set empty name variable:
      clName = vector()
      
      #generate a matrix of values of each cluster for each gene:
      m=AverageExpression(data1,features=goi,slot='scale.data') #Seurat has a function for this
      
      for (j in 1:length(levels(data1@active.ident))) #for each cluster set
      {
        clName[j]=as.integer(which.max(m$RNA[,j])) #choose the highest value
      }
      sort(clName) #this prints which IDs were selected; sanity check that it did what you wanted... should have no unwanted replicated IDs
      clusterNames$ID[clName] #you can also look at the names to see where the problem might have been
      
      
      #first order the identities..
      data1@active.ident = factor(data1@active.ident,
                                  levels(data1@active.ident)[order(clName)])
      #set the names to your IDs; once the gene list works, it is really easy to then update your cluster names with your excel spreadsheet.
      levels(data1@active.ident) = clusterNames$ID[clName][order(clName)]
      #save the IDs in metadata:
      data1@meta.data$IDs = data1@active.ident
      
      
      DotPlot(data1.subsets$`retractor muscle`,'RNA',unique(c('SoxC','AshC','INSM1-like-1','AshA',
                                                              'NKx2.2D','ZNF845','Nem64','Nem24','Hand1','Hand2')),
                                                              # goi[sort(clName)])),
              cols=gene.cp[c(2,11)],scale.by = 'size',col.min = 0
      )&RotatedAxis()
      
      #Then update markers with 
      j=11
    }
    
  }
  DimPlot(data1,cols=clust.cp.separate)
  data1.subsets$`retractor muscle`=data1
  
  ## cluster Ident pharynx ----
  data1=data1.subsets$pharyngeal.ect
  data1<-SetIdent(data1,value = 'IDs.coarse')
  
  DimPlot(data1,cols=clust.cp.separate)
  {
    clusterNames<- read_excel("SupplementalData2.xlsx",sheet = 'AllData.Pharyngeal')
    goi = clusterNames$Marker
    
    DimPlot(data1,label=T,cols = clust.cp.separate)+NoAxes()
    #assign cluster ID to the individual libraries
    data1<-ScaleData(data1,features = goi, split.by = 'orig.ident') #works with scaled gene values, so must scale all of these first (in case they are not within the VarFeatures set!)
    
    #set empty name variable:
    clName = vector()
    
    #generate a matrix of values of each cluster for each gene:
    m=AverageExpression(data1,features=goi,slot='scale.data') #Seurat has a function for this
    
    for (j in 1:length(levels(data1@active.ident))) #for each cluster set
    {
      clName[j]=as.integer(which.max(m$RNA[,j])) #choose the highest value
    }
    sort(clName) #this prints which IDs were selected; sanity check that it did what you wanted... should have no unwanted replicated IDs
    clusterNames$ID[clName] #you can also look at the names to see where the problem might have been
    
    
    #first order the identities..
    data1@active.ident = factor(data1@active.ident,
                                levels(data1@active.ident)[order(clName)])
    #set the names to your IDs; once the gene list works, it is really easy to then update your cluster names with your excel spreadsheet.
    levels(data1@active.ident) = clusterNames$ID[clName][order(clName)]
    #save the IDs in metadata:
    data1@meta.data$IDs = data1@active.ident
    data1.subsets$pharyngeal.ect=data1
    
    
    
  }
  DimPlot(data1,cols=clust.cp.separate)
  
  ## cluster Ident gastrodermis ----
  data1=data1.subsets$gastrodermis
  data1<-SetIdent(data1,value = 'IDs.coarse')
  
  DimPlot(data1,cols=clust.cp.separate)
  {
    clusterNames<- read_excel("SupplementalData2.xlsx",sheet = 'AllData.Gastrodermis')
    goi = clusterNames$Marker
    
    DimPlot(data1,label=T,cols = clust.cp.separate)+NoAxes()
    #assign cluster ID to the individual libraries
    data1<-ScaleData(data1,features = goi, split.by = 'orig.ident') #works with scaled gene values, so must scale all of these first (in case they are not within the VarFeatures set!)
    
    #set empty name variable:
    clName = vector()
    
    #generate a matrix of values of each cluster for each gene:
    m=AverageExpression(data1,features=goi,slot='scale.data') #Seurat has a function for this
    
    for (j in 1:length(levels(data1@active.ident))) #for each cluster set
    {
      clName[j]=as.integer(which.max(m$RNA[,j])) #choose the highest value
    }
    sort(clName) #this prints which IDs were selected; sanity check that it did what you wanted... should have no unwanted replicated IDs
    clusterNames$ID[clName] #you can also look at the names to see where the problem might have been
    
    
    #first order the identities..
    data1@active.ident = factor(data1@active.ident,
                                levels(data1@active.ident)[order(clName)])
    #set the names to your IDs; once the gene list works, it is really easy to then update your cluster names with your excel spreadsheet.
    levels(data1@active.ident) = clusterNames$ID[clName][order(clName)]
    #save the IDs in metadata:
    data1@meta.data$IDs = data1@active.ident
  }
  DimPlot(data1,cols=clust.cp.separate)
  data1.subsets$gastrodermis=data1
  
  ## cluster Ident immune ----
  data1=data1.subsets$unchar.immune
  data1<-SetIdent(data1,value = 'IDs.coarse')
  # not used...
  DimPlot(data1,cols=clust.cp.separate)
  {
    clusterNames<- read_excel("SupplementalData2.xlsx",sheet = 'Alldata.immune.Nv2')
    goi = clusterNames$Marker
    
    DimPlot(data1,label=T,cols = clust.cp.separate)+NoAxes()&DimPlot(data1,label=F,cols = LibCP,group.by='orig.ident')+NoAxes()
    #assign cluster ID to the individual libraries
    data1<-ScaleData(data1,features = goi, split.by = 'orig.ident') #works with scaled gene values, so must scale all of these first (in case they are not within the VarFeatures set!)
    
    #set empty name variable:
    clName = vector()
    
    #generate a matrix of values of each cluster for each gene:
    m=AverageExpression(data1,features=goi,slot='scale.data') #Seurat has a function for this
    
    for (j in 1:length(levels(data1@active.ident))) #for each cluster set
    {
      clName[j]=as.integer(which.max(m$RNA[,j])) #choose the highest value
    }
    sort(clName) #this prints which IDs were selected; sanity check that it did what you wanted... should have no unwanted replicated IDs
    clusterNames$ID[clName] #you can also look at the names to see where the problem might have been
    
    
    #first order the identities..
    data1@active.ident = factor(data1@active.ident,
                                levels(data1@active.ident)[order(clName)])
    #set the names to your IDs; once the gene list works, it is really easy to then update your cluster names with your excel spreadsheet.
    levels(data1@active.ident) = clusterNames$ID[clName][order(clName)]
    #save the IDs in metadata:
    data1@meta.data$IDs = data1@active.ident
    
  }
  DimPlot(data1,cols=clust.cp.separate)
  data1.subsets$unchar.immune=data1
  
  if (!exists("alldata.markers"))
    load(file='AlldataDEGS.Robj')
  
  # generate comparative gene lists:
  immune.deg.ect.all=alldata.markers$gene[alldata.markers$cluster=='immune']
  immune.deg.ect=markers$neurogland.all$gene[markers$neurogland.all$cluster=='immune']
  immune.deg.end=all.markers$gene[all.markers$cluster=='immune']
  a1=(immune.deg.ect.all)
  a2=(immune.deg.ect)
  a3=(immune.deg.end)
  plot.new()
  # FigS2.2 center
VennDiagram::draw.triple.venn(length(a1),length(a2),length(a3),
                                 length(intersect(a1,a2)),length(intersect(a1,a3)),
                                 length(intersect(a2,a3)), 
                                 length(intersect(a1,(intersect(a2,a3)))),
                                 category = c('ect.all','ect.neur','end'),
                                 fill =c("darkred","steelblue3","darkgreen"),
                                 cex=2)
  DotPlot(Alldata.Nv2,'RNA',unique(c(a1,a2,a3)),group.by = 'ID.separate')&RotatedAxis()
  
  Yehulist=genes$gene_short_name[match(c('NVE4535','NVE7179','NVE8855','NVE21594','NVE7180','NVE12407','NVE695','NVE14728','NVE10479','NVE21992','NVE19650','NVE11443','NVE16835','NVE10529'),genes$NVE,nomatch = 0)]
  
  
  ## cluster Ident cnidocyte ----
  data1=data1.subsets$cnidocyte
  data1<-SetIdent(data1,value = 'IDs.coarse')

  DimPlot(data1,cols=clust.cp.separate)+
  DimPlot(data1,group.by = 'orig.ident',cols=LibCP)
  {
    clusterNames<- read_excel("SupplementalData2.xlsx",sheet = 'Alldata.cnidocytes.Nv2')
    goi = clusterNames$Marker
    
    DimPlot(data1,label=T,cols = clust.cp.separate)+NoAxes()
    #assign cluster ID to the individual libraries
    data1<-ScaleData(data1,features = goi, split.by = 'orig.ident') #works with scaled gene values, so must scale all of these first (in case they are not within the VarFeatures set!)
    
    #set empty name variable:
    clName = vector()
    
    #generate a matrix of values of each cluster for each gene:
    m=AverageExpression(data1,features=goi,slot='scale.data') #Seurat has a function for this
    
    for (j in 1:length(levels(data1@active.ident))) #for each cluster set
    {
      clName[j]=as.integer(which.max(m$RNA[,j])) #choose the highest value
    }
    sort(clName) #this prints which IDs were selected; sanity check that it did what you wanted... should have no unwanted replicated IDs
    clusterNames$ID[clName] #you can also look at the names to see where the problem might have been
    
    
    #first order the identities..
    data1@active.ident = factor(data1@active.ident,
                                levels(data1@active.ident)[order(clName)])
    #set the names to your IDs; once the gene list works, it is really easy to then update your cluster names with your excel spreadsheet.
    levels(data1@active.ident) = clusterNames$ID[clName][order(clName)]
    #save the IDs in metadata:
    data1@meta.data$IDs = data1@active.ident
  }
  
  DimPlot(data1,cols=clust.cp.separate)
  data1.subsets$cnidocyte=data1
  
  # Neural and Gland altogether ----

  data1=data1.subsets$neurogland.all
  re-do=F
  if(re-do)
  {
    data1@active.assay = 'RNA'
    data1 <- FindVariableFeatures(data1,nfeatures = 3000)#,
    
    data1 <- ScaleData(data1,split.by = 'orig.ident')
    
    # run different reduction algorythms:
    # #PCA
    data1 <- RunPCA(data1, pcs.compute = 50)
    ElbowPlot(object = data1, ndims = 50)
    
    # set dimensions
    d=as.integer(which(data1@reductions$pca@stdev>2))
    if (length(d) <10 || length(d)>=30)
      d=1:30
    
    data1 <- RunUMAP(data1, reduction ="pca",
                     n.neighbors = 5L,spread =1,
                     dims = d,reduction.name ='umap',
                     reduction.key ='umap',min.dist = 0.4,
                     local.connectivity = 1)#
    data1 <- FindNeighbors(object = data1,
                           reduction ="pca",dims = d,
                           nn.method = 'annoy',
                           annoy.metric = 'cosine',
                           k.param = 10)
    data1 <- FindClusters(object = data1,resolution = 1,random.seed = 0)
    # #look at relationship between clusters
    data1 <- BuildClusterTree(object = data1, reorder = TRUE,
                              reorder.numeric = TRUE)#, dims = c(1:d))
    data1$IDs.fine = data1@active.ident
    data1 <- FindClusters(object = data1,resolution = 0.2,random.seed = 0)
    # #look at relationship between clusters
    data1 <- BuildClusterTree(object = data1, reorder = TRUE,
                              reorder.numeric = TRUE)#, dims = c(1:d))
    data1$IDs.coarse = data1@active.ident
  } 
  {
    clusterNames<- read_excel("SupplementalData2.xlsx",sheet = 'NeurGlanAll.Nv2')
    goi = clusterNames$Marker
    goi[18]='GLWa'
    data1 <- SetIdent(data1, value = 'IDs.fine') #if cluster names are not integers
    data1 <- BuildClusterTree(data1,reorder = T,reorder.numeric = T) #this re-sets everything to numeric; required by script below
    # DotPlot(data1,features= goi)+RotatedAxis()+NoLegend() #checks that all your genes are available and actually express specifically; useful for building the gene list
    DimPlot(data1,label=T,cols = clust.cp.separate)+NoAxes()
    #assign cluster ID to the individual libraries
    data1<-ScaleData(data1,features = goi, split.by = 'orig.ident') #works with scaled gene values, so must scale all of these first (in case they are not within the VarFeatures set!)
    
    #set empty name variable:
    clName = vector()
    
    #generate a matrix of values of each cluster for each gene:
    m=AverageExpression(data1,features=goi,slot='scale.data') #Seurat has a function for this
    
    for (j in 1:length(levels(data1@active.ident))) #for each cluster set
    {
      clName[j]=as.integer(which.max(m$RNA[,j])) #choose the highest value
    }
    sort(clName) #this prints which IDs were selected; sanity check that it did what you wanted... should have no unwanted replicated IDs
    clusterNames$ID[clName] #you can also look at the names to see where the problem might have been
    
    #first order the identities..
    data1@active.ident = factor(data1@active.ident,
                                levels(data1@active.ident)[order(clName)])
    #set the names to your IDs; once the gene list works, it is really easy to then update your cluster names with your excel spreadsheet.
    levels(data1@active.ident) = clusterNames$ID[clName][order(clName)]
    #save the IDs in metadata:
    data1@meta.data$IDs = data1@active.ident
    data1.subsets$neurogland.all = data1
  }
  DimPlot(data1,cols=c(clust.cp.separate,clust.cp.graded))
  data1$orig.ident = droplevels(as.factor(data1$orig.ident))
  ids.cluster.library = as.data.frame(table(Idents(data1), data1@meta.data$orig.ident))
  colnames(ids.cluster.library) = c('ID','Library','CellCount')
  
  dist.clust=
    ggplot(ids.cluster.library, aes(fill=ID, y= CellCount,
                                    x=Library)) +
    geom_bar(mapping =aes(fill=ID, y= (CellCount),
                          x=(Library)),
             position="fill", stat="identity", width = 0.5)+
    scale_fill_manual(values = clust.cp.separate)+
    theme(axis.text.x = element_text(colour=LibCP,face="bold", 
                                     size=8, angle=-45,hjust=0,vjust = 0.5))+
    geom_area(mapping =aes(fill=ID, y= (CellCount),
                           x=as.integer(Library)),
              position="fill", stat="identity",alpha=0.2 , size=.5, colour="white") +
    geom_bar(mapping =aes(fill=ID, y= (CellCount),#this re-plots the bars over the area
                          x=(Library)),
             position="fill", stat="identity", width = 0.5)+
    labs(title="B   Distribution of cell types in time and space",subtitle='relative contribution')+
    theme(legend.position = 'right',legend.title = element_text(size = 6),legend.text = element_text(size=6))
  DimPlot(data1,cols=clust.cp.separate)+NoLegend()+NoAxes()+labs(title='A')+dist.clust
  
  {  if (!exists("markers"))
    load (file = "DEG.subsets.RData")
    
    all.markers=markers$neurogland.all
    list = NULL
    for (i in 1:length(levels(data1@active.ident)))
    {
      x = all.markers[as.numeric(all.markers$cluster) == i, ][1:min(3, length(which(as.numeric(all.markers$cluster) == i))), 7]
      
      list = c(list, x)
    }
    
    #Image the list
    DEG.plot = DotPlot(
      data1,'RNA',
      features = unique(c(list)),
      scale.by = 'radius' ,
      col.min = 0,
      col.max = 3,
      cols = c('lightgrey', 'red')
    ) +
      RotatedAxis() + FontSize(7, 7) + #coord_flip()+
      # NoLegend()+
      labs(title='C', subtitle = 'Top 5 DEGs')+theme(legend.position = 'bottom')
    
    #also for transcription factors only:
    all.markers_TF = markers.TF$neurogland.all
    
    list_TF = NULL
    for (i in 1:length(levels(data1@active.ident)))
    {
      x = all.markers_TF[as.numeric(all.markers_TF$cluster) == i, ][1:min(5, length(which(as.numeric(all.markers_TF$cluster) == i))), 7]
      list_TF = c(list_TF, x)
    }
    list_TF = unique(list_TF)
    TF.plot = DotPlot(
      data1,
      'RNA',
      features = unique(c('SoxC','SoxB.2',list_TF)),
      scale.by = 'radius',
      col.min = 0,
      col.max = 3,
      split.by = 'lifehistory',
      cols = c('slateblue4', 'darkorange')
    ) +
      RotatedAxis() + FontSize(7, 7) +theme(legend.position = 'bottom')+
      labs(title='D',subtitle = 'Top 5 TFs')
  }
  DEG.plot+TF.plot+plot_layout(ncol=1)
  
# Epithelia: In and Out ----
  
  ## Gastrodermis----
  {
  data1=subset(Alldata.Nv2,idents=levels(Alldata.Nv2)[11:12])
  data1@active.assay = 'RNA'
  data1 <- FindVariableFeatures(data1,nfeatures = 2000)#,
  
  
  coi = as.vector(data1@assays$RNA@counts@Dimnames[[2]])
  t=ScaleData(Alldata.Nv2,
              split.by = 'orig.ident', features = data1@assays$RNA@var.features)
  data1@assays$RNA@scale.data = t@assays$RNA@scale.data[,coi]
 
  # run different reduction algorythms:
  # #PCA
  data1 <- RunPCA(data1, pcs.compute = 50)
  ElbowPlot(object = data1, ndims = 50)
  
  # set dimensions
  d=as.integer(which(data1@reductions$pca@stdev>2))
  if (length(d) <10 || length(d)>=30)
    d=1:20
  
  data1 <- RunUMAP(data1, reduction ="pca", 
                   n.neighbors = 10L,spread =1,
                   dims = d,reduction.name ='umap',
                   reduction.key ='umap',min.dist = 0.3,
                   local.connectivity = 10)#
  data1 <- FindNeighbors(object = data1,
                         reduction ="pca",dims = d,
                         nn.method = 'annoy',  
                         annoy.metric = 'cosine',
                         k.param = 20)
  data1 <- FindClusters(object = data1,resolution = 0.6,random.seed = 0)
  # #look at relationship between clusters
  data1 <- BuildClusterTree(object = data1, reorder = TRUE,
                            reorder.numeric = TRUE)#, dims = c(1:d))
  data1$IDs.fine = data1@active.ident
  data1 <- FindClusters(object = data1,resolution = 0.2,random.seed = 0)
  # #look at relationship between clusters
  data1 <- BuildClusterTree(object = data1, reorder = TRUE,
                            reorder.numeric = TRUE)#, dims = c(1:d))
  data1$IDs.coarse = data1@active.ident
  data1@meta.data$orig.ident = as.factor(data1@meta.data$orig.ident)
  data1@meta.data$orig.ident = factor(data1@meta.data$orig.ident,levels(data1@meta.data$orig.ident)[c(2,4,3,5:13,1,14:21)])
  }
  
  {
    data1=SetIdent(data1,value='IDs.fine')
    clusterNames<- read_excel("ClusterAnnotations_Nv2.vs1.xlsx",sheet = 'Alldata.Gastrodermis.merged (2)')
    goi = clusterNames$Marker
    
     DimPlot(data1,label=T,cols = clust.cp.separate)+NoAxes()#&DimPlot(data1,cols = LibCP,group.by = 'orig.ident')+NoAxes()
    # #assign cluster ID to the individual libraries
    data1<-ScaleData(data1,features = goi, split.by = 'orig.ident') #works with scaled gene values, so must scale all of these first (in case they are not within the VarFeatures set!)
    
    #set empty name variable:
    clName = vector()
    
    #generate a matrix of values of each cluster for each gene:
    m=AverageExpression(data1,features=goi,slot='scale.data') #Seurat has a function for this
    
    for (j in 1:length(levels(data1@active.ident))) #for each cluster set
    {
      clName[j]=as.integer(which.max(m$RNA[,j])) #choose the highest value
    }
    sort(clName) #this prints which IDs were selected; sanity check that it did what you wanted... should have no unwanted replicated IDs
    clusterNames$ID[clName] #you can also look at the names to see where the problem might have been
    
    
    #first order the identities..
    data1@active.ident = factor(data1@active.ident,
                                levels(data1@active.ident)[order(clName)])
    #set the names to your IDs; once the gene list works, it is really easy to then update your cluster names with your excel spreadsheet.
    levels(data1@active.ident) = clusterNames$ID[clName][order(clName)]
    #save the IDs in metadata:
    data1@meta.data$IDs = data1@active.ident
  }
  data1.subsets$all.gastrodermis=data1

  endoderm=data1.subsets$all.gastrodermis
    {
      #generate marker lists for each population (cluster)
      endoderm@active.assay='RNA'
      endoderm=FindVariableFeatures(endoderm,nfeatures = 5000)
      all.markers.endoderm <- FindAllMarkers(endoderm,logfc.threshold = 1,
                                             features = endoderm@assays$RNA@var.features,
                                             return.thresh = 0.001,
                                             min.pct = 0.2,
                                             only.pos = TRUE, 
                                             max.cells.per.ident = 200)
    
        
        # add annotations associated with this list:
        
        l=c(names(all.markers.endoderm),colnames(genes)[c(1,4:6,3)])
        all.markers.endoderm[,8:12]<-'NA'
        all.markers.endoderm=setNames(all.markers.endoderm,l) 
        # 
        for (i in 1:length(levels(data1@active.ident))) # 
        {
          x=all.markers.endoderm[as.numeric(all.markers.endoderm$cluster)==i,][1:length(which(as.numeric(all.markers.endoderm$cluster)==i)),7]
          anInd = match(genes[match(x,genes$gene_short_name),1],genes$geneID)
          all.markers.endoderm[as.numeric(all.markers.endoderm$cluster)==i,][1:length(which(as.numeric(all.markers.endoderm$cluster)==i)),8:12]<-genes[anInd,c(1,4:6,3)]
        }  
      
      #generate a collated list of unique DE genes
      list = NULL
      for (i in 1:length(levels(endoderm@active.ident)))
      {
        x=all.markers.endoderm[as.numeric(all.markers.endoderm$cluster)==i,][1:min(5,length(which(as.numeric(all.markers.endoderm$cluster)==i))),7]
        if (is.na (x) ==F)
          list=c(list,x)
      }
      
      #Image the list
      DEG.plot=DotPlot(endoderm, features = unique(c(list)), 
                       scale.by='size' , col.min = 0, col.max = 3, dot.min = 0.1, 
                       cols = c('grey90','red')) + 
        RotatedAxis() +FontSize(10,10) +#coord_flip()+
        NoLegend()+
        labs(title = 'Top 5 DEGs',subtitle = 'p-val < 0.001')+theme(legend.position = 'bottom')
      
      
      #also for transcription factors only:
      all.markers_TF.endoderm <- FindAllMarkers(endoderm,logfc.threshold = 0.6,
                                                features = intersect(TF_list,endoderm@assays$RNA@var.features),
                                                min.pct = 0.05,
                                                only.pos = TRUE, 
                                                max.cells.per.ident = 200,
                                                return.thresh = 0.001)
      
      # add annotations associated with this list:
      
      l=c(names(all.markers_TF.endoderm),colnames(genes)[c(1,4:6,3)])
      all.markers_TF.endoderm[,8:12]<-'NA'
      all.markers_TF.endoderm=setNames(all.markers_TF.endoderm,l) 
      # 
      for (i in 1:length(levels(data1@active.ident))) # 
      {
        x=all.markers_TF.endoderm[as.numeric(all.markers_TF.endoderm$cluster)==i,][1:length(which(as.numeric(all.markers_TF.endoderm$cluster)==i)),7]
        anInd = match(genes[match(x,genes$gene_short_name),1],genes$geneID)
        all.markers_TF.endoderm[as.numeric(all.markers_TF.endoderm$cluster)==i,][1:length(which(as.numeric(all.markers_TF.endoderm$cluster)==i)),8:12]<-genes[anInd,c(1,4:6,3)]
      } 
      
      list_TF = NULL
      for (i in 1:length(levels(endoderm@active.ident)))
      {
        x=all.markers_TF.endoderm[as.numeric(all.markers_TF.endoderm$cluster)==i,][1:min(5,length(which(as.numeric(all.markers_TF.endoderm$cluster)==i))),7]
        if (is.na (x) ==F)
          list_TF=c(list_TF,x)
      }
      list_TF=unique(list_TF)
      TF.plot=DotPlot(endoderm, features = unique(c(list_TF)), 
                      scale.by='size', 
                      col.min = 0, col.max = 3,  dot.min = 0.1,
                      cols = c('grey90','red')) + 
        RotatedAxis() +FontSize(10,10) +NoLegend()+#coord_flip()+
        labs(title = 'Top 5 TFs',subtitle = 'p-val < 0.001')+theme(legend.position = 'bottom')
      DEG.plot+TF.plot+plot_layout(ncol = 1)
      
      
      if(run.save)
      {
        save(all.markers.endoderm,all.markers_TF.endoderm, file = 'endoderm_DEGenes.RData')
      }
    }
    
    list = NULL
    for (i in  1:length(levels(endoderm@active.ident)))
    {
      x=all.markers.endoderm[as.numeric(all.markers.endoderm$cluster)==i,][1:min(50,length(which(as.numeric(all.markers.endoderm$cluster)==i))),7]
      if (is.na (x) ==F)
        list=c(list,x)
    }
    list_TFlong = NULL
    for (i in 1:length(levels(endoderm@active.ident)))
    {
      x=all.markers_TF.endoderm[as.numeric(all.markers_TF.endoderm$cluster)==i,][1:min(50,length(which(as.numeric(all.markers_TF.endoderm$cluster)==i))),7]
      if (is.na (x) ==F)
        list_TFlong=c(list_TFlong,x)
    }
    list_TFlong=unique(list_TFlong)
  
  1## Ectoderm ---- 
  {
    data1=subset(Alldata.Nv2,idents=levels(Alldata.Nv2)[c(9,10,11)])
    data1@meta.data$orig.ident = factor(data1@meta.data$orig.ident,levels(data1@meta.data$orig.ident)[c(1,3,4,2,5:21)])
    
  data1@active.assay = 'RNA'
  data1 <- FindVariableFeatures(data1,nfeatures = 2000)#,
  
  coi = as.vector(data1@assays$RNA@counts@Dimnames[[2]])
  t=ScaleData(Alldata.Nv2,
              split.by = 'orig.ident', features = data1@assays$RNA@var.features)
  data1@assays$RNA@scale.data = t@assays$RNA@scale.data[,coi]
  # run different reduction algorythms:
  # #PCA
  data1 <- RunPCA(data1, pcs.compute = 50)
  ElbowPlot(object = data1, ndims = 50)
  
  # set dimensions
  d=as.integer(which(data1@reductions$pca@stdev>2))
  if (length(d) <=10 || length(d)>=30)
    d=1:20
  
  data1 <- RunUMAP(data1, reduction ="pca", 
                   n.neighbors = 10L,spread =1,
                   dims = d,reduction.name ='umap',
                   reduction.key ='umap',min.dist = 0.1,
                   local.connectivity = 10)#

  data1 <- FindNeighbors(object = data1,
                         reduction ="pca",dims = d,
                         nn.method = 'annoy',  
                         annoy.metric = 'cosine',
                         k.param = 10)
  data1 <- FindClusters(object = data1,resolution = 1,random.seed = 0)
  # #look at relationship between clusters
  data1 <- BuildClusterTree(object = data1, reorder = TRUE,
                            reorder.numeric = TRUE)#, dims = c(1:d))
  data1$IDs.fine = data1@active.ident
  data1 <- FindClusters(object = data1,resolution = 0.045,random.seed = 0)
  # #look at relationship between clusters
  data1 <- BuildClusterTree(object = data1, reorder = TRUE,
                            reorder.numeric = TRUE)#, dims = c(1:d))
  data1$IDs.coarse = data1@active.ident
  data1.subsets$all.ectoderm=data1
  }
  
  {
    data1=SetIdent(data1.subsets$all.ectoderm,value='IDs.coarse')
    clusterNames<- read_excel("ClusterAnnotations_Nv2.vs1.xlsx",sheet = 'all.ectoderm')
    goi = clusterNames$Marker
    DimPlot(data1,label=T,cols = clust.cp.separate)+NoAxes()
    # #assign cluster ID to the individual libraries
    data1<-ScaleData(data1,features = goi, split.by = 'orig.ident') #works with scaled gene values, so must scale all of these first (in case they are not within the VarFeatures set!)
    
    #set empty name variable:
    clName = vector()
    
    #generate a matrix of values of each cluster for each gene:
    m=AverageExpression(data1,features=goi,slot='scale.data') #Seurat has a function for this
    
    for (j in 1:length(levels(data1@active.ident))) #for each cluster set
    {
      clName[j]=as.integer(which.max(m$RNA[,j])) #choose the highest value
    }
    sort(clName) #this prints which IDs were selected; sanity check that it did what you wanted... should have no unwanted replicated IDs
    clusterNames$ID[clName] #you can also look at the names to see where the problem might have been
    
    
    #first order the identities..
    data1@active.ident = factor(data1@active.ident,
                                levels(data1@active.ident)[order(clName)])
    #set the names to your IDs; once the gene list works, it is really easy to then update your cluster names with your excel spreadsheet.
    levels(data1@active.ident) = clusterNames$ID[clName][order(clName)]
    #save the IDs in metadata:
    data1@meta.data$IDs = data1@active.ident
    
  }
  data1.subsets$all.ectoderm=data1

  if (ao.subset)
  {
    ao=subset(data1.subsets$all.ectoderm,cells=WhichCells(data1.subsets$all.ectoderm,idents = levels(data1.subsets$all.ectoderm)[1]))
    data1=ao
    #standard workflow 
    data1 <- NormalizeData(data1, scale.factor = 5000)
    #calculate variable genes
    data1 <- FindVariableFeatures(data1, nfeatures = 2000)

    #use the full dataset scaling:
    coi = as.vector(data1@assays$RNA@counts@Dimnames[[2]])
    t=ScaleData(Alldata.Nv2,split.by = 'orig.ident', features = data1@assays$RNA@var.features)
    data1@assays$RNA@scale.data = t@assays$RNA@scale.data[,coi]
    # data1<-ScaleData(data1,split.by = 'orig.ident')
    data1 <- RunPCA(data1, pcs.compute = 50)
    ElbowPlot(data1,ndims = 50)
    d=10 #20 new
    data1 <- RunUMAP(data1,reduction = 'pca',reduction.name = 'umap',reduction.key = 'umap',
                     n.neighbors = 10L,spread = 1,seed.use = 42,
                     dims = 1:d,min.dist = 0.2,#saved: d=11
                     metric = 'cosine', local.connectivity = 1)#1000 (?!!?)
    o=DimPlot(data1,label=F,group.by = 'IDs')+NoAxes()
    # u=FeaturePlot(data1,'nFeature_RNA')
    l=DimPlot(data1,group.by = 'orig.ident',cols = LibCP,reduction = 'umap')+NoAxes()
    l+o
    
    data1 <- FindNeighbors(object = data1,reduction ="pca",dims = 1:d,
                           nn.method = 'annoy',  
                           annoy.metric = 'cosine', k.param = 10)
    #k=40 for dev
    data1 <- FindClusters(object = data1,resolution = 0.2,random.seed = 0)#RM only 0.05 for end/ectosplit
    data1 <- BuildClusterTree(object = data1, reorder = TRUE, 
                              reorder.numeric = T, dims = c(1:d))
    DimPlot(data1, label = T,label.size = 5, repel = T,order=(levels(data1@active.ident)),
            cols = clust.cp)+NoAxes()+labs(title = 'cluster check')
    levels(data1@active.ident)=c('ect.AO.spot','ect.AO.ring','ect.AO.early')
    data1@active.ident = factor(data1@active.ident,
                                levels(data1@active.ident)[c(3,1,2)])#1,3,4,5,9,10,12,11,13,15,14,2,6,7,8,16
    # run DEG
    ao=data1
    list = NULL
    for (i in  1:length(levels(ao@active.ident)))
    {
      x=all.markers[as.numeric(all.markers$cluster)==i,][1:min(500,length(which(as.numeric(all.markers$cluster)==i))),7]
        list=c(list,x)
    }
    list_TF = NULL
    for (i in 1:length(levels(ao@active.ident)))
    {
      x=all.markers_TF[as.numeric(all.markers_TF$cluster)==i,][1:min(500,length(which(as.numeric(all.markers_TF$cluster)==i))),7]
       list_TF=c(list_TF,x)
    }
   ao.cp=c(stepped(12)[c(9)],stepped3(12)[c(11,9)])
    DimPlot(ao,group.by = 'orig.ident',cols = LibCP)+NoAxes()+NoLegend()+
    DimPlot(ao,cols =ao.cp)+NoAxes()#+NoLegend()
    FeaturePlot(data1,c('POXA-like-1','ISX-like-1'),order=T)&NoAxes()&NoLegend()&scale_color_gradientn(colours = gene.cp)
    save(ao,file = 'ecto.ao.RObj')
    data1=data1.subsets$all.ectoderm
    data2=ao
    cellID = NULL
    cellsp = levels(data2@active.ident)
    
    for (i in 1:length(cellsp))#unique(data2$seurat_clusters)))
      cellID[[i]] = WhichCells(data2, idents = cellsp[i])
    names(cellID)=levels(data2)
    
    #run the plot one on another
    Fig3Ca =  DimPlot(data1,cells.highlight = rev(cellID),
                      cols.highlight = ao.cp,
                      label = F,order = (levels(data2)),pt.size = 1)+NoAxes()+NoLegend()
    Fig3Cb=DimPlot(data2,order = sort(levels(data2@active.ident)),
                   cols=ao.cp,pt=3,
                   #glasbey(length(levels(data2@active.ident)))
    )+NoAxes()
    Fig3Ca+Fig3Cb+plot_layout(ncol =1)
  }
  {
    load(file = 'ecto.ao.RObj')
    data1=ao
    #generate marker lists for each population (cluster)
    data1@active.assay = 'RNA'
    all.markers <- FindAllMarkers(
      data1,
      logfc.threshold = 1,
      features = data1@assays$RNA@var.features,
      #this is faster but restricted to variable gene set.
      return.thresh = 0.001,
      min.pct = 0.2,
      only.pos = TRUE,
      max.cells.per.ident = 200 #downsamples cells in clusters to this size
    )
    # add GO terms associated with this list:
    all.markers$OG.desc <- 'NA'
    all.markers$go.annotation <- 'NA'
    for (i in 1:length(levels(data1@active.ident)))
    {
      x = all.markers[as.numeric(all.markers$cluster) == i, ][1:length(which(as.numeric(all.markers$cluster) == i)), 7]
      anInd = match(x, genes$gene_short_name)
      all.markers[as.numeric(all.markers$cluster) == i, ][1:length(which(as.numeric(all.markers$cluster) ==i)), 8] <- genes$best_OG_desc[anInd]
      all.markers[as.numeric(all.markers$cluster) == i, ][1:length(which(as.numeric(all.markers$cluster) == i)), 9] <- genes$gene_ontology_Pfam[anInd]
    }
    
    #generate a collated list of unique DE genes
    list = NULL
    for (i in 1:length(levels(data1@active.ident)))
    {
      x = all.markers[as.numeric(all.markers$cluster) == i, ][1:min(5, length(which(as.numeric(all.markers$cluster) == i))), 7]
      
      list = c(list, x)
    }
    
    #Image the list
    DEG.plot = DotPlot(
      data1,
      features = unique(c(list)),
      scale.by = 'size' ,
      col.min = 0,
      col.max = 3,
      cols = c('lightgrey', 'darkred')
    ) +
      RotatedAxis() + FontSize(8, 8) + #coord_flip()+
      # NoLegend()+
      labs(title = 'Top 5 DEGs', subtitle = 'p-val < 0.001')
    
    #also for transcription factors only:
    all.markers_TF <- FindAllMarkers(
      data1,
      logfc.threshold = 0.4,
      features = intersect(TF_list,rownames(data1)), #include only TFs that are present
      min.pct = 0.05,
      only.pos = TRUE,
      max.cells.per.ident = 200,
      return.thresh = 0.001,
      grouping.var = 'orig.ident'
    )
    
    # add GO terms associated with this list:
    all.markers_TF$OG_desc <- 'NA'
    all.markers_TF$go.annotation <- 'NA'
    
    for (i in 1:length(levels(data1@active.ident)))
      #
    {
      x = all.markers_TF[as.numeric(all.markers_TF$cluster) == i, ][1:length(which(as.numeric(all.markers_TF$cluster) == i)), 7]
      anInd = match(x, genes$gene_short_name)
      all.markers_TF[as.numeric(all.markers_TF$cluster) == i, ][1:length(which(as.numeric(all.markers_TF$cluster) == i)), 8] <- genes$best_OG_desc[anInd]
      all.markers_TF[as.numeric(all.markers_TF$cluster) == i, ][1:length(which(as.numeric(all.markers_TF$cluster) == i)), 9] <- genes$gene_ontology_Pfam[anInd]
    }
    
    list_TF = NULL
    for (i in 1:length(levels(data1@active.ident)))
    {
      x = all.markers_TF[as.numeric(all.markers_TF$cluster) == i, ][1:min(5, length(which(as.numeric(all.markers_TF$cluster) == i))), 7]
      list_TF = c(list_TF, x)
    }
    list_TF = unique(list_TF)
    TF.plot = DotPlot(
      data1,
      features = unique(c(list_TF)),
      scale.by = 'size',
      col.min = 0,
      col.max = 3,
      cols = c('lightgrey', 'darkred')
    ) +
      RotatedAxis() + FontSize(8, 8) + NoLegend() + #coord_flip()+
      labs(title = 'Top 5 TFs', subtitle = 'p-val < 0.001')
    
    #* **optionally** you can **save these gene lists** to re-call them later:
    # write.csv(all.markers, file = 'XXXXX_DEGenes.csv')
    # write.csv(all.markers_TF, file = 'XXXXX_DEGenesTF.csv')
  }
   ao.markers=all.markers
  ao.tfs=all.markers_TF
  save(ao.markers, ao.tfs,file='ao.markers.RData')
  
  ## PGC & pSC ---- 
  {
# 
    data1=subset(Alldata.Nv2,idents=levels(Alldata.Nv2)[1:2])
    data1@meta.data$orig.ident = factor(data1@meta.data$orig.ident,levels(data1@meta.data$orig.ident)[c(1,3,4,2,5:21)])
    
    data1@active.assay = 'RNA'
    data1 <- FindVariableFeatures(data1,nfeatures = 2000)#,
    
    coi = as.vector(data1@assays$RNA@counts@Dimnames[[2]])
    t=ScaleData(Alldata.Nv2,
                split.by = 'orig.ident', features = data1@assays$RNA@var.features)
    data1@assays$RNA@scale.data = t@assays$RNA@scale.data[,coi]
    # data1 <- ScaleData(data1,split.by = 'orig.ident')
    
    # run different reduction algorythms:
    # #PCA
    data1 <- RunPCA(data1, pcs.compute = 50)
    ElbowPlot(object = data1, ndims = 50)
    
    # set dimensions
    d=as.integer(which(data1@reductions$pca@stdev>2))
    if (length(d) <=10 || length(d)>=30)
      d=1:20
    
    data1 <- RunUMAP(data1, reduction ="pca", 
                     n.neighbors = 10L,spread =0.6,
                     dims = d,reduction.name ='umap',
                     reduction.key ='umap',min.dist = 0.3,
                     local.connectivity = 1)#
    
    data1 <- FindNeighbors(object = data1,
                           reduction ="pca",dims = d,
                           nn.method = 'annoy',  
                           annoy.metric = 'cosine',
                           k.param = 5)
    data1 <- FindClusters(object = data1,resolution = 1,random.seed = 0)
    # #look at relationship between clusters
    data1 <- BuildClusterTree(object = data1, reorder = TRUE,
                              reorder.numeric = TRUE)#, dims = c(1:d))
    data1$IDs.fine = data1@active.ident
    DimPlot(data1,cols=clust.cp.separate,label=T)&NoAxes()
    data1 <- FindClusters(object = data1,resolution = 0.4,random.seed = 0)
    # #look at relationship between clusters
    data1 <- BuildClusterTree(object = data1, reorder = TRUE,
                              reorder.numeric = TRUE)#, dims = c(1:d))
    data1$IDs.coarse = data1@active.ident
  }
  
  FeaturePlot(data1,c('PCNA','SoxC','Neurogenin1','Six1-2'),cols=gene.cp,order=T)&NoAxes()&NoLegend()
  
  
  {
    data1=SetIdent(data1,value='IDs.coarse')
    clusterNames<- read_excel("SupplementalData2.xlsx",sheet = 'Alldata.PGCs_pSC.Nv2')
    goi = clusterNames$Marker
    
     DimPlot(data1,label=T,cols = clust.cp.separate)+NoAxes()#&DimPlot(data1,cols = LibCP,group.by = 'orig.ident')+NoAxes()
    # #assign cluster ID to the individual libraries
    data1<-ScaleData(data1,features = goi, split.by = 'orig.ident') #works with scaled gene values, so must scale all of these first (in case they are not within the VarFeatures set!)
    
    #set empty name variable:
    clName = vector()
    
    #generate a matrix of values of each cluster for each gene:
    m=AverageExpression(data1,features=goi,slot='scale.data') #Seurat has a function for this
    
    for (j in 1:length(levels(data1@active.ident))) #for each cluster set
    {
      clName[j]=as.integer(which.max(m$RNA[,j])) #choose the highest value
    }
    sort(clName) #this prints which IDs were selected; sanity check that it did what you wanted... should have no unwanted replicated IDs
    clusterNames$ID[clName] #you can also look at the names to see where the problem might have been
    
    
    #first order the identities..
    data1@active.ident = factor(data1@active.ident,
                                levels(data1@active.ident)[order(clName)])
    #set the names to your IDs; once the gene list works, it is really easy to then update your cluster names with your excel spreadsheet.
    levels(data1@active.ident) = clusterNames$ID[clName][order(clName)]
    #save the IDs in metadata:
    data1@meta.data$IDs = data1@active.ident
  }
  data1.subsets$pSC.PGC=data1

  #save all ----
  save(data1.subsets, file = 'data1.subsets.updated.Robj')
  save (Alldata.Nv2, file = 'Alldata.Nv2.updated.Robj')  
 