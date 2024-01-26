#Load & update data -----
library(easypackages)
libraries("Seurat", "Matrix", "readxl","RColorBrewer",'Rmagic',
          'patchwork','dplyr','viridis','ggplot2','pals','SeuratWrappers')
load(file = 'data1.subsets.updated.Robj') #* load *data subsets*
load(file= 'Alldata.Nv2.publish.Robj') #* load *Alldata* object
load (file = 'Genes.Nv2.RData') #* load *workspace* with genes and colour palettes

### update to general sub.cluster names ----
for (j in 1:length(names (data1.subsets)))
{
  data1=data1.subsets[[j]]
  clName=NULL
  for (cl in 1:length(levels(data1)))
    clName[cl] = paste(names(data1.subsets)[j],levels(data1)[cl],sep = '.')
  levels(data1@active.ident) = clName
  # DimPlot(data1,cols = clust.cp.separate)
  data1.subsets[[j]]=data1
}

## Add labels to Alldata ----
{
  cl.order=NULL
  Alldata.Nv2$ID.separate = as.character(Alldata.Nv2$IDs)
  for (i in 1:length(names (data1.subsets)))
  {
    coi=NULL
    coi=colnames(data1.subsets[[i]])
    Alldata.Nv2$ID.separate[coi] = as.character(data1.subsets[[i]]@active.ident[coi])
    cl.order = c(cl.order,levels(data1.subsets[[i]]))
  }
  cl.ind = match(unique(cl.order),levels(as.factor(Alldata.Nv2$ID.separate)))
  # data1@XXXXX= factor(data1@XXXXX,levels(data1@XXXXX)[order])
  Alldata.Nv2$ID.separate = as.factor(Alldata.Nv2$ID.separate)
  Alldata.Nv2$ID.separate = factor(Alldata.Nv2$ID.separate,levels(Alldata.Nv2$ID.separate)[cl.ind])
  DimPlot(Alldata.Nv2,group.by = 'ID.separate',cols=c(clust.cp.graded,clust.cp.separate))&NoLegend()
  levels(Alldata.Nv2$ID.separate)
}

# generate colour palette for clusters ----
clust.cp.all.IDs = NULL
for(i in 1:length(data1.subsets))
  clust.cp.all.IDs = c(clust.cp.all.IDs,rep(clust.cp.separate[i],length(levels(data1.subsets[[i]]))))
names(clust.cp.all.IDs)=cl.order
clust.cp.all.IDs=clust.cp.all.IDs[unique(names(clust.cp.all.IDs))]
pal.bands(clust.cp.all.IDs)


### generate DEG lists ----

markers = NULL
markers.TF = NULL#
for (j in 1:length(names(data1.subsets)))
{
  data1=data1.subsets[[j]]
  data1@active.assay='RNA'

  all.markers <- FindAllMarkers(data1,
                                logfc.threshold = 1,
                                return.thresh = 0.001,
                                min.pct = 0.2,
                                only.pos = TRUE, 
                                max.cells.per.ident = 200)
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

  markers[[j]] = all.markers
  
  #also for transcription factors only: this did not work...
  data1=FindVariableFeatures(data1,nfeatures = 5000)
  
  all.markers_TF <- FindAllMarkers(data1,
                                   logfc.threshold = 0.4,
                                   features = intersect(TF_list,data1@assays$RNA@var.features),
                                   min.pct = 0.05,
                                   only.pos = TRUE, 
                                   max.cells.per.ident = 200,
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

  markers.TF[[j]] = all.markers_TF
}
names(markers) = names(data1.subsets)
names(markers.TF) = names(data1.subsets)

save (markers,markers.TF,file = "DEG.subsets.updated.RData")

#generate a list for ID.separate
{
  data1=SetIdent(Alldata.Nv2,value = 'ID.separate')
  data1@active.assay='RNA'
  all.markers <- FindAllMarkers(data1,
                                logfc.threshold = 1,
                                return.thresh = 0.001,
                                min.pct = 0.2,
                                only.pos = TRUE, 
                                max.cells.per.ident = 200)
  
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
  
  #also for transcription factors only: this did not work...
  data1=FindVariableFeatures(data1,nfeatures = 5000)
  
  all.markers_TF <- FindAllMarkers(data1,
                                   logfc.threshold = 0.4,
                                   features = intersect(TF_list,data1@assays$RNA@var.features),
                                   min.pct = 0.05,
                                   only.pos = TRUE, 
                                   max.cells.per.ident = 200,
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
save (all.markers,all.markers_TF,file='DEG.Alldata.Nv2.ID.separate.updated.RData')

# save DEG list supplement ----
load(file='AlldataDEGS.updated.Robj')
# AllData. partitions:
xlsx::write.xlsx(alldata.markers, file="markers.data1.subsets.updated.xlsx", sheetName='PartitionMarkers', row.names=FALSE, append=T)
xlsx::write.xlsx(alldata.markers.TF, file="markers.data1.subsets.updated.xlsx", sheetName='PartitionDETFs', row.names=FALSE, append=T)
# rm(alldata.markers,alldata.markers.TF )
gc()

# Alldata all clusters:
load(file='DEG.Alldata.Nv2.ID.separate.updated.RData')
xlsx::write.xlsx(all.markers, file="markers.data1.subsets.updated.xlsx", sheetName='AllDEGs', row.names=FALSE, append=T)
xlsx::write.xlsx(all.markers_TF, file="markers.data1.subsets.updated.xlsx", sheetName='AllDETFs', row.names=FALSE, append=T)
rm(all.markers,all.markers_TF)
gc()

load(file='ao.markers.RData')
xlsx::write.xlsx(ao.markers, file="markers.data1.subsets.updated.xlsx", sheetName='ao.DEGs', row.names=FALSE, append=T)
xlsx::write.xlsx(ao.tfs, file="markers.data1.subsets.updated.xlsx", sheetName='ao.DETFs', row.names=FALSE, append=T)
rm(ao.markers,ao.tfs)
gc()

# partition subsets:  
load(file = "DEG.subsets.updated.RData")
for (j in 1:length(names (markers)))
{xlsx::write.xlsx(markers[[j]], file="markers.data1.subsets.updated.xlsx", sheetName=names(markers)[j], row.names=FALSE, append=T)
  gc()
  xlsx::write.xlsx(markers.TF[[j]], file="markers.data1.subsets.updated.xlsx", sheetName=paste(names(markers.TF)[j],"TFs",sep = '.'), row.names=FALSE, append=T)
  gc()}

# Generate Figures Supplement ----
#*Plot Cells on Main Dataset
c=1
fig.allscale=NULL
# fig=NULL
for (c in 1:length(data1.subsets))
{
  data1=data1.subsets[[c]]
  #generate idents table:
  #summarize the data of interest:
  data1$orig.ident = droplevels(as.factor(data1$orig.ident))
  ids.cluster.library = as.data.frame(table(Idents(data1), data1@meta.data$orig.ident))
  colnames(ids.cluster.library) = c('ID','Library','CellCount')
  clust.cp=clust.cp.separate
  #barplot of cluster identities in each library:
  
  dist.clust2=
    ggplot(ids.cluster.library, aes(fill=ID, y= log(CellCount),
                                    x=Library)) +
    geom_bar(mapping =aes(fill=ID, y= log(CellCount),
                          x=(Library)),
             position="stack", stat="identity", width = 0.5)+
    scale_fill_manual(values = clust.cp)+
    theme(axis.text.x = element_text(colour=LibCP,face="bold", 
                                     size=8, angle=-45,hjust=0,vjust = 0.5))+
    geom_area(mapping =aes(fill=ID, y= log(CellCount),
                           x=as.integer(Library)),
              position="stack", stat="identity",alpha=0.2 , size=.5, colour="white") +
    geom_bar(mapping =aes(fill=ID, y= log(CellCount),#this re-plots the bars over the area
                          x=(Library)),
             position="stack", stat="identity", width = 0.5)+
    labs(title="C    Distribution of cell types in time and space", subtitle='absolute cell numbers | log scale')+
    theme(legend.position = 'right',legend.title = element_text(size = 6),legend.text = element_text(size=6))
  
  
  
  p1=DimPlot(Alldata.Nv2,cells.highlight=colnames(data1),cols=c('grey80',clust.cp.separate[c]),raster=T)+NoAxes()+NoLegend()+labs(title='A',subtitle=names(data1.subsets[c]))
  # +
  p2= DimPlot(data1,group.by='orig.ident',cols=LibCP)&NoAxes()&NoLegend()&labs(tag = 'B')
  
  # Plot subset with clustering
  p3=DimPlot(data1,cols=clust.cp)+NoAxes()+NoLegend()+labs(tag = 'D')
  # Plot Top 5 DEG and DETFs from markers file
  #generate a collated list of unique DE genes
  {
    all.markers=markers[[c]]
    all.markers_TF=markers.TF[[c]]
    list = NULL
    for (i in 1:length(levels(data1@active.ident)))
    {
      x=all.markers[as.numeric(all.markers$cluster)==i,][1:min(5,length(which(as.numeric(all.markers$cluster)==i))),7]
      if (is.na (x) ==F)
        list=c(list,x)
    }
    
    #Image the list
    DEG.plot=DotPlot(data1, 'RNA', features = unique(c(alldata.markers$gene[alldata.markers$cluster==names(data1.subsets)[c]][1:10],list)), 
                     scale.by='radius' , dot.min = 0.1,
                     col.min = 0, col.max = 3,
                     cols = c('slateblue3','darkorange3'),split.by = 'lifehistory') + 
      RotatedAxis() +FontSize(8,8) +theme(legend.position = 'bottom',legend.title = element_text(size = 6),legend.text = element_text(size=6),legend.box = 'horizontal')+
      labs(title = 'E',subtitle = 'Top 5 DEGs')
    list_TF = NULL
    for (i in 1:length(levels(data1@active.ident)))
    {
      x=all.markers_TF[as.numeric(all.markers_TF$cluster)==i,][1:min(5,length(which(as.numeric(all.markers_TF$cluster)==i))),7]
      if (is.na (x) ==F)
        list_TF=c(list_TF,x)
    }
    list_TF=unique(list_TF)
    TF.plot=
      DotPlot(data1,'RNA', features = unique(c(list_TF)), 
              scale.by='radius',  dot.min = 0.1,
              col.min = 0, col.max = 3,
              cols = c('slateblue3','darkorange3'),split.by = 'lifehistory') + 
      RotatedAxis() +FontSize(8,8) +theme(legend.position = 'bottom',
                                          # axis.text.y = element_text(colour = rep(clust.cp.separate,each=2)),
                                          legend.title = element_text(size = 6),legend.text = element_text(size=6))+
      labs(title = 'F',subtitle = 'Top 5 DETFs')
  }
  p4=DEG.plot+TF.plot+plot_layout(ncol=1)
  
  p4
  #generate a plot area setup:
  S.layout=c(area(1,1,1,1),area(1,2,1,2),area(2,1,2,1),area(2,2,2,2),area(3,1,4,2),area(5,1,6,2))
  fig.allscale[[c]]=
    p1+p2+dist.clust2+p3+DEG.plot+TF.plot+plot_layout(design = S.layout)
}  

pdf(file='D:/Alison/manuscripts/Nv2Atlas/RDataFigures/supplement.updated.pdf',width=16,height = 24,onefile=T)
for (c in c(9,8,7,12,11,10,6,1,2,3,4,5))
  print(fig.allscale[[c]])
dev.off()  





# Generate RawData Figures ----
### Figure 1 ----
load(file = "DEG.subsets.updated.RData")

##**Fig1a** cell plot
Fig.1a=DimPlot(Alldata.Nv2,label=F,split.by='lifehistory',
               cols = clust.cp.separate)&NoAxes()&NoLegend()

#**fig1.b**
data2=SetIdent(Alldata.Nv2,value='lifehistory')
for (i in 1:2)
{
  data1=subset(data2,idents=levels(data2)[i])
  data1=SetIdent(data1,value='IDs')
  data1<-droplevels(data1)
  data1$orig.ident = droplevels(as.factor(data1$orig.ident))
  ids.cluster.sample = as.data.frame(table(Idents(data1), data1@meta.data$orig.ident))
  colnames(ids.cluster.sample) = c('ID', 'sample', 'CellCount')
  {
    dist.clust2 =
      ggplot(ids.cluster.sample, aes(fill = ID, y = CellCount,
                                     x = sample)) +
      geom_bar(
        mapping = aes(
          fill = ID,
          y = (CellCount),
          x = (sample)
        ),
        position = "fill",
        stat = "identity",
        width = 0.5
      ) +
      scale_fill_manual(values = clust.cp.separate) +
      theme(axis.text.x = element_text(
        #face="bold", color="#993333",
        size = 8,
        angle = -45,
        hjust = 0,
        vjust = 0.5
      )) +
      geom_area(
        mapping = aes(
          fill = ID,
          y = (CellCount),
          x = as.integer(sample)
        ),
        position = "fill",
        stat = "identity",
        alpha = 0.2 ,
        size = .5,
        colour = "white"
      ) +
      geom_bar(
        mapping = aes(
          fill = ID,
          y = (CellCount),
          #this re-plots the bars over the area
          x = (sample)
        ),
        position = "fill",
        stat = "identity",
        width = 0.5
      ) +
      ggtitle("Distribution of cell types in time and space")
  }
  if (i==1)
    Dev.dist=dist.clust2
  else 
    Tissue.dist=dist.clust2
}
Dev.dist+Tissue.dist+NoLegend()+plot_layout(ncol=2)

#**fig 1cd**
{
  load(file='AlldataDEGS.Robj')
  
  Alldata.Nv2<-SetIdent(Alldata.Nv2,value='IDs')
  #generate a collated list of unique DE genes
  list = NULL
  for (i in 1:length(levels(Alldata.Nv2@active.ident)))
  {
    x = alldata.markers[as.numeric(alldata.markers$cluster) == i, ][1:min(5, length(which(as.numeric(alldata.markers$cluster) == i))), 7]
    
    list = c(list, x)
  }
  
  #Image the list
  DEG.plot = DotPlot(
    Alldata.Nv2,
    features = unique(c(list)),
    scale.by = 'size' ,
    col.min = 0,
    col.max = 3,
    cols = c('slateblue3','darkorange3'),split.by = 'lifehistory'
  ) +
    RotatedAxis() + FontSize(8, 8) +theme(legend.position = 'bottom')+ #coord_flip()+
    # NoLegend()+
    labs(title = 'Top 5 DEGs', subtitle = 'p-val < 0.001')
  
  list_TF = NULL
  for (i in 1:length(levels(Alldata.Nv2@active.ident)))
  {
    x = alldata.markers.TF[as.numeric(alldata.markers.TF$cluster) == i, ][1:min(5, length(which(as.numeric(alldata.markers.TF$cluster) == i))), 7]
    list_TF = c(list_TF, x)
  }
  list_TF = unique(list_TF)
  TF.plot = DotPlot(
    Alldata.Nv2,
    features = unique(c(list_TF)),
    scale.by = 'size',
    col.min = 0,
    col.max = 3,
    cols = c('slateblue3','darkorange3'),split.by = 'lifehistory'
  ) +
    RotatedAxis() + FontSize(8, 8) + NoLegend() + #coord_flip()+
    labs(title = 'Top 5 TFs', subtitle = 'p-val < 0.001')+theme(legend.position = 'bottom')
  
  
  DEG.plot+TF.plot+plot_layout(ncol=1)
}

#**Figure 1**
fig1.layout=c(area(1,1,1,2),area(2,1),area(2,2),area(3,1,3,2),area(4,1,4,2))
pdf(file='D:/Alison/manuscripts/Nv2Atlas/RDataFigures/Fig1.raw.pdf',width=16,height = 24,onefile=T)
print(Fig.1a+Dev.dist+Tissue.dist+NoLegend()+DEG.plot+TF.plot+plot_layout(design = fig1.layout))
dev.off()  



### Figure 2 ----
load(file = "DEG.subsets.updated.RData")
load(file = 'ecto.ao.RObj')
load(file='ao.markers.RData')

ectoderm=data1.subsets$all.ectoderm
all.markers.ectoderm=markers$all.ectoderm

ect.clust.cp = c('darkred',stepped(8)[c(5,7,6)],stepped2(20)[c(16,18)],stepped3(4)[c(2,1,3)],'darkblue',stepped(12)[c(12,10,11,9)])
clust.cp = ect.clust.cp

pal.bands(clust.cp)
Fig2A.1=DimPlot(Alldata.Nv2,cells.highlight = ectoderm@assays$RNA@counts@Dimnames[[2]],
                cols.highlight = clust.cp[7])+NoLegend()+NoAxes()#.split$developmental.series
Fig2A=DimPlot(ectoderm,group.by = 'orig.ident',cols = LibCP)+NoAxes()+NoLegend()
ectoderm.ids.cluster.library = as.data.frame(table(Idents(ectoderm), ectoderm@meta.data$orig.ident))
colnames(ectoderm.ids.cluster.library) = c('ID','Library','CellCount')
Fig2E=
  ggplot(ectoderm.ids.cluster.library, aes(fill=ID, y= CellCount,x=Library)) +
  geom_bar(mapping =aes(fill=ID, y=CellCount, x=Library),
           position="fill", stat="identity", width = 0.5)+
  scale_fill_manual(values = clust.cp)+
  theme(axis.text.x = element_text(size=8, angle=-45,hjust=0,vjust = 0.5))+
  geom_area(mapping =aes(fill=ID, y= (CellCount),x=as.integer(Library)),
            position="fill", stat="identity",alpha=0.2 , size=.5, colour="white") +
  geom_bar(mapping =aes(fill=ID, y= (CellCount), x=(Library)),
           position="fill", stat="identity", width = 0.5)+
  ggtitle("Distribution of cell types in time and space")+theme(legend.position = 'bottom')


Fig2.clusterplot = DimPlot(ectoderm,cols = clust.cp)+NoAxes()#+NoLegend()

leg <- ggpubr::get_legend(Fig2E+theme(legend.position = 'right',legend.text.align = 1))
# Convert to a ggplot and print
leg=ggpubr::as_ggplot(leg)
leg
FigB1=FeaturePlot(ectoderm,'FoxA',order=T,cols=gene.cp)+NoAxes()+NoLegend()
FigB2=FeaturePlot(ectoderm,'Wnt2',order=T,cols=gene.cp)+NoAxes()+NoLegend()
FigB3=FeaturePlot(ectoderm,'Six3-6',order=T,cols=gene.cp)+NoAxes()+NoLegend()

#DEGs
{     #generate a collated list of unique DE genes
  list = NULL
  for (i in 1:length(levels(ectoderm@active.ident)))
  {
    x=all.markers.ectoderm[as.numeric(all.markers.ectoderm$cluster)==i,][1:min(5,length(which(as.numeric(all.markers.ectoderm$cluster)==i))),7]
    if (is.na (x) ==F)
      list=c(list,x)
  }
  
  #Image the list
  DEG.plot=DotPlot(ectoderm, 'RNA',features = unique(c(list)), 
                   scale.by='size' , col.min = 0, col.max = 3, dot.min = 0.1, 
                   cols = c('slateblue3','darkorange3'),split.by = 'lifehistory'
  ) + 
    RotatedAxis() +FontSize(0,10) +coord_flip()+
    NoLegend()+
    labs(title = 'Top 5 DEGs',subtitle = 'p-val < 0.001')+theme(legend.position = 'bottom')
  
  list_TF = NULL
  for (i in 1:length(levels(ectoderm@active.ident)))
  {
    x=all.markers_TF.ectoderm[as.numeric(all.markers_TF.ectoderm$cluster)==i,][1:min(5,length(which(as.numeric(all.markers_TF.ectoderm$cluster)==i))),7]
    if (is.na (x) ==F)
      list_TF=c(list_TF,x)
  }
  list_TF=unique(list_TF)
  TF.plot=DotPlot(ectoderm, 'RNA',features = unique(c(list_TF)), 
                  scale.by='size', 
                  col.min = 0, col.max = 3,  dot.min = 0.1,
                  cols = c('slateblue3','darkorange3'),split.by = 'lifehistory') + 
    RotatedAxis() +FontSize(0,10) +NoLegend()+coord_flip()+
    labs(title = 'Top 5 TFs',subtitle = 'p-val < 0.001')+theme(legend.position = 'bottom')
  
}

#ao.degs
{
  {     #generate a collated list of unique DE genes
    list = NULL
    for (i in 1:length(levels(ao@active.ident)))
    {
      x=ao.markers[as.numeric(ao.markers$cluster)==i,][1:min(5,length(which(as.numeric(ao.markers$cluster)==i))),7]
      if (is.na (x) ==F)
        list=c(list,x)
    }
    
    #Image the list
    DEG.ao=DotPlot(ao, 'RNA',features = unique(c(list)), 
                   scale.by='size' , col.min = 0, col.max = 3, dot.min = 0.1, 
                   cols = c('grey90','darkred')) + 
      RotatedAxis() +FontSize(10,10) +#coord_flip()+
      NoLegend()+
      labs(title = 'Top 5 DEGs',subtitle = 'p-val < 0.001')+theme(legend.position = 'bottom')
    
    list_TF = NULL
    for (i in 1:length(levels(ao@active.ident)))
    {
      x=ao.tfs[as.numeric(ao.tfs$cluster)==i,][1:min(5,length(which(as.numeric(ao.tfs$cluster)==i))),7]
      if (is.na (x) ==F)
        list_TF=c(list_TF,x)
    }
    list_TF=unique(list_TF)
    TF.ao=DotPlot(ao, 'RNA',features = unique(c(list_TF)), 
                  scale.by='size', 
                  col.min = 0, col.max = 3,  dot.min = 0.1,
                  cols = c('grey90','darkred')) + 
      RotatedAxis() +FontSize(10,10) +NoLegend()+#coord_flip()+
      labs(title = 'Top 5 TFs',subtitle = 'p-val < 0.001')+theme(legend.position = 'bottom')
    
  }
}

#**Figure 2**
#* needs to be build independently in illustrator; generate the components here

fig2.layout=c(area(1,1),area(1,2),area(1,3),area(1,4),area(2,1,3,2),area(2,3,2,4),area(3,3,3,4),area(4,1,6,1),area(4,2,6,2),area(4,3,6,4))
pdf(file='D:/Alison/manuscripts/Nv2Atlas/RDataFigures/Fig2.raw.pdf',width=20,height = 28,onefile=T)
print(Fig2A.1+FigB1+FigB2+FigB3+
        Fig2A+DEG.ao+TF.ao+
        DEG.plot+TF.plot+Fig2E+
        plot_layout(design = fig2.layout))
dev.off()

# add 
ao.cp=c(stepped(12)[c(9)],stepped3(12)[c(11,9)])
pdf(file='D:/Alison/manuscripts/Nv2Atlas/RDataFigures/Fig2C.raw.pdf',width=15,height = 8,onefile=T)
DimPlot(ao,cols=ao.cp)+NoAxes()+DEG.ao+TF.ao+plot_layout(design = c(area(1,1,2,1),area(1,2),area(2,2)))
dev.off()

### Figure 3 ----
{
  #needs to be build independently in illustrator; generate the components here
endoderm=data1.subsets$all.gastrodermis
  clust.cp = c(rev(stepped2(12))[c(1:3,6:7)],rev(stepped(4)),rev(stepped3(4)),rev(stepped2(4)))
  Fig3Ai=DimPlot(Alldata.Nv2,cells.highlight = endoderm@assays$RNA@counts@Dimnames[[2]],
                 cols.highlight = LibCP[5] )+NoLegend()+NoAxes()#.split$developmental.series
  Fig3A=DimPlot(endoderm,group.by = 'orig.ident',cols = LibCP)+NoAxes()+NoLegend()
  endoderm.ids.cluster.library = as.data.frame(table(Idents(endoderm), endoderm@meta.data$orig.ident))
  colnames(endoderm.ids.cluster.library) = c('ID','Library','CellCount')
  Fig3B=
    ggplot(endoderm.ids.cluster.library, aes(fill=ID, y= CellCount,x=Library)) +
    geom_bar(mapping =aes(fill=ID, y=CellCount, x=Library),
             position="fill", stat="identity", width = 0.5)+
    scale_fill_manual(values = clust.cp)+
    theme(axis.text.x = element_text(size=8, angle=-45,hjust=0,vjust = 0.5))+
    geom_area(mapping =aes(fill=ID, y= (CellCount),x=as.integer(Library)),
              position="fill", stat="identity",alpha=0.2 , size=.5, colour="white") +
    geom_bar(mapping =aes(fill=ID, y= (CellCount), x=(Library)),
             position="fill", stat="identity", width = 0.5)+
    ggtitle("Distribution of cell types in time and space")
  
  Fig3i = DimPlot(endoderm,cols = clust.cp)+NoAxes()+NoLegend()
  names(LibCP)=unique(endoderm$orig.ident)
  pal.bands(LibCP) 
  
  all.markers.endoderm=markers$all.gastrodermis
  {     #generate a collated list of unique DE genes
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
                     cols = c('slateblue3','darkorange3'),split.by = 'lifehistory'
    ) + 
      RotatedAxis() +FontSize(10,10) +#coord_flip()+
      NoLegend()+
      labs(title = 'Top 5 DEGs',subtitle = 'p-val < 0.001')+theme(legend.position = 'bottom')
    all.markers_TF.endoderm =markers.TF$all.gastrodermis
    
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
                    cols = c('slateblue3','darkorange3'),split.by = 'lifehistory') + 
      RotatedAxis() +FontSize(10,10) +NoLegend()+#coord_flip()+
      labs(title = 'Top 5 TFs',subtitle = 'p-val < 0.001')+theme(legend.position = 'bottom')
    
  }
  
  
}
#**Figure 3**
fig3.layout=c(area(1,1),area(1,2),area(2,1),area(2,2),area(4,1,4,2),area(3,1,3,2))
pdf(file='D:/Alison/manuscripts/Nv2Atlas/RDataFigures/Fig3.raw.pdf',width=16,height = 24,onefile=T)
print(Fig3Ai+Fig3i+Fig3A+Fig3B+TF.plot+DEG.plot+NoLegend()+DEG.plot+TF.plot+plot_layout(design = fig3.layout))
dev.off()

### Figure 4 ----

data1=data1.subsets$`retractor muscle`
c=10
#generate idents table:
#summarize the data of interest:
data1$orig.ident = droplevels(as.factor(data1$orig.ident))
ids.cluster.library = as.data.frame(table(Idents(data1), data1@meta.data$orig.ident))
colnames(ids.cluster.library) = c('ID','Library','CellCount')
clust.cp=clust.cp.separate
#barplot of cluster identities in each library:

dist.clust2=
  ggplot(ids.cluster.library, aes(fill=ID, y= log(CellCount),
                                  x=Library)) +
  geom_bar(mapping =aes(fill=ID, y= log(CellCount),
                        x=(Library)),
           position="stack", stat="identity", width = 0.5)+
  scale_fill_manual(values = clust.cp)+
  theme(axis.text.x = element_text(colour=LibCP,face="bold", 
                                   size=8, angle=-45,hjust=0,vjust = 0.5))+
  geom_area(mapping =aes(fill=ID, y= log(CellCount),
                         x=as.integer(Library)),
            position="stack", stat="identity",alpha=0.2 , size=.5, colour="white") +
  geom_bar(mapping =aes(fill=ID, y= log(CellCount),#this re-plots the bars over the area
                        x=(Library)),
           position="stack", stat="identity", width = 0.5)+
  labs(title="B    Distribution of cell types in time and space", subtitle='absolute cell numbers')+
  theme(legend.position = 'NULL',legend.title = element_text(size = 6),legend.text = element_text(size=6))

dist.clust=
  ggplot(ids.cluster.library, aes(fill=ID, y= CellCount,
                                  x=Library)) +
  geom_bar(mapping =aes(fill=ID, y= (CellCount),
                        x=(Library)),
           position="fill", stat="identity", width = 0.5)+
  scale_fill_manual(values = clust.cp)+
  theme(axis.text.x = element_text(colour=LibCP,face="bold", 
                                   size=8, angle=-45,hjust=0,vjust = 0.5))+
  geom_area(mapping =aes(fill=ID, y= (CellCount),
                         x=as.integer(Library)),
            position="fill", stat="identity",alpha=0.2 , size=.5, colour="white") +
  geom_bar(mapping =aes(fill=ID, y= (CellCount),#this re-plots the bars over the area
                        x=(Library)),
           position="fill", stat="identity", width = 0.5)+
  labs(title="XX   Distribution of cell types in time and space",subtitle='relative contribution')+
  theme(legend.position = 'NULL',legend.title = element_text(size = 6),legend.text = element_text(size=6))

# }

p1=DimPlot(Alldata.Nv2,cells.highlight=colnames(data1),cols=c('grey80',clust.cp.separate[c]))+NoAxes()+NoLegend()+labs(title='A`',subtitle=names(data1.subsets[c]))
p3=DimPlot(data1,cols=clust.cp)+NoAxes()+theme(legend.position = 'bottom',legend.text = element_text(size=8))+labs(title='A')
# Plot Top 5 DEG and DETFs from markers file
#generate a collated list of unique DE genes
{
  all.markers=markers[[c]]
  all.markers_TF=markers.TF[[c]]
  list = NULL
  for (i in 1:length(levels(data1@active.ident)))
  {
    x=all.markers[as.numeric(all.markers$cluster)==i,][1:min(5,length(which(as.numeric(all.markers$cluster)==i))),7]
    if (is.na (x) ==F)
      list=c(list,x)
  }
  list_TF = NULL
  for (i in 1:length(levels(data1@active.ident)))
  {
    x=all.markers_TF[as.numeric(all.markers_TF$cluster)==i,][1:min(5,length(which(as.numeric(all.markers_TF$cluster)==i))),7]
    if (is.na (x) ==F)
      list_TF=c(list_TF,x)
  }
  list_TF=unique(list_TF)
  
  
  #Image the lists
  DEG.plot=DotPlot(data1, 'RNA', features = unique(c(alldata.markers$gene[alldata.markers$cluster==names(data1.subsets)[c]][1:10],list,list_TF)), 
                   scale.by='radius' , dot.min = 0.1,
                   col.min = 0, col.max = 3,
                   cols = c('slateblue3','darkorange3'),split.by = 'lifehistory') + 
    RotatedAxis() +FontSize(8,8) +theme(legend.position = 'bottom',axis.text.y = element_text(colour = rep(clust.cp.separate,each=2)),legend.title = element_text(size = 6),legend.text = element_text(size=6))+
    labs(title = 'C',subtitle = 'Top 5 DEGs: Allgene vs Regulatory genes')
  
}
#generate a plot area setup:
#**Figure 4**
fig4.layout=c(area(1,1),area(1,2,2,2),area(2,1),area(3,1,3,2))
pdf(file='D:/Alison/manuscripts/Nv2Atlas/RDataFigures/Fig4.raw.pdf',width=16,height = 24,onefile=T)
print(p1+dist.clust2+p3+DEG.plot+plot_layout(design = fig4.layout))
dev.off()

### Figure 5 ----

{
  data1=data1.subsets$cnidocyte
  #generate idents table:
  #summarize the data of interest:
  data1$orig.ident = droplevels(as.factor(data1$orig.ident))
  ids.cluster.library = as.data.frame(table(Idents(data1), data1@meta.data$orig.ident))
  colnames(ids.cluster.library) = c('ID','Library','CellCount')
  clust.cp=clust.cp.separate
  #barplot of cluster identities in each library:
  
  dist.clust2=
    ggplot(ids.cluster.library, aes(fill=ID, y= log(CellCount),
                                    x=Library)) +
    geom_bar(mapping =aes(fill=ID, y= log(CellCount),
                          x=(Library)),
             position="stack", stat="identity", width = 0.5)+
    scale_fill_manual(values = clust.cp)+
    theme(axis.text.x = element_text(colour=LibCP,face="bold", 
                                     size=8, angle=-45,hjust=0,vjust = 0.5))+
    geom_area(mapping =aes(fill=ID, y= log(CellCount),
                           x=as.integer(Library)),
              position="stack", stat="identity",alpha=0.2 , size=.5, colour="white") +
    geom_bar(mapping =aes(fill=ID, y= log(CellCount),#this re-plots the bars over the area
                          x=(Library)),
             position="stack", stat="identity", width = 0.5)+
    labs(title="C    Distribution of cell types in time and space", subtitle='absolute cell numbers')+
    theme(legend.position = 'right',legend.title = element_text(size = 6),legend.text = element_text(size=6))
  
  dist.clust=
    ggplot(ids.cluster.library, aes(fill=ID, y= CellCount,
                                    x=Library)) +
    geom_bar(mapping =aes(fill=ID, y= (CellCount),
                          x=(Library)),
             position="fill", stat="identity", width = 0.5)+
    scale_fill_manual(values = clust.cp)+
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
  
  # }
  
  p1=DimPlot(Alldata.Nv2,cells.highlight=colnames(data1),cols=c('grey80',clust.cp.separate[c]),raster = T)+NoAxes()+NoLegend()+labs(title='A`',subtitle=names(data1.subsets[c]))
  # +
  # DimPlot(data1,group.by='orig.ident',cols=LibCP)&NoAxes()&NoLegend()
  # p2=pal.bands(LibCP)+coord_flip() #unecessary
  # Plot subset with clustering
  p2=DimPlot(data1,cols=clust.cp)+NoAxes()+NoLegend()+labs(tag='A')
  # Plot Top 5 DEG and DETFs from markers file
  #generate a collated list of unique DE genes
  {
    all.markers=markers$cnidocyte
    all.markers_TF=markers.TF$cnidocyte
    list = NULL
    for (i in 1:length(levels(data1@active.ident)))
    {
      x=all.markers[as.numeric(all.markers$cluster)==i,][1:min(5,length(which(as.numeric(all.markers$cluster)==i))),7]
      if (is.na (x) ==F)
        list=c(list,x)
    }
    
    #Image the list
    DEG.plot=DotPlot(data1, 'RNA', features = unique(c(alldata.markers$gene[alldata.markers$cluster=='cnidocyte'][1:10],list)), 
                     scale.by='radius' , dot.min = 0.1,
                     col.min = 0, col.max = 3,
                     cols = c('slateblue3','darkorange3'),split.by = 'lifehistory') + 
      RotatedAxis() +FontSize(8,8) +theme(legend.position = 'bottom',
                                          # axis.text.y = element_text(colour = rep(clust.cp.separate,each=2)),
                                          legend.title = element_text(size = 6),legend.text = element_text(size=6))+
      labs(tag = 'D',title = 'Up-regulated marker genes',subtitle='p.value <=0.001')
    list_TF = NULL
    for (i in 1:length(levels(data1@active.ident)))
    {
      x=all.markers_TF[as.numeric(all.markers_TF$cluster)==i,][1:min(5,length(which(as.numeric(all.markers_TF$cluster)==i))),7]
      if (is.na (x) ==F)
        list_TF=c(list_TF,x)
    }
    list_TF=unique(list_TF)
    TF.plot=
      DotPlot(data1,'RNA', features = unique(c(list_TF)), 
              scale.by='radius',  dot.min = 0.1,
              col.min = 0, col.max = 3,
              cols = c('slateblue3','darkorange3'),split.by = 'lifehistory') + 
      RotatedAxis() +FontSize(8,8) +theme(legend.position = 'bottom',
                                          # axis.text.y = element_text(colour = rep(clust.cp.separate,each=2)),
                                          legend.title = element_text(size = 6),legend.text = element_text(size=6))+
      labs(title = 'F',subtitle = 'Top 5 DETFs')
    cnido.toxin=c("MAC2", "MAC1", "CO6-like-5", 'CO6-like-4')
    data1<-RunALRA(data1,genes.use = unique(c(cnido.toxin,markers$cnidocyte$gene,markers.TF$cnidocyte$gene)),setDefaultAssay = T)
    pD=FeaturePlot(data1,cnido.toxin,order = T,cols=gene.cp,ncol = 4)&NoAxes()&NoLegend()&labs(tag='C')
  }
  
} 
#**Figure 5**
fig5.layout=c(area(1,1),area(2,1),area(1,2,2,3),area(3,1,3,3),area(4,1,5,3))
pdf(file='D:/Alison/manuscripts/Nv2Atlas/RDataFigures/Fig5.raw.pdf',width=16,height = 24,onefile=T)
print(p1+p2+dist.clust+pD+DEG.plot+plot_layout(design = fig5.layout))
dev.off()

### Figure 6 ----
data1=data1.subsets$pSC.PGC

clust.cp=c('grey','black','orange4','green4','orange', 'pink','steelblue2','blue','blue4')

names(clust.cp)=levels(data1)
pal.bands(rev(clust.cp))
p1=DimPlot(Alldata.Nv2,cells.highlight=colnames(data1),cols=c('grey80',clust.cp.separate[1]),raster = T)+NoAxes()+NoLegend()+labs(tag='A',subtitle='pSC & PGCs')
data1$orig.ident = droplevels(as.factor(data1$orig.ident))
ids.cluster.library = as.data.frame(table(Idents(data1), data1@meta.data$orig.ident))
colnames(ids.cluster.library) = c('ID','Library','CellCount')
dist.clust2=
  ggplot(ids.cluster.library, aes(fill=ID, y= CellCount,
                                  x=Library)) +
  geom_bar(mapping =aes(fill=ID, y= log(CellCount),
                        x=(Library)),
           position="stack", stat="identity", width = 0.5)+
  scale_fill_manual(values = clust.cp)+
  theme(axis.text.x = element_text(colour=LibCP,face="bold", 
                                   size=8, angle=-45,hjust=0,vjust = 0.5))+
  geom_area(mapping =aes(fill=ID, y= log(CellCount),
                         x=as.integer(Library)),
            position="stack", stat="identity",alpha=0.2 , size=.5, colour="white") +
  geom_bar(mapping =aes(fill=ID, y= log(CellCount),#this re-plots the bars over the area
                        x=(Library)),
           position="stack", stat="identity", width = 0.5)+
  labs(tag='C',title="Distribution of cell types in time and space", subtitle='absolute cell numbers')+
  theme(legend.position = 'right',legend.title = element_text(size = 6),legend.text = element_text(size=6))

fig6c=FeaturePlot(data1,c('NUSAP-like-1','PCNA','SoxC','Nanos1'),cols=gene.cp,order=T,raster=T)&NoAxes()&NoLegend()&labs(tag='D')
#,'KI67-like-1','RFA2-like-1'

p3=DimPlot(data1,cols=clust.cp,raster=F,label=T)+NoAxes()+NoLegend()+labs(tag='B',subtitle='clusters')

all.markers=markers$pSC.PGC
all.markers_TF=markers.TF$pSC.PGC

{list = NULL
  for (i in 1:length(levels(data1@active.ident)))
  {
    x=all.markers[as.numeric(all.markers$cluster)==i,][1:min(5,length(which(as.numeric(all.markers$cluster)==i))),7]
    if (is.na (x) ==F)
      list=c(list,x)
  }
  
  #Image the list
  DEG.plot=DotPlot(data1, 'RNA', features = unique(list), 
                   scale.by='radius' , dot.min = 0.1,
                   col.min = 0, col.max = 3,
                   cols = c('slateblue3','darkorange3'),split.by = 'lifehistory') + 
    RotatedAxis() +FontSize(8,8) +theme(legend.position = 'bottom',
                                        axis.text.y = element_text(colour = rep(clust.cp,each=2)),legend.title = element_text(size = 6),legend.text = element_text(size=6))+
    labs(tag = 'E',title='Differentially expressed genes',subtitle = 'p.val <= 0.001')
  list_TF = NULL
  for (i in 1:length(levels(data1@active.ident)))
  {
    x=all.markers_TF[as.numeric(all.markers_TF$cluster)==i,][1:min(5,length(which(as.numeric(all.markers_TF$cluster)==i))),7]
    if (is.na (x) ==F)
      list_TF=c(list_TF,x)
  }
  list_TF=unique(list_TF)
  TF.plot=
    DotPlot(data1,'RNA', features = unique(c(list_TF)), 
            scale.by='radius',  dot.min = 0.1,
            col.min = 0, col.max = 3,
            cols = c('slateblue3','darkorange3'),split.by = 'lifehistory') + 
    RotatedAxis() +FontSize(8,8) +theme(legend.position = 'bottom',
                                        axis.text.y = element_text(colour = rep(clust.cp,each=2)),
                                        legend.title = element_text(size = 6),legend.text = element_text(size=6))+
    labs(tag = 'F',title='Differentially expressed transcription factors',subtitle = 'p.val <= 0.001')
}
DEG.plot+TF.plot+plot_layout(ncol=1)
#generate a plot area setup:
S.layout=c(area(1,1,1,1),area(1,2,1,2),area(2,1,2,1),area(2,2,2,2),area(3,1,3,2),area(4,1,4,2))

pdf(file='D:/Alison/manuscripts/Nv2Atlas/RDataFigures/Fig6.raw.revise.pdf',width=16,height = 26,onefile=T)
p1+p3+dist.clust2+fig6c+DEG.plot+TF.plot+plot_layout(design = S.layout)
dev.off()

###Figure 7 ----

neur.cp=c('blue',"green","darkgreen","#ff7f00","#FFD700" ,"orange2",'#f2db00','#e7a300','orange4','yellow3',"#33a02c",stepped2(8)[6:8],stepped(12)[c(10:12)],'limegreen', "#6A33C2", brewer.purples(5)[c(2:5)],'grey','grey30','grey50',brewer.blues(10)[c(6,9)],"#ff0000",brewer.reds(20)[c(6,15,10)],'orange3',stepped(8)[5:8],'darkred',"#FB6496",brewer.rdpu(7)[2:7],'purple','#f9ffa1',stepped3(8)[7:8])


data1=data1.subsets$neurogland.all
names(neur.cp)=levels(data1)
pal.bands(neur.cp)
DimPlot(data1,cols=neur.cp)
data1$orig.ident = droplevels(as.factor(data1$orig.ident))
ids.cluster.library = as.data.frame(table(Idents(data1), data1@meta.data$orig.ident))
colnames(ids.cluster.library) = c('ID','Library','CellCount')

p1=DimPlot(Alldata.Nv2,cells.highlight = colnames(data1),raster = T)&NoAxes()&NoLegend()
p2=DimPlot(data1,cols=neur.cp)&NoAxes()&NoLegend()&labs(tag='A')
dist.clust=
  ggplot(ids.cluster.library, aes(fill=ID, y= CellCount,
                                  x=Library)) +
  geom_bar(mapping =aes(fill=ID, y= (CellCount),
                        x=(Library)),
           position="fill", stat="identity", width = 0.5)+
  scale_fill_manual(values = neur.cp)+
  theme(axis.text.x = element_text(colour=LibCP,face="bold", 
                                   size=8, angle=-45,hjust=0,vjust = 0.5))+
  geom_area(mapping =aes(fill=ID, y= (CellCount),
                         x=as.integer(Library)),
            position="fill", stat="identity",alpha=0.2 , size=.5, colour="white") +
  geom_bar(mapping =aes(fill=ID, y= (CellCount),#this re-plots the bars over the area
                        x=(Library)),
           position="fill", stat="identity", width = 0.5)+
  labs(tag='B',title="Distribution of cell types in time and space",subtitle='relative contribution')+
  theme(legend.position = 'right',legend.title = element_text(size = 6),legend.text = element_text(size=6))
p2+dist.clust

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
    split.by = 'lifehistory',
    cols = c('slateblue4', 'darkorange')
    # cols = c('lightgrey', 'red')
  ) +
    RotatedAxis() + FontSize(6, 7) + #coord_flip()+
    # NoLegend()+
    labs(tag='C', title = 'Top 5 differentially expressed genes', subtitle = 'p-val < 0.001')+theme(legend.position = 'bottom')
  
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
    RotatedAxis() + FontSize(6, 7) +theme(legend.position = 'bottom')+
    labs(tag='D',title = 'Top 5 upregulated TFs', subtitle = 'p-val < 0.001')
}
DEG.plot+TF.plot+plot_layout(ncol=1)


#generate a plot area setup:
Fig7.layout=c(area(1,1),area(1,2,1,3),area(1,4,1,5),area(2,1,3,5),area(4,1,5,5))

pdf(file='D:/Alison/manuscripts/Nv2Atlas/RDataFigures/Fig7.raw.pdf',width=18,height = 26,onefile=T)
p1+p2+dist.clust+DEG.plot+TF.plot+plot_layout(design = Fig7.layout)
dev.off()

###Figure 8 ----
pdf(file='D:/Alison/manuscripts/Nv2Atlas/RDataFigures/Fig8_raw.pdf',width=16,height = 24,onefile=T)
  print(fig.allscale[[8]])
dev.off()  

###Figure 9 ----
{
  c=5
  data1=data1.subsets$unchar.immune
  #generate idents table:
  data1$orig.ident = droplevels(as.factor(data1$orig.ident))
  ids.cluster.library = as.data.frame(table(Idents(data1), data1@meta.data$orig.ident))
  colnames(ids.cluster.library) = c('ID','Library','CellCount')
  clust.cp=clust.cp.separate[1:length(levels(data1))]
  names(clust.cp)=levels(data1)
  #barplot of cluster identities in each library:
  
  dist.clust2=
    ggplot(ids.cluster.library, aes(fill=ID, y= log(CellCount),
                                    x=Library)) +
    geom_bar(mapping =aes(fill=ID, y= log(CellCount),
                          x=(Library)),
             position="stack", stat="identity", width = 0.5)+
    scale_fill_manual(values = clust.cp.separate)+
    theme(axis.text.x = element_text(colour=LibCP,face="bold", 
                                     size=8, angle=-45,hjust=0,vjust = 0.5))+
    geom_area(mapping =aes(fill=ID, y= log(CellCount),
                           x=as.integer(Library)),
              position="stack", stat="identity",alpha=0.2 , size=.5, colour="white") +
    geom_bar(mapping =aes(fill=ID, y= log(CellCount),#this re-plots the bars over the area
                          x=(Library)),
             position="stack", stat="identity", width = 0.5)+
    labs(tag='C',title="Distribution of cell types in time and space", subtitle='absolute cell numbers | log scale')+
    theme(legend.position = 'right',legend.title = element_text(size = 6),legend.text = element_text(size=8))
  
  
  p1=DimPlot(Alldata.Nv2,cells.highlight=colnames(data1),cols=c('grey80',clust.cp.separate[c]))+NoAxes()+NoLegend()+labs(title='A',subtitle=names(data1.subsets[c]))
  # +
  p2= DimPlot(data1,group.by='orig.ident',cols=LibCP)&NoAxes()&NoLegend()&labs(tag = 'B')
  
  # Plot subset with clustering
  p3=DimPlot(data1,cols=clust.cp)+NoAxes()+NoLegend()+labs(tag = 'D')
  # Plot Top 5 DEG and DETFs from markers file
  
  goi=alldata.markers.TF$gene[alldata.markers.TF$cluster=='unchar.immune' & alldata.markers.TF$p_val_adj ==0.000000e+00]
  goi2=markers.TF$all.gastrodermis$gene[markers.TF$all.gastrodermis$cluster=='all.gastrodermis.immune' & markers.TF$all.gastrodermis$p_val_adj <=0.0001]
  goi2=goi2[c(1:3,5)]
  data1=SetIdent(Alldata.Nv2,value='ID.separate')
  E=DotPlot(data1,'RNA',c(goi,goi2),split.by = 'lifehistory',cols=c('slateblue4','darkorange'))&RotatedAxis()&coord_flip()&labs(tag='E')&FontSize(5,10)
  
}  
Fig9.layout=c(area(1,1),area(1,2),area(2,1),area(2,2),area(3,1,3,2))
pdf(file='D:/Alison/manuscripts/Nv2Atlas/RDataFigures/Fig9.raw.pdf',width=18,height = 26,onefile=T)
p1+p2+dist.clust2+p3+E+plot_layout(design = Fig9.layout)
dev.off()
