sci_scrna = sci_scrna[, Idents(sci_scrna) %in% c("uninjured","7d")]
sci_scrna$time <- factor(x = sci_scrna$time, levels = c("uninjured","4hr","1d","3d","7d","14d","38d"))
Idents(sci_scrna) <- sci_scrna$time

colour = c("#FB8072", "#1965B0", "#7BAFDE", "#882E72","#B17BA6", 
           "#FDB462", "#E7298A", "#E78AC3","#33A02C", 
           "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D","#E6AB02","#FF7F00")
##
DimPlot(sci_scrna,reduction = "umap",label = F,split.by = "time",group.by = "celltype",cols=colour,pt.size = 0.3)

##
sample_table <- as.data.frame(table(sci_scrna@active.ident,sci_scrna@meta.data$celltype))
names(sample_table) <- c("Time","celltype","CellNumber")

plot_sample<-ggplot(sample_table,aes(x=Time,weight=CellNumber,fill=celltype))+
  geom_bar(position="fill")+
  scale_fill_manual(values=colour) + 
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        axis.line.x = element_line(colour = "black") ,
        axis.line.y = element_line(colour = "black") ,
        plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =16)
  )+labs(y="Percentage")+RotatedAxis()
plot_sample

##
FeaturePlot(sci_scrna,features = c("P2rx7","Ifngr1","Naip6","Nlrc4"),cols = c("lightblue", "red"),pt.size = 0.3,order = T,raster = T)
FeaturePlot(sci_scrna,features = c("P2rx7","Ifngr1","Naip6","Nlrc4"),cols = c("lightblue", "red"),pt.size = 0.3,order = T,raster = F)


##
Idents(sci_scrna) <- sci_scrna$celltype
DotPlot(sci_scrna, features = c("P2rx7","Ifngr1","Naip6","Nlrc4")) + RotatedAxis()

##
library(ggsci)
cors<-pal_igv()(12)
Idents(sci_scrna) <- sci_scrna$time
VlnPlot(sci_scrna,features = c("P2rx7","Ifngr1","Naip6","Nlrc4"),cols = cors,pt.size = 0)

##
#####GSEA analysis####
library(devtools)
library(fgsea)
library(msigdbr)
library(fgsea)
library(dplyr)
library(tidyverse)
library(Seurat)
library(ggplot2)
library(clusterProfiler)
Microglia <- subset(sci_scrna, subset = celltype == "Microglia")
Idents(Microglia) <- Microglia$time
DEGU7d <- FindMarkers(Microglia, ident.1 ="7d", ident.2 = "uninjured", logfc.threshold = 0,min.pct = 0.25)
DEGU7d$genes<-rownames(DEGU7d)

cluster.genes<- DEGU7d %>% arrange(desc(avg_log2FC)) %>% dplyr::select(genes,avg_log2FC) #基因按logFC排序
ranks<- deframe(cluster.genes)
head(ranks)

mdb_c2 <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "BP" )## 定义基因集，选取C2
??msigdbr

m_list = m_df %>% split(x = mdb_c2$gene_symbol, f = mdb_c2$gs_name)

fgsea_sets = mdb_c2 [grep("^GOBP_",mdb_c2$gs_name),] %>% split(x = .$gene_symbol, f = .$gs_name)
length(fgsea_sets)


geneList = cluster.genes$avg_log2FC #把foldchange按照从大到小提取出来
names(geneList) <- cluster.genes$genes #给上面提取的foldchange加上对应上ENTREZID
head(geneList)

#运行fgsea
??GSEA
geneset <- read.gmt("m5.go.bp.v2023.2.Mm.symbols.gmt")  
gsea <- GSEA(geneList, TERM2GENE=geneset,verbose=F)
#运行fgsea
fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)

###条形图
p <-ggplot(fgseaRes %>% as_tibble() %>% arrange(desc(NES)) %>% filter(pval < 0.05) %>% head(n= 20), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES)) +
  coord_flip() +
  labs(x="BP", y="Normalized Enrichment Score",title="") ##输出差异排秩前20的条目
pdf('GSEA-BP.pdf',width=12,height=5)
print(p)
dev.off()



##基线图
library(enrichplot)
geneset_plot <- c("GOBP_POSITIVE_REGULATION_OF_IMMUNE_EFFECTOR_PROCESS",
                  "GOBP_INFLAMMATORY_RESPONSE",
                  "GOBP_PEPTIDE_BIOSYNTHETIC_PROCESS",
                  "GOBP_ACTIVATION_OF_IMMUNE_RESPONSE")


geneset_plot <- c("GOBP_IMMUNE_RESPONSE")
mycol <- pal_nejm()(8)

??gseaplot2
gseaplot2(gsea,  
          geneSetID = geneset_plot,
          color = mycol[c(1:4)],
          title ="Specific GO_BP", 
          rel_heights = c(1.3, 0.3, 0.6),
          ES_geom = "line",
          pvalue_table = FALSE)

####整理数据
fgseaResTidy <- fgseaRes %>%  as_tibble() %>% arrange(desc(NES))
fgseaResTidy %>% dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% arrange(padj) %>% head()
a <- fgseaResTidy[ ,c(1,2,3,4,5,6)]
write.csv(a,"KEGG.csv")



pdf('fgsea_KEGG_PRIMARY_IMMUNODEFICIENCY.pdf',width=8,height=5)
plotEnrichment(fgsea_sets[["KEGG_PRIMARY_IMMUNODEFICIENCY"]],ranks) + labs(title="KEGG_PRIMARY_IMMUNODEFICIENCY") #对某一特定通路分析
dev.off()