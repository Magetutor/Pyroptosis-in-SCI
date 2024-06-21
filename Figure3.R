#Figure3_A
library(ggplot2)
library(ggrepel)
data <- allDiff_35d
data <- cbind("gene" = rownames(data),data)
logFCfilter = 1.5

index = data$adj.P.Val <0.05 & abs(data$logFC) > logFCfilter
data$group <- 0
data$group[index & data$logFC>0] = 1
data$group[index & data$logFC<0] = -1
data$group <- factor(data$group,levels = c(1,0,-1),labels =c("Up","NS","Down") )

library(ggplot2) 
library(ggprism) 
library(ggrepel) 
ggplot(data, aes(x=logFC, y =-log10(adj.P.Val),color=group))+
  geom_point(alpha=0.85, size=1.5) + 
  scale_color_manual(values=c('brown','gray','steelblue')) + 
  xlim(c(-9, 9)) +  
  geom_vline(xintercept=c(-1.5,1.5),lty=4,col="black",lwd=0.8)+ 
  geom_hline(yintercept = -log10(0.05), lty=4,col="black",lwd=0.8) + 
  labs(x="log2 (fold change)",y="-log10 (adj.P.Val)") +
  ggtitle("2d after SCI") + 
  theme(plot.title = element_text(hjust = 0.5),legend.position="right",legend.title = element_blank())+
  theme_prism(border = T) 

#Figure3_B
library(readxl)
exprSet <- tpm_7d

load("E:/R stuation/py/public_mice/output/7d_2/diffgene_7d_limma.Rdata")
Pyroptpsis_geneset <- read_excel("data/Pyroptpsis_geneset.xlsx",sheet = "基因总数")
library(tidyr)
library(dplyr)
library(tibble)
py_2d <- diffgene_7d %>%
  rownames_to_column("gene_symbol") %>%
  filter(gene_symbol %in% Pyroptpsis_geneset$all)

heatdata <- exprSet[c(py_2d$gene_symbol,"Cycs","Prkaca"),]
heatdata <-na.omit(heatdata)

identical(colnames(heatdata),metadata_merge$sample)
group <- metadata_merge$group
annotation_col <- data.frame(group)
annotation_col$group <- factor(annotation_col$group,levels = c("sham","7d"))
rownames(annotation_col) <- colnames(heatdata)


library(pheatmap)
library(viridisLite)
pheatmap(heatdata)
p=pheatmap(heatdata, 
           cluster_rows = TRUE,
           cluster_cols = T,
           annotation_col =annotation_col,
           annotation_legend=TRUE,
           show_rownames = T,
           scale = "row", 
           color = viridis(10, alpha = 1, begin = 0.5, end = 1, direction = 1),
           cellwidth = 60, 
           cellheight = 10,
           fontsize = 10 
)

#Figure3_C
ggplot(aa,aes(x=group, y=Description),
       label_format = 1) + #多个分组时需要选取Cluster
  geom_point(aes(size = richFactor,color = -log10(pvalue))) + # 气泡大小及颜色设置
  facet_grid(~cluster) +
  labs(x = "Rich Factor",
       y = "Description",
       title = "GO:BP Enrichment Dotplot", # 设置坐标轴标题及图标题
       size = "Rich Factor") +
  theme_bw() + 
  scale_color_distiller(palette = "YlOrBr",direction = 1) + 
  theme (text = element_text (size = 15))


#Figure3_D
library(enrichplot)
gseaNb(object= y,geneSetID="GOBP_MACROPHAGE_ACTIVATION",
       subPlot=3,
       addPval=T,
       pvalX=0.95,
       pvalY=0.8)

#Figure3_F
annotation <- data.frame("group" = group$group)
rownames(annotation) <- rownames(group)
p <- pheatmap(result,
              show_colnames = F,
              scale = "row", #以行来标准化，这个功能很不错
              annotation_col = annotation,
              color = colorRampPalette(c("navy", "white", "firebrick3"))(200),#调色
              cellwidth = 25, # 格子宽度
              cellheight = 20,
              fontsize = 12)
graph2pdf(p,file="output/heatmap.pdf")

#Figure3_G/H/I
ggscatter(result_t, x = "Nfe2l2", y = "Microglia.Differentiation",
          size = 1.5,
          add = "reg.line",  # 添加回归线
          add.params = list(color = "#77C034", fill = "#C5E99B",size = 1),  # 自定义回归线的颜色
          conf.int = TRUE) +  # 添加置信区间
  stat_cor(method = "spearman") +
  xlab("Nfe2l2") +
  ylab("Microglia Differentiation") +
  ggtitle("The correlation between Nfe2l2 and Microglia Differentiation") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text = element_text(size = 12),
        axis.title  = element_text(size = 13))