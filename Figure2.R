#Figure2_A
library(tidyr)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggrepel)
library(readxl)
Pyroptpsis_geneset <- read_excel("data/Pyroptpsis_geneset.xlsx",sheet = "基因总数")
py_2d <- allDiff_2d %>%
  rownames_to_column("gene_symbol") %>%
  mutate(cluster = rep("2d",times = nrow(allDiff))) %>%
  filter(gene_symbol %in% Pyroptpsis_geneset$all)
py_3d <- allDiff_3d %>%
  rownames_to_column("gene_symbol") %>%
  mutate(cluster = rep("3d",times = nrow(allDiff_3d))) %>%
  filter(gene_symbol %in% Pyroptpsis_geneset$all)
py_35d <- allDiff_35d %>%
  rownames_to_column("gene_symbol") %>%
  mutate(cluster = rep("35d",times = nrow(allDiff_35d))) %>%
  filter(gene_symbol %in% Pyroptpsis_geneset$all)
py_7d <- allDiff_7d %>%
  rownames_to_column("gene_symbol") %>%
  mutate(cluster = rep("7d",times = nrow(allDiff_7d))) %>%
  filter(gene_symbol %in% Pyroptpsis_geneset$all)
py_7d_2 <- allDiff_7d_2 %>%
  rownames_to_column("gene_symbol") %>%
  mutate(cluster = rep("7d",times = nrow(allDiff_7d_2))) %>%
  filter(gene_symbol %in% Pyroptpsis_geneset$all)
allDiff_py <- rbind(py_2d,py_3d,py_7d,py_7d_2,py_35d)
save(allDiff_py,file = "output/allDiff_py.Rdata")
df <- allDiff_py
head(df)
df$gene <- ifelse(df$gene_symbol %in% Pyroptpsis_geneset$all,"py","none")
df$label <- ifelse(df$adj.P.Val<0.05,"adjust P-val<0.05","adjust P-val>=0.05")
head(df)
colnames(df)[1] <- "gene"

top10sig2d_up <- filter(df,cluster=="2d",label=="adjust P-val<0.05",gene=="py") %>% distinct(gene_symbol,.keep_all = T) %>% top_n(10,logFC)
top10sig2d_down <- filter(df,cluster=="2d",label=="adjust P-val<0.05",gene=="py") %>% distinct(gene_symbol,.keep_all = T) %>% top_n(-1,logFC)
top10sig2d <- rbind(top10sig2d_up,top10sig2d_down)

top10sig3d_up <- filter(df,cluster=="3d",label=="adjust P-val<0.05",gene=="py") %>% distinct(gene_symbol,.keep_all = T) %>% top_n(10,logFC)
top10sig3d_down <- filter(df,cluster=="3d",label=="adjust P-val<0.05",gene=="py") %>% distinct(gene_symbol,.keep_all = T) %>% top_n(0,logFC)
top10sig3d <- rbind(top10sig3d_up,top10sig3d_down)

top10sig7d_up <- filter(df,cluster=="7d",label=="adjust P-val<0.05",gene=="py") %>% distinct(gene_symbol,.keep_all = T) %>% top_n(10,logFC)
top10sig7d_down <- filter(df,cluster=="7d",label=="adjust P-val<0.05",gene=="py") %>% distinct(gene_symbol,.keep_all = T) %>% top_n(-1,logFC)
top10sig7d <- rbind(top10sig7d_up,top10sig7d_down)

top10sig35d_up <- filter(df,cluster=="35d",label=="adjust P-val<0.05",gene=="py") %>% distinct(gene_symbol,.keep_all = T) %>% top_n(10,logFC)
top10sig35d_down <- filter(df,cluster=="35d",label=="adjust P-val<0.05",gene=="py") %>% distinct(gene_symbol,.keep_all = T) %>% top_n(-1,logFC)
top10sig35d <- rbind(top10sig35d_up,top10sig35d_down)

top10sig <- rbind(top10sig2d,top10sig3d,top10sig7d,top10sig35d)

df$size <- case_when(!(df$gene_symbol %in% top10sig$gene_symbol)~ 1,
                     df$gene_symbol %in% top10sig$gene_symbol ~ 2)

dt <- filter(df,size==1)

dt <- filter(df,abs(logFC) >= 1)

dt<-dt[complete.cases(dt),]
head(dt)

dt1 <- dt
dt1$cluster <- factor(dt1$cluster,levels = c("2d","3d","7d","35d"))
top10sig$cluster <- factor(top10sig$cluster,levels = c("2d","3d","7d","35d"))

p <- ggplot()+
  geom_jitter(data =dt1,
              aes(x = cluster, y = logFC, color = label),
              size = 1,
              width =0.4)+
  geom_jitter(data = top10sig,
              aes(x = cluster, y = logFC, color = label),
              size = 2,
              width =0.4)
p

dfbar<-data.frame(x=c("2d","3d","7d","35d"),
                  y=c(8.85,7.5,9.2,8.8))
dfbar1<-data.frame(x=c("2d","3d","7d","35d"),
                   y=c(-4.9,-4,-7.2,-4.8))

p1 <- ggplot()+
  geom_col(data = dfbar,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)+
  geom_col(data = dfbar1,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)
p1


dfbar$x<- factor(dfbar$x,levels = c("2d","3d","7d","35d")) 
dfbar1$x<- factor(dfbar1$x,levels = c("2d","3d","7d","35d")) 

p2 <- ggplot()+
  geom_col(data = dfbar,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)+
  geom_col(data = dfbar1,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)+
  geom_jitter(data =dt1,
              aes(x = cluster, y = logFC, color = label,shape=gene),
              size = 1,
              width =0.4)+
  geom_jitter(data = top10sig,
              aes(x = cluster, y = logFC, color = label),
              size = 2,
              width =0.4)
p2

dfcol<-data.frame(x=c(1:4),
                  y=0,
                  label=c(1:4))
mycol <- c("#E64B357F","#4DBBD57F","#00A0877F","#3C54887F")
p3 <- p2 + geom_tile(data = dfcol,
                     aes(x=x,y=y),
                     height=1.9,
                     color = "black",
                     fill = mycol,
                     alpha = 0.6,
                     show.legend = F)
p3

p4 <- p3+
  geom_text_repel(
    data=top10sig,
    aes(x=cluster,y=logFC,label=gene_symbol),
    force = 1.2,
    arrow = arrow(length = unit(0.008, "npc"),
                  type = "open", ends = "last")
  )
p4

p5 <- p4 +
  scale_color_manual(name=NULL,
                     values = c("red","black"))
p5

p6 <- p5+
  labs(x="SCI_TIME",y="average logFC")+
  geom_text(data=dfbar1,
            aes(x=x,y=0,label=x),
            size =5,
            color ="white")
p6

p7 <- p6+
  theme_minimal()+
  theme(
    axis.title = element_text(size = 13,
                              color = "black",
                              face = "bold"),
    axis.line.y = element_line(color = "black",
                               size = 1.2),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.direction = "vertical",
    legend.justification = c(1,0),
    legend.text = element_text(size = 15)
  )
p7

#Figure2_B
library(UpSetR)
library(openxlsx)
library(RColorBrewer)
library(dplyr)
library(tidyverse)
library(tibble)
py_2d <- diffgene_2d %>%
  rownames_to_column("gene_symbol") %>%
  filter(gene_symbol %in% Pyroptpsis_geneset$all)
py_3d <- diffgene_3d %>%
  rownames_to_column("gene_symbol") %>%
  filter(gene_symbol %in% Pyroptpsis_geneset$all)
py_7d <- diffgene_7d %>%
  rownames_to_column("gene_symbol") %>%
  filter(gene_symbol %in% Pyroptpsis_geneset$all)
py_35d <- diffgene_35d %>%
  rownames_to_column("gene_symbol") %>%
  filter(gene_symbol %in% Pyroptpsis_geneset$all)

upsetvenn <- list(SCI_2d = py_2d$gene_symbol,SCI_3d=py_3d$gene_symbol,SCI_7d = py_7d$gene_symbol,SCI_35d = py_35d$gene_symbol)

upset(fromList(upsetvenn), order.by = "freq")

p1 <- upset(fromList(upsetvenn),     
            nsets=length(upsetvenn),#显示数据集的所有数据,nsets = 数值调整可视化数据集数量
            nintersects=30,#显示前多少个
            #sets=c("Control","0min","1h","6h","12h","1d","3d","5d","7d","14d","1m","2m","3m",order.by = "freq",keep.order = TRUE), # 参数可以选择想要展示的左下方的分类变量
            # order.by：调整个分类变量的顺序。freq/degree
            # keep.order：保留sets中给定的顺序。
            #group.by = "sets", #按数据集为单位分组排列
            number.angles = 0, #调整顶部柱状图上数字的角度。
            point.size=4, #调整矩阵中点的大小。
            line.size=1.5, #调整矩阵中点的大小。
            mainbar.y.label="Gene Intersections", #y轴的标签
            main.bar.color = 'black', #y轴柱状图颜色
            matrix.color="black", #x轴点的颜色
            sets.x.label="Number of DEG",   #x轴的标签
            # sets.bar.color=brewer.pal(12,"Set3"),#x轴柱状图的颜色;Set1中只有9个颜色，Set3中有12个颜色，Paired中有12个颜色
            mb.ratio = c(0.7, 0.3), #调整柱状图与矩阵的比例。以一个和为1的向量给出，每个值需要在0.3-0.7的范围内。
            order.by = "freq", #y轴矩阵排序,如"freq"频率，"degree"程度
            text.scale=c(1.5,1.5,1.5,1.5,1.5,1.3), #6个参数intersection size title（y标题大小）,intersection size tick labels（y刻度标签大小）, set size title（set标题大小）, set size tick labels（set刻度标签大小）, set names（set 分类标签大小）, numbers above bars（柱数字大小）的设置
            #shade.color="red",#图中阴影部分的颜色
            #empty.intersections = "on" #展示交集为空的柱状图
            queries=list(list(query=intersects,params=list("SCI_2d","SCI_3d","SCI_35d","SCI_7d"),color="#6C1319",active=T),
                         list(query=intersects,params=list("SCI_7d"),color="#1A2D75",active=T),
                         list(query=intersects,params=list("SCI_2d","SCI_3d","SCI_7d"),color="#15522B",active=T))
)
p1
library(export)
graph2pdf(file="output/venn/upset_venn.pdf")