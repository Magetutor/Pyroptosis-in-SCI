#Figure5_A
t <- result1 %>% t() %>% as.data.frame() %>%
  rownames_to_column("sample") %>%
  gather(key = cell_type,
         value = value, -sample)
head(dt)

dtt <- dt %>%
  group_by(sample) %>%
  mutate(proportion = round(value/sum(value),3))
head(dtt)
dtt <- dtt %>%
  filter(!(cell_type == "CD8 T cells")) %>%
  filter(!(cell_type == "T cells")) %>%
  filter(!(cell_type == "Mast cells")) %>%
  filter(!(cell_type == "Neutrophils")) %>%
  filter(!(cell_type == "Granulocytes")) %>%
  filter(!(cell_type == "B derived")) %>%
  filter(!(cell_type == "Eosinophils")) %>%
  filter(!(cell_type == "Neutrophil")) %>%
  filter(!(cell_type == "Microglia"))

dtt_arrange <- dtt %>%
  group_by(cell_type) %>%
  summarise(de = median(proportion)) %>%
  arrange(desc(de)) %>%
  pull(cell_type)

dtt$cell_type <- factor(dtt$cell_type,levels = unique(dtt_arrange))

dtt$group <- group$group[match(dtt$sample,group$sample)]
dtt$group <- factor(dtt$group,levels = c("sham","7d"))

t <- t_test(group_by(dtt, cell_type), proportion ~ group)
tj <- adjust_pvalue(t, method = 'fdr')
tj

tj <- add_significance(tj, 'p.adj')
tj
lab <- add_xy_position(tj, x = 'cell_type', dodge = 0.65)
p3 <- ggboxplot(dtt, 
                x = "cell_type", 
                y = "proportion",
                fill = "group",
                alpha = 0.8,
                color = "black") +
  scale_fill_manual(values = c("navy","firebrick3")) +
  labs(x = "", y = "proportion") +
  theme_bw() + 
  mytheme + 
  theme(axis.text.x = element_text(angle = 45,size = 11),
        axis.title.y = element_text(size = 11),
        axis.title  = element_text(size = 12)) +
  stat_pvalue_manual(lab, label = 'p.adj.signif', label.size=4, bracket.size=0.5, tip.length = 0.02)
p3
graph2pdf(file="output/boxplot_propotion_group.pdf")

#Figure5_B
annotation <- data.frame("group" = group$group)
rownames(annotation) <- rownames(group)
p <- pheatmap(result,
              show_colnames = F,
              cluster_cols = F,
              scale = "row", 
              annotation_col = annotation,
              color = colorRampPalette(c("#fdebac","white" ,"#582e8c"))(200),
              cellwidth = 25, 
              cellheight = 20,
              fontsize = 12)
graph2pdf(p,file="output/heatmap.pdf")

#Figure5_C
resmcor <- cor(t(result1), method = "pearson")
View(resmcor)
corrplot(resmcor,
         method = "square",
         order = "hclust",
         tl.cex = 0.6,
         tl.col = "black")
resmorp <- cor.mtest(resmcor, method = "pearson",conf.level = 0.95) 
p.mat <- resmorp$p
View(p.mat)
col2 <-  colorRampPalette(c("darkblue", "white", "darkred"))
corrplot(resmcor,
         method = "color",
         order = "hclust",
         tl.cex = 1.1,
         tl.col = "black",
         tl.srt = 45,
         col = col2(1000),
         p.mat = resmorp$p, sig.level = c(.001, .01, .05),outline="white",
         insig = "label_sig",pch.cex = 0.8, pch.col = "white")
graph2pdf(file="output/corrplot.pdf")

#Figure5_D/E/F
result_t <- t(result)
result_t <- data.frame(result_t)
ggscatter(result_t, x = "Pyroptosis", y = "Monocyte",
          size = 1.5,
          add = "reg.line",  # 添加回归线
          add.params = list(color = "#0AA1FF", fill = "#a5dff9",size = 1),  # 自定义回归线的颜色
          conf.int = TRUE) +  # 添加置信区间
  stat_cor(method = "pearson") +
  xlab("Pyroptosis") +
  ylab("Monocyte") +
  ggtitle("The correlation between Pyroptosis and Monocyte") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text = element_text(size = 12),
        axis.title  = element_text(size = 13))