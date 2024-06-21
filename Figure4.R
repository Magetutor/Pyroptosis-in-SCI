#Figure4_A
library(glmnet)
library(tidyr)
mydata <- data.frame(t(gene_merge))
identical(rownames(mydata), rownames(metadata_merge))
mydata <- cbind("group" = metadata_merge$group,mydata)
py <- as.data.frame(py)
py <- py[grep("Pyroptosis",py$Cell_type),]

mydata$group <- ifelse(mydata$group == "7d",1,0)
y <- as.matrix(mydata[,1]) 
x <- as.matrix(mydata[,colnames(mydata) %in% py$Metagene]) 
set.seed(12345)
lasso_model <- glmnet(x,
                      y,
                      family = "binomial",
                      alpha = 1)
print(lasso_model)

plot(lasso_model,
     xvar = "lambda",
     label = F)

coef_lasso <- coef(lasso_model,
                   s = 0.054410) 
coef_lasso

cv_model <- cv.glmnet(x,y,family = "binomial",alpha = 1,nfolds = 10)

plot(cv_model)

lambda_min <- cv_model$lambda.min
lambda_min
lambda_1se <- cv_model$lambda.1se
lambda_1se

coef_cv <- coef(lasso_model,s = lambda_min)
coef_cv
coef_cv <- coef(lasso_model,s = lambda_1se)
coef_cv
exp(coef_cv)

coef_cv <- as.matrix(coef_cv)
coef_cv <- data.frame(coef_cv)

coef_cv$OR <- exp(coef_cv$s1)
nonzero_vars <- rownames(coef_cv[coef_cv$OR != 1,])
nonzero_vars <- nonzero_vars[2:6]

lasso_data <- mydata[,nonzero_vars]
save
set.seed(123)
train_index <- sample(1:nrow(mydata),nrow(mydata)*0.7) 
train_data <- mydata[train_index,]
text_data <- mydata[-train_index,]

#Figure4_B
load("E:/R stuation/py/public_mice/output/7d_2/exprSet_vst.Rdata")
load("E:/R stuation/py/public_mice/output/7d_2/metadata.Rdata")
load("E:/R stuation/py/public_mice/output/7d_2/diffgene_7d_limma.Rdata")
Pyroptpsis_geneset <- read_excel("data/Pyroptpsis_geneset.xlsx",sheet = "基因总数")
mydata <- gene_merge
mydata <- as.data.frame(t(mydata))
group <- metadata_merge$group
mydata <- cbind(group,mydata)
py_7d <- diffgene_7d %>%
  rownames_to_column("gene_symbol") %>%
  filter(gene_symbol %in% Pyroptpsis_geneset$all)

mydata$group <- factor(ifelse(mydata$group == "7d",1,0))

table(mydata$group)
grepl(colnames(mydata),py_7d$gene_symbol)
X <- mydata[,colnames(mydata) %in% py_7d$gene_symbol] 
Y <- as.numeric(as.factor(mydata$group)) 

control <- trainControl(method = "repeatedcv",   
                        number = 5,     
                        repeats = 5,     
                        search = "random"   
)
set.seed(12345)   
svm_rfe <- rfe(X,
               Y,
               sizes = 1:34,
               rfeControl = rfeControl(functions = caretFuncs,
                                       method = "repeatedcv",
                                       number = 5,
                                       repeats = 5,
                                       verbose = FALSE),
               method = "svmLinear",
               trControl = control,
               preProc = c("center", "scale")
)

save(svm_rfe,file = "output/svm_rfe.Rdata")
svm_rfe_ranking <- svm_rfe$variables
head(svm_rfe_ranking)

varImp(svm_rfe)

varImp_dataframe <- data.frame(Gene = row.names(varImp(svm_rfe))[1:34],
                               importance = varImp(svm_rfe)[1:34, 1])

varImp_dataframe <- na.omit(varImp_dataframe)

mycolors <- c('#D4E2A7','#88D7A4','#A136A1','#BAE8BC','#C757AF',
              '#DF9FCE','#D5E1F1','#305691','#B6C2E7','#E8EFF7',
              '#9FDFDF','#EEE0F5','#267336','#98CEDD','#CDE2EE',
              '#DAD490','#372E8A','#4C862D','#81D5B0','#BAE8C9',
              '#A7DCE2','#AFDE9C')

ggplot(varImp_dataframe, aes(x = reorder(Gene, -importance), y = importance , fill = Gene)) + 
  geom_col() +
  ggtitle("Hub Genes") +
  theme(panel.border = element_blank(),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        plot.title = element_text(margin = margin(b = 20)), 
        panel.grid.major = element_line(color = "grey", size = 0.2)) +
  xlab("Gene") + ylab("Importance") +
  scale_fill_manual(values = mycolors)
library(export)
graph2ppt(file="output/Genes_improtance.pptx")
graph2pdf(file="output/Genes_improtance.pdf")

top_10_vars <- svm_rfe_ranking$var[1:10]
top_10_vars

top10_SVM_data <- mydata[,top_10_vars]

X_plot = svm_rfe$results$Variables
Y_plot = svm_rfe$results$RMSE
plot(X_plot, Y_plot,  
     xlab="Variable Number",  
     ylab="RMSE (Cross-Validation)",  
     col="#7DEE44",    
     pch=16,               
     cex=1.5,         
     lwd=2,           
     type="b",
     ylim=c(0.08, 0.17)) 

abline(h=min(Y_plot), col="skyblue")   
grid(col="grey",lwd=1,lty=3) 

legend("topright",c("Training RMSE","Cross-Validation RMSE"),
       col=c("#7DEE44","#DF294C"),pch=c(16,NA),lwd=2,bg="white")  
wmin <- which.min(svm_rfe$results$RMSE)
wmin

points(wmin, svm_rfe$results$RMSE[wmin], col = "orange", pch = 16, cex=2)  
text(wmin, svm_rfe$results$RMSE[wmin],  
     paste0("N=", wmin), pos = 2, col = "orange", cex=2)  

Target_Genes <- svm_rfe$optVariables
Target_Genes

Best_SVM_data <- mydata[,Target_Genes]

save(Best_SVM_data,file = "output/Best_SVM_data.Rdata")

#Figure4_C
library(ggvenn)
p1 <- ggvenn(
  data = xx,         
  columns = NULL,          
  show_elements = F ,       
  label_sep = "\n",         
  show_percentage = T,      
  digits = 0,               
  fill_color = c("#6C1319", "#1A2D75"), 
  fill_alpha = 0.8,         
  stroke_color = "white",   
  stroke_alpha = 0.5,       
  stroke_size = 0.5,         
  stroke_linetype = "solid", 
  set_name_color = "black", 
  set_name_size = 12,        
  text_color = "black",     
  text_size = 10             
)
p1
ggsave(p1, file="venn.pdf", width=12, height=8)

#Figure4_D
spinal_Exp_vst <- data.frame(t(gene_merge))
identical(rownames(spinal_Exp_vst),rownames(metadata_merge))
spinal_Exp_vst <- cbind(metadata_merge$group,spinal_Exp_vst)
colnames(spinal_Exp_vst)[1] <- "group"
library(dplyr)
library(tidyr)
genelist <- c("P2rx7","Ifngr1","Naip6","Nlrc4")
data <- spinal_Exp_vst[,c("group",genelist)]
library(tidyr)
data <- data %>% 
  pivot_longer(cols=-1,
               names_to= "gene",
               values_to = "expression")
data$group <- factor(data$group,levels = "sham","7d")

t <- t_test(group_by(data, gene), expression ~ group)
tj <- adjust_pvalue(t, method = 'fdr') 
tj

tj <- add_significance(tj, 'p.adj')
tj

lab <- add_xy_position(tj, x = 'cell_type', dodge = 0.65)
ggviolin(data, x = "gene", y = "expression", 
         color = "group",
         fill="group",
         palette =c("#4DBBD57F","#E64B357F"),
         add = c("boxplot"),
         add.params = list(color="white"),
         xlab = F,
         legend = "right"
)+
  stat_compare_means(aes(group =group),label = "p.signif")

ggsave(Ncbp3_plot, file="output/sci/ggviolin/m7G/Ncbp3_plot.pdf", width=6, height=8)

