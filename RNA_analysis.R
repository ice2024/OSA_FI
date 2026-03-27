ExprMatrix=expr.matrix
fenzu=read.table("group.GSE111762.txt",header=T,row.names=1,stringsAsFactors = FALSE,check.names = FALSE)
ExprMatrix=ExprMatrix[,row.names(fenzu)]
all(colnames(ExprMatrix)==row.names(fenzu))
library(limma)
group=factor(fenzu$group)
design <- model.matrix(~-1+group)
colnames(design)<-levels(group)
table(fenzu$group)
contrasts <- makeContrasts(AG-Control,levels=design)
fit<-lmFit(ExprMatrix, design)
fit1<-contrasts.fit(fit, contrasts)
fit2<-eBayes(fit1)
diff.matrix1 <-topTable(fit2,coef=1,number= nrow(ExprMatrix), adjust.method="BH",sort.by="B",resort.by="M")

write.csv(diff.matrix1, file = "diff.matrix.csv")

diff.matrix1=read.csv("diff.matrix.csv",header = T,stringsAsFactors = FALSE,check.names = FALSE,row.names = 1)
s=which(diff.matrix1$P.Value<0.05&abs(diff.matrix1$logFC)>0.5)
s1=diff.matrix1[s,]
write.csv(s1, file = "diff.p0.05.lfc0.5.csv")
devtools::install_github("BioSenior/ggVolcano")
library(ggVolcano)
deg_data=read.csv("GSE56363_diff.matrix.csv",header = T,stringsAsFactors = FALSE,check.names = FALSE,row.names = 1)
colnames(deg_data)
deg_data$label=""
p=gradual_volcano(deg_data, x = "logFC", y = "P.Value",
                x_lab = "log2(FoldChange)",y_lab = "-log10(P.Value)",legend_title = "-log10(P.Value)",
                log2FC_cut = 0.585,FDR_cut = 0.05,
                add_label = TRUE,label = "label",label_number = nrow(deg_data),
                output = TRUE,filename = "volcano_plot"
                #output = FALSE
)+labs(title="NCR vs CR")
ggsave("volcano_plot.pdf",p,width = 10,height = 10)
library(pheatmap)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)
library(dendextend)
library(dendsort)
counts <- read.csv("GSE56363_exp.csv",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
deg <- read.table("GSE56363_diff.p0.05.lfc0.585.csv",sep = ",",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
degree=read.delim("degree.txt",header = T,stringsAsFactors = F,check.names = F)
deg=deg[which(row.names(deg)%in%degree$gene),]
logFC_cutoff=0.5
type1 = (deg$P.Value<0.05)&(deg$logFC < -logFC_cutoff)
type2 = (deg$P.Value<0.05)&(deg$logFC > logFC_cutoff)
deg$change = ifelse(type1,"Down",ifelse(type2,"Up","Not"))
table(deg$change)


deg$Label = ""   
  deg <- deg[order(deg$P.Value), ]   
deg$Gene <- rownames(deg)
up.genes <- head(deg$Gene[which(deg$change == "Up")], 50)

down.genes <- head(deg$Gene[which(deg$change == "Down")], 50)

deg.top5.genes <- c(as.character(up.genes), as.character(down.genes))
deg$Label[match(deg.top5.genes, deg$Gene)] <- deg.top5.genes


deg <- deg[deg.top5.genes,]
deg$id <- rownames(deg)
ann_row=as.data.frame(deg[,c(7:8)])
ann_row <- ann_row[which(ann_row$change!="Not"),]
ann_row <- ann_row[order(ann_row$change,decreasing = T),]
ann_row1 <- as.data.frame(ann_row[,1])
rownames(ann_row1) <- rownames(ann_row)
colnames(ann_row1) <-"Change"

group<- read.delim("group.GSE56363.txt",check.names = F,stringsAsFactors = F,header = T)
ann_col <- group[order(group$group,decreasing = T),]
row.names(ann_col)=ann_col$sample
ann_col$sample=NULL
ann_col$id <- rownames(ann_col)
ann_col1 <- as.data.frame(ann_col[,1])
rownames(ann_col1) <- ann_col$id
colnames(ann_col1) <- "Group"
ha1 = HeatmapAnnotation(df = ann_col1,col = list(Group=c(CR="#0072b5",NCR="#bc3c29")))


data_pheat <- t(scale(t(as.matrix(counts[rownames(ann_row1),rownames(ann_col1)]))))
pdf("heatmap_plot.pdf", width = 6, height = 6)
densityHeatmap(data_pheat,top_annotation = ha1,title = "Distribution as heatmap", ylab = "Values") %v%
  Heatmap(data_pheat, col = colorRampPalette(c("navy", "white", "firebrick3"))(50),name = "Exp", 
          height = unit(8, "cm"),show_row_names = T,show_column_names=F,cluster_rows =T)
dev.off()


library(WGCNA)   
dataExpr=read.delim("after_merge.txt",header=T,row.names=1,stringsAsFactors = FALSE)
group=read.delim("group.wgcna.txt",header=T,row.names=1,stringsAsFactors = FALSE)
dataExpr=dataExpr[,row.names(group)]
all(colnames(dataExpr)==row.names(group))

dataExprVar=dataExpr
diff_exp=dataExprVar
mydata <- as.data.frame(t(dataExprVar))

#mydata=as.data.frame(t(diff_exp))  
nGenes = ncol(mydata)
nSamples = nrow(mydata)
mytree=hclust(dist(mydata), method="complete") 
pdf(file="1.sampleClustering.pdf",width = 20, height = 9)
par(mar = c(5,5,5,5))
plot(mytree, main = "Sample clustering", sub="", xlab="")
dev.off()

powers = seq(1,20)
sft = pickSoftThreshold(mydata, powerVector=powers, verbose=5)  
pdf(file = "2.softPower.pdf", width = 9, height = 5)
par(mfrow = c(1,2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylim = c(-1,1), ylab="Scale Free Topology Model Fit,signed R^2",type="n", main = "Scale independence")
abline(h=0.85,col="blue")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=0.6,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n",main = "Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=0.6,col="red")
dev.off()
softPower=sft$powerEstimate

adjacency=adjacency(mydata,power=softPower)

TOM=TOMsimilarity(adjacency)
save(TOM,file = "TOM.rda")

dissTOM=1-TOM

geneTree = hclust(as.dist(dissTOM), method = "average");
minModuleSize =70

dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize) #建立模块
#table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)

MEList = moduleEigengenes(mydata, colors = dynamicColors)
MEs = MEList$eigengenes 
MEDiss = 1-cor(MEs)
METree = hclust(as.dist(MEDiss), method = "average")

pdf(file = "4.eigengeneClustering.pdf", width = 12, height = 9)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "") 
dev.off()

pdf("5.Eigengene_adjacency_heatmap.pdf")
plotEigengeneNetworks(MEs, 
                      "Eigengene adjacency heatmap", 
                      marHeatmap = c(3,4,2,2), 
                      plotDendrograms = FALSE, 
                      xLabelsAngle = 90)
dev.off()

MEDissThres=0.3
# 合并模块???
merge=mergeCloseModules(mydata,dynamicColors,cutHeight=MEDissThres,verbose=3)
# 合并后的颜色
mergedColors=merge$colors

mergedMEs=merge$newMEs
pdf("6.ModuleTree.pdf", width = 12, height = 9)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()
pdf("7.mergeModuleTree.pdf", width = 12, height = 9)
plotDendroAndColors(geneTree,mergedColors, "Mergeddynamic",dendroLabels=FALSE,hang=0.03,addGuide=TRUE,guideHang=0.05)
dev.off()
moduleColors = mergedColors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder) - 1
MEs = mergedMEs
write.table(paste(colnames(mydata), moduleColors, sep = "\t"), file = "4.netcolor2gene.xls", row.names=FALSE, col.names=FALSE, quote=FALSE)
gene_color=read.table("4.netcolor2gene.xls",row.names=1)
all(row.names(gene_color)==colnames(mydata))
color=unique(gene_color[,1])
color=color[color!="grey"]
traitColors = numbers2colors(group, signed = TRUE,centered=TRUE);
pdf("9.Sample_dendrogram_and_trait_heatmap.pdf",width = 20, height = 9)
plotDendroAndColors(mytree, 
                    traitColors, 
                    groupLabels = names(group), 
                    rowTextAlignment = "right-justified",
                    addTextGuide = TRUE ,
                    hang = 0.03,
                    dendroLabels = NULL, 
                    addGuide = FALSE,  
                    guideHang = 0.05,
                    main = "Sample dendrogram and trait heatmap")
dev.off()
moduleTraitCor_noFP <- cor(mergedMEs, group, use = "p");
moduleTraitPvalue_noFP = corPvalueStudent(moduleTraitCor_noFP, nSamples); 
textMatrix_noFP <-paste(signif(moduleTraitCor_noFP, 2), "\n(", signif(moduleTraitPvalue_noFP, 1), ")", sep = ""); 
dim(textMatrix_noFP) = dim(moduleTraitCor_noFP)
pdf("10.Module-trait_relationships.pdf",width = 6, height = 8)
par(mar = c(4, 6, 5, 3)); 
labeledHeatmap(Matrix = moduleTraitCor_noFP, 
               xLabels = colnames(group), 
               yLabels = names(mergedMEs), 
               ySymbols = names(mergedMEs), 
               colorLabels = FALSE, 
               colors = blueWhiteRed(50), 
               textMatrix = textMatrix_noFP,
               setStdMargins = FALSE, 
               cex.text = 0.65, 
               zlim = c(-1,1), 
               main = paste("Module-trait relationships"))
dev.off()


library(clusterProfiler)
library(enrichplot)
library(tidyverse)
library(ggstatsplot)
library(ggnewscale)
library(org.Hs.eg.db)
library(DOSE)
diff=read.csv("venn221.diff.csv",row.names = 1,header = T,stringsAsFactors = FALSE,check.names = FALSE)
gene_up=unique(row.names(diff))
gene_diff_entrez <- as.character(na.omit(bitr(gene_up,
                                              fromType="SYMBOL", 
                                              toType="ENTREZID",
                                              OrgDb="org.Hs.eg.db")[,2])) 


#####kegg###
kegg_enrich_results <- enrichKEGG(gene  = gene_diff_entrez,
                                  organism  = "hsa", 
                                  pvalueCutoff = 0.05,
                                  qvalueCutoff = 0.2,pAdjustMethod = "none")

kk_read <- DOSE::setReadable(kegg_enrich_results, 
                             OrgDb="org.Hs.eg.db", 
                             keyType='ENTREZID')

a=kk_read@result
a=a[which(a$pvalue<0.05),]
write.csv(a,'KEGG_gene_enrichresults.csv')
save(kk_read, file ='KEGG_gene_enrichresults.Rdata')
go_enrich_results <- enrichGO(gene = gene_diff_entrez,
                              OrgDb = "org.Hs.eg.db",
                              ont   = "MF"  ,   
                              pvalueCutoff  = 0.05,
                              qvalueCutoff  = 0.2,pAdjustMethod = "none",
                              readable      = TRUE)
b=go_enrich_results@result

b=b[which(b$pvalue<0.05),]
write.csv(b, 'GO_gene_MF_enrichresults.csv') 
save(go_enrich_results, file ='GO_gene_MF_enrichresults.Rdata')

pt <- pairwise_termsim(go_enrich_results)
treep <- treeplot(pt,color = "pvalue",
                  showCategory = 30)
ggsave(treep, filename = 'treeplot_MF.pdf', width=15, height=8)
diff <- read.csv("venn.diff.csv",header = T,stringsAsFactors = FALSE,check.names = FALSE,row.names = 1)
genelist=as.numeric(diff$logFC)
names(genelist) <- row.names(diff)
pdf("kegg_cnetplot.pdf",width = 16,height = 12)
cnetplot(kk_read,  foldChange = genelist,
         showCategory = 10,
         node_label = 'all',
         circular = T, 
         colorEdge = T)
dev.off()


library(tidyverse)
library(glmnet)
source('msvmRFE.R') 
library(e1071)
library(caret)
library(randomForest)
exp=read.delim("../GSE135917_exp.txt",header = T,stringsAsFactors = FALSE,check.names = FALSE,row.names = 1)
venn=read.csv("../diff_snp.gene.venn.wgcna.csv",header = T,stringsAsFactors = FALSE,check.names = FALSE,row.names = 1)
s=which(row.names(exp)%in%row.names(venn))
exp1=exp[s,]
exp1=as.data.frame(t(exp1))
exp1$sample=row.names(exp1)
group=read.delim("../group.GSE135917.txt",header = T,stringsAsFactors = FALSE,check.names = FALSE)
exp2=merge(group,exp1,by="sample")
row.names(exp2)=exp2$sample
exp2$sample=NULL
write.csv(exp2,"lasso_svm_RF_input-1.csv")
train <- read.csv("lasso_svm_RF_input-1.csv",row.names = 1, 
                  as.is = F) 
dim(train)
x <- as.matrix(train[,-1])  
y <- ifelse(train$group == "Control", 0,1)
library(glmnet)
set.seed(123)
fit = glmnet(x, y, family = "binomial", alpha = 1, lambda = NULL)
pdf("A_lasso_model.pdf", width = 5, height = 5)
plot(fit, xvar = "lambda", label = TRUE)
dev.off()
cvfit = cv.glmnet(x, y, 
                  nfold=5, 
                  family = "binomial", type.measure = "class")


pdf("A_lasso_cvfit1.pdf", width = 6, height = 4)
plot(cvfit)
dev.off()
cvfit$lambda.min
myCoefs <- coef(cvfit, s="lambda.min");
lasso_fea <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
(lasso_fea <- lasso_fea[-1])
write.csv(lasso_fea,"feature_lasso_0.01428362.csv")
predict <- predict(cvfit, newx = x[1:nrow(x),], s = "lambda.min", type = "class")
table(predict,y)

input=train
input$group=ifelse(input$group == "Normal", 0,1)
set.seed(123456)
svmRFE(input, k = 5, halve.above = 100)
nfold = 5
nrows = nrow(input)
folds = rep(1:nfold, len=nrows)[sample(nrows)]
folds = lapply(1:nfold, function(x) which(folds == x))
results = lapply(folds, svmRFE.wrap, input, k=5, halve.above=100) 

top.features = WriteFeatures(results, input, save=F) 
head(top.features)

write.csv(top.features,"feature_svm.csv")

featsweep = lapply(1:16, FeatSweep.wrap, results, input) 
save(featsweep,file = "featsweep.RData")

no.info = min(prop.table(table(input[,1])))
errors = sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))
pdf("B_svm-error.pdf",width = 5,height = 5)
PlotErrors(errors, no.info=no.info) 
dev.off()
pdf("B_svm-accuracy.pdf",width = 5,height = 5)
Plotaccuracy(1-errors,no.info=no.info) 
dev.off()

library(randomForest)
df_train <- read.csv("lasso_svm_RF_input-1.csv",sep=",",header=T,check.names=F,row.names = 1)

set.seed(123)
df_train$group <- factor(df_train$group)

train.forest <- randomForest(group ~ ., 
                             data = df_train, 
                             ntree = 500,
                             importance = TRUE,
                             proximity = TRUE)
train.forest

pdf("RF_trees.pdf",width = 5,height = 5)
plot(train.forest,main="Random Forest")  
dev.off()
summary(train.forest)
importance_otu <- train.forest$importance
head(importance_otu)

importance_otu <- data.frame(importance(train.forest))

importance_otu <- importance_otu[order(importance_otu$MeanDecreaseAccuracy,decreasing = T),]
head(importance_otu)
write.csv(importance_otu,"RF_importance_otu.csv")


pdf("RF_top10.pdf",width = 6.5,height = 6.5)
varImpPlot(train.forest,n.var = min(10, nrow(train.forest$importance)), main = 'Top 10 - variable importance')
dev.off()
importance_otu=read.csv("RF_importance_otu.csv",header = T,stringsAsFactors = FALSE,check.names = FALSE,row.names = 1)
importance_otu$Feature=row.names(importance_otu)
importance_otu=importance_otu[1:10,]

colnames(importance_otu)
ggplot(importance_otu, aes(x= reorder( Feature,MeanDecreaseAccuracy), y=MeanDecreaseAccuracy,fill=Feature)) +
  geom_bar(stat="identity") +
  theme_classic() +
  guides(fill=FALSE)+
  scale_fill_manual(values=c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#D9D9D9","#BC80BD","#CCEBC5","#FFED6F","#377EB8","#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#D9D9D9","#BC80BD","#CCEBC5","#FFED6F","#377EB8"))+
  coord_flip()+
  theme_bw()+
  ggtitle('Top 10 - variable importance')+
  theme(plot.title = element_text(size=18,color='black'),
        axis.title.x =element_text(size=16,color='black'),
        axis.text.x =element_text(size=14, color='black'),
        axis.title.y =element_blank(),
        axis.text.y=element_text(size=14,   color='black'),
        legend.title=element_text(size=16, color='black' ),
        legend.text=element_text(size=14, color='black'),
        title=element_text(size=20, color='black'),
        strip.text = element_text(size = 14))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  labs(x="gene",y="MeanDecreaseAccuracy",fill="")
ggsave('RF_importance_TOP10.pdf',w=10,h=10)

library(VennDiagram)
lasso=read.csv("feature_lasso_0.01428362.csv",header = T,stringsAsFactors = FALSE,check.names = FALSE,row.names = 1)
lasso_fea=lasso$x
svm=read.csv("feature_svm.csv",header = T,stringsAsFactors = FALSE,check.names = FALSE,row.names = 1)
svm=svm[1:8,]
svm_fea=svm$FeatureName
rf=read.csv("RF_importance_otu.csv",header = T,stringsAsFactors = FALSE,check.names = FALSE,row.names = 1)
rf=rf[1:10,]
rf_fea=row.names(rf)
pdf("Venn_lasso_svm_RF.pdf", width = 6, height = 3)
venn.plot <- venn.diagram(list(LASSO = lasso_fea, #»?ͼ
                               SVM=svm_fea,
                               RF=rf_fea), NULL, 
                          fill = c("#E7B800",
                                   "#3E8E93",
                                   "#E31A1C"), 
                          cex = 2, cat.fontface=4, 
                          category.names = c("LASSO", 
                                             "SVM",
                                             "RandomForest"), 
                          main = "Overlap")
grid.draw(venn.plot)
dev.off()


library(foreign)
library(survival)
library(rms)
bc=read.csv("GSE135917_group_gene3.csv",header = T,stringsAsFactors = FALSE,row.names = 1)
bc$group=ifelse(bc$group == "Control", 0,1)
bc <- na.omit(bc)
dd <- datadist(bc)
options(datadist="dd")
colnames(bc)
formula1<-as.formula(group~ MSANTD1+PPP6C+PRCP)
fit1<-lrm(formula1,data = bc,x=T,y=T)
summary(fit1)
nom1<-nomogram(fit1,
               fun=function(x)1/(1+exp(-x)),
               lp=F,
               fun.at = c(0.1,0.5,0.9),
              funlabel = "Risk of OSA")
pdf("nomogram.pdf",height = 4,width = 6)
plot(nom1)
dev.off() 

cal1<-calibrate(fit1,method = "boot",B=1000)
pdf("nomogram_predict.pdf")
plot(cal1,xlim=c(0,1.0),ylim=c(0,1.0),
     xlab = "Nomogram Predicted Probability", ylab = "Actual Probability")
dev.off()

library(corrplot)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(psych)
uniq.symbol_sva=read.delim("GSE135917_exp.txt",header = T,stringsAsFactors = FALSE,check.names = FALSE,row.names = 1)
group=read.delim("group.GSE135917.txt",header = T,stringsAsFactors = FALSE,check.names = FALSE)
group1=group[which(group$group=="OSA"),]
data_FPKM=uniq.symbol_sva[,group1$sample]
rt <- data_FPKM

for (ii in c( "PPP6C",
"MSANTD1",
"PRCP")){  
  tar.exp <- rt[ii,]
  y <- as.numeric(tar.exp)
  data1 <- data.frame()
  for (i in rownames(rt)) {
    dd  <- corr.test(as.numeric(rt[i,]), y , method="pearson",adjust = "fdr")
    data1 = rbind(data1,data.frame(gene=i,cor=dd$r,p.value=dd$p))
  }
  data1 <- data1[order(data1$cor,decreasing = T),]
  write.csv(data1,file = paste0(ii,".pearson.cor.csv"))
}
  ii="PRCP"
  data1=read.csv(file = paste0(ii,".pearson.cor.csv"),header = T,stringsAsFactors = FALSE,check.names = FALSE)
  gene_cor=data1
  geneList <- gene_cor$cor

  names(geneList) = gene_cor$gene

  geneList = sort(geneList, decreasing = TRUE)
  head(geneList)
  tail(geneList)
  library(clusterProfiler)

  marks <- read.gmt("c5.go.bp.v7.4.symbols.gmt")

  y <- GSEA(geneList,TERM2GENE =marks)
  save(y,file=paste0(ii,"-GSEA-GOBP.rda"))
  write.table(y,file=paste0(ii,"-GSEA-GOBP.txt"),sep="\t",quote=F,row.names = F)   
  

  library(GseaVis)
  library(ggplot2)
  pdf(file = paste0(ii, "-GSEA-GOBP-TOP10.pdf"),width = 10, height = 8)
  gseaNb(y, geneSetID = y@result$Description[1:10],
         curveCol=c("#0000FF","#00FFFF","#04FF8A","#7CFC00","#AEFFB3",
                    "#F5DEB3","#F4A460","#FF6803","#FF6969","#FF0000"),
         pvalX = 0.9,pvalY = 0.9,subPlot = 3,rmHt = T)
  dev.off()
  

  marks <- read.gmt("c2.cp.kegg.v7.4.symbols.gmt")

  y1 <- GSEA(geneList,TERM2GENE =marks)
  save(y1,file=paste0(ii,"-GSEA-KEGG.rda"))
  write.table(y1,file=paste0(ii,"-GSEA-KEGG.txt"),sep="\t",quote=F,row.names = F)
  
  
  pdf(file = paste0(ii, "-GSEA-KEGG-TOP10.pdf"),width = 10, height = 8)
  gseaNb(y1, geneSetID = y1@result$Description[1:10],
         curveCol=c("#0000FF","#00FFFF","#04FF8A","#7CFC00","#AEFFB3",
                    "#F5DEB3","#F4A460","#FF6803","#FF6969","#FF0000"),
         pvalX = 0.9,pvalY = 0.9,subPlot = 3,rmHt = T)
  dev.off()


library(GSVA)
library(Biobase)
library(stringr)
uni_matrix=read.delim("GSE135917_exp.txt",header = T,stringsAsFactors = FALSE,check.names = FALSE,row.names = 1)
gene_immune=read.delim("cell.txt",header=T,stringsAsFactors = FALSE)
list<- split(as.matrix(gene_immune)[,1], gene_immune[,2])
gsva_matrix<- gsva(as.matrix(uni_matrix), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
immune_score=as.data.frame(t(gsva_matrix))
phenotype_file=read.delim("group.GSE135917.txt",header = T,stringsAsFactors = FALSE,check.names = FALSE)

immune_score$sample=row.names(immune_score)
imm=merge(phenotype_file,immune_score,by="sample")
row.names(imm)=imm$sample
imm$sample=NULL
write.csv(imm,"immune.cells.group.csv")
data=imm
test=as.list(data[,-1])
height<-stack(test)
Group=rep(data$group,ncol(data)-1)     
df=as.data.frame(cbind(Group,height))
colnames(df)=c("group","Infiltration_level","cells")
library(ggplot2)
library(ggpubr)
library(ggsignif)
pdf("ssgsea_boxplot.pdf",height=6,width=12)
ggplot(data=df,aes(x=cells,y=Infiltration_level,fill=group))+
  geom_boxplot(width=0.6)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_fill_manual(values = c( "#0073c2","#bc3c29"))+
  stat_compare_means(aes(group=group),
                     label="p.signif",
                     method = "wilcox.test",
                     hide.ns = T)+
  theme(axis.text.x = element_text(size=12,colour="black",angle=45,hjust = 1),axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=12),axis.title.y = element_text(size=12))

dev.off() 
score=uni_matrix
gene=read.csv("lasso_svm_rf/venn_gene3.csv",header = T,stringsAsFactors = FALSE,check.names = FALSE)
score=score[gene$x,]
score=as.data.frame(t(score))
cell_exp=read.csv("immune.cells.group.csv",header=T,row.names=1,stringsAsFactors=FALSE,check.names=FALSE)
cell_exp$group=NULL
cell_exp=cell_exp[row.names(score),]
all(row.names(score)==row.names(cell_exp))
data_m=score
data_mi=cell_exp
score <- c( )
immune_cells=c()
r <- c( )
p.value <- c( )
for(i in c(1:ncol(data_m))){
  for(j in c(1:ncol(data_mi))){
    x=as.numeric(data_m[,i])
    y =as.numeric(data_mi[ ,j])
    test = cor.test(x, y,method = "spearman")
    score = c(score, colnames(data_m)[i])
    immune_cells=c(immune_cells,colnames(data_mi)[j])
    r <- c(r, as.numeric(test$estimate))
    p.value <- c(p.value, test$p.value)
  }}
result <- data.frame(score,immune_cells, r, p.value)
write.csv(result, "gene_cells_correlation.csv",row.names = FALSE)
library(ggplot2)
table(result$score)
mtcars=result[which(result$score=="PCP2"),]
mtcars <- mtcars[order(mtcars$r), ]
mtcars$inflammation <- factor(mtcars$immune_cells, levels = mtcars$immune_cells)
pdf("GSE18897.spearman_cells.pdf",height=7,width=8)
ggplot(mtcars, aes(x=r,y=inflammation)) +
  geom_point(aes(color=r, size=abs(r)))  +
  scale_colour_gradient(low="#0073c2",high="#bc3c29")+
  scale_size_continuous(range = c(3, 8))+
  geom_segment(aes(x= 0, 
                   y =inflammation, 
                   yend = inflammation, 
                   xend = r), 
               color = "black") +
  xlab("Correlation Coefficient(r)")+ylab("")+
  theme_bw() +
  labs(title="PCP2",colour="Correlation")+
  theme(plot.title = element_text(hjust = 0.5,size = 16))+
  theme(axis.text.x = element_text(size=12,colour="black"),axis.text.y = element_text(size=12,colour="black"),
        axis.title.x = element_text(size=14),axis.title.y = element_text(size=14))+
  theme(#panel.grid.major =element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank())+
  geom_text(aes(label=ifelse(mtcars$p.value<0.001, "<0.001", round(mtcars$p.value,3)),x=max(r)+0.1))
dev.off() 

data=read.csv("gene_cells_correlation.csv",header=T,stringsAsFactors = FALSE)
data$pstar <- ifelse(data$p.value < 0.05,
                     ifelse(data$p.value < 0.01,"**","*"),
                     "")
library(ggplot2)
library(dplyr)
ggplot(data, aes( score,immune_cells)) + 
  geom_tile(aes(fill = r), colour = "white",size=1)+
  scale_fill_gradient2(low = "#2b8cbe",mid = "white",high = "#e41a1c")+
  geom_text(aes(label=paste(round(r,2),"\n",pstar)),col ="black",size = 3)+
  theme_minimal()+# 
  theme(axis.title.x=element_blank(),#
        axis.ticks.x=element_blank(),#
        axis.title.y=element_blank(),#
        axis.text.x = element_text(angle = 45, hjust = 1,size = 12,colour = "black"),# 
        axis.text.y = element_text(size = 12,colour = "black"))+#
  labs(fill =paste0(" * p < 0.05","\n\n","** p < 0.01","\n\n","Correlation"))
ggsave("gene_cells_correlation.pdf", width = 6, height = 8)

