
#sva去批次
library(sva)
library(FactoMineR)
library(factoextra)
library(limma)
rt1=read.table("GSE63060_exp.txt",row.names = 1,header = T)
GSE63060=t(rt1[,-5])
rt2=read.table("GSE140829_exp.txt",row.names = 1,header = T)
GSE140829=t(rt2[,-5])
all.data=cbind(GSE63060,GSE140829)
group=read.table("group.txt",header = T,check.names = F)
library(FactoMineR)
pre.pca <- PCA(t(all.data),graph = FALSE)
fviz_pca_ind(pre.pca,
             geom= "point",
             col.ind = group$batch,
             addEllipses = TRUE,
             legend.title="Group"  )
design = model.matrix(~group,data = group)
combat_Expr <- ComBat(dat = all.data,batch = group$batch, mod = design)
pre.pca <- PCA(t(combat_Expr),graph = FALSE)
fviz_pca_ind(pre.pca,
             geom= "point",
             col.ind = group$batch,
             addEllipses = TRUE,
             legend.title="Group"  )
GSE63060=combat_Expr[,1:ncol(GSE63060)]
GSE140829=combat_Expr[,(ncol(GSE63060)+1):ncol(combat_Expr)]
GSE63060=as.data.frame(t(GSE63060))
GSE140829=as.data.frame(t(GSE140829))
rownames(group)=group[,1]
same1=intersect(rownames(GSE63060),rownames(group))
group1=group[same1,]
GSE63060=GSE63060[same1,]
GSE63060$group=group1$group
same2=intersect(rownames(GSE140829),rownames(group))
group2=group[same2,]
GSE140829=GSE140829[same2,]
GSE140829$group=group2$group
################################################################
#构建模型
## 导入包
library(neuralnet)
library(NeuralNetTools)
library(ggplot2)
library(plotROC)
library(pROC)
library(ggpol)
library(RSNNS)
library(caret)
traindata = GSE63060
traindata$group = ifelse(traindata$group=="Control",0,1)
traindata_target = decodeClassLabels(traindata$group)

testdata = GSE140829
testdata$group = ifelse(testdata$group=="Control",0,1)
testdata_target = decodeClassLabels(testdata$group)

train_x = traindata[,1:4]
test_x = testdata[,1:4]

train_x = normalizeData(train_x,type = "0_1")
test_x = normalizeData(test_x,type = "0_1")
set.seed(77)  
model = mlp(train_x, traindata_target, size=5, learnFuncParams=c(0.001), 
             maxit=60,inputsTest=test_x, targetsTest=testdata_target)


pdf("net_plot.pdf",5,5,family = "serif")
par(cex = 0.6)
plotnet(model,pos_col = "red", neg_col = "grey")
dev.off()

traindata$probe=predict(model,train_x)[,1]
trainROC<-roc(response=traindata$group,predictor=traindata$probe)
trainAUC=round(auc(trainROC),3)
trainCI=ci.auc(trainROC)


testdata$probe=predict(model,test_x)[,1]
testROC<-roc(response=testdata$group,predictor=testdata$probe)
testAUC=round(auc(testROC),3)
testCI=ci.auc(testROC)

pdf("TrainROC.pdf",5,5)
plot(trainROC,
  print.auc=TRUE,
     grid=F,
     legacy.axes=T,
     auc.polygon=F,
     max.auc.ploygon=F,
     main="Train group",
     col="red")
text(x=0.35,y=0.42,labels=paste0("95%CI: ",sprintf("%0.3f", trainCI[1]),
                           "-",sprintf("%0.3f", trainCI[3])),
     cex=1,
     col="red")
dev.off()

pdf("TestROC.pdf",5,5,family = "Arial")
plot(testROC,
     print.auc=TRUE,
     grid=F,
     legacy.axes=T,
     auc.polygon=F,
     max.auc.ploygon=F,
     main="Test group",
     col="red")
text(x=0.35,y=0.42,labels=paste0("95%CI: ",sprintf("%0.3f", testCI[1]),
                                 "-",sprintf("%0.3f", testCI[3])),
     cex=1,
     col="red")
dev.off()

