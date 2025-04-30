libraries <- c("limma", "CorLevelPlot", "impute", "tidyr", "dplyr", "GEOquery", "DESeq2", "WGCNA", "tidyverse", "gridExtra", "preprocessCore", "httpuv", "devtools", "flashClust", "caret", "xgboost", "pROC")
lapply(libraries, library, character.only = TRUE)

# For dataset GSE42568
data<- read.csv("GSE42568.csv",row.names = 1)
# Data preprocessing
# Principal component analysis (PCA)
pca<- prcomp(t(data))
pca.dat<- pca$x
pca.var<- pca$sdev^2
pca.var.percent<-round(pca.var/sum(pca.var)*100,digits=2)
pca.dat<-as.data.frame(pca.dat)
ggplot(pca.dat,aes(PC1,PC2))+
  geom_point()+
  geom_text(label=rownames(pca.dat))+
  labs(x=paste0("PC1: ",pca.var.percent[1],"%"),y=paste0("PC2: ",pca.var.percent[2],"%"))
#Cluster dendrogram (CD)
gsg<-goodSamplesGenes(t(data))
summary(gsg)
gsg$allOK
table(gsg$goodGenes)
table(gsg$goodSamples)
data<- data[gsg$goodGenes==TRUE,]
tdata<- t(data)
dist(tdata)
htree<- hclust(dist(t(data)),method="average")
plot(htree)
# Differentially expressed genes (DEGs) identification
design <- cbind(Normal=1,cancervsNormal=rep(c(1,0),times=c(17,104)))
fit <- lmFit(data,design)
fit <- eBayes(fit)
top_gene<- topTable(fit,number=Inf, coef="cancervsNormal", adjust="BH")
de_genes <- topTable(fit,number=Inf, coef = "cancervsNormal", adjust.method = "BH", sort.by = "none")
de_genes_1 <- de_genes[abs(de_genes$logFC) >= 1 & de_genes$adj.P.Val <= 0.05, ]

# For dataset GSE61304
data<- read.csv("GSE61304.csv",row.names = 1)
# Data preprocessing
# Principal component analysis (PCA)
pca<- prcomp(t(data))
pca.dat<- pca$x
pca.var<- pca$sdev^2
pca.var.percent<- round(pca.var/sum(pca.var)*100,digits=2)
pca.dat<- as.data.frame(pca.dat)
ggplot(pca.dat,aes(PC1,PC2))+
  geom_point()+
  geom_text(label=rownames(pca.dat))+
  labs(x=paste0("PC1: ",pca.var.percent[1],"%"),y=paste0("PC2: ",pca.var.percent[2],"%"))
#Cluster dendrogram (CD)
gsg<- goodSamplesGenes(t(data))
summary(gsg)
gsg$allOK
table(gsg$goodGenes)
table(gsg$goodSamples)
data<- data[gsg$goodGenes==TRUE,]
tdata<- t(data)
dist(tdata)
htree<- hclust(dist(t(data)),method="average")
plot(htree)
# Differentially expressed genes (DEGs) identification
design <- cbind(Normal=1,cancervsNormal=rep(c(1,0),times=c(4,58)))
fit <- lmFit(data,design)
fit <- eBayes(fit)
top_gene<- topTable(fit,number=Inf, coef="cancervsNormal", adjust="BH")
de_genes <- topTable(fit,number=Inf, coef ="cancervsNormal", adjust.method = "BH", sort.by = "none")
de_genes_2 <- de_genes[abs(de_genes$logFC) >= 1 & de_genes$adj.P.Val <= 0.05, ]

# For dataset GSE29431
data<- read.csv("GSE29431.csv",row.names = 1)
# Data preprocessing
# Principal component analysis (PCA)
pca<- prcomp(t(data))
pca.dat<- pca$x
pca.var<- pca$sdev^2
pca.var.percent<- round(pca.var/sum(pca.var)*100,digits=2)
pca.dat<- as.data.frame(pca.dat)
ggplot(pca.dat,aes(PC1,PC2))+
  geom_point()+
  geom_text(label=rownames(pca.dat))+
  labs(x=paste0("PC1: ",pca.var.percent[1],"%"),y=paste0("PC2: ",pca.var.percent[2],"%"))
#Cluster dendrogram (CD)
gsg<- goodSamplesGenes(t(data))
summary(gsg)
gsg$allOK
table(gsg$goodGenes)
table(gsg$goodSamples)
data<-data[gsg$goodGenes==TRUE,]
tdata<- t(data)
dist(tdata)
htree<- hclust(dist(t(data)),method="average")
plot(htree)
# Differentially expressed genes (DEGs) identification
design <- cbind(Normal=1,cancervsNormal=rep(c(1,0),times=c(12,54)))
fit <- lmFit(data,design)
fit <- eBayes(fit)
top_gene<-  topTable(fit,number=Inf, coef="cancervsNormal", adjust="BH")
de_genes <- topTable(fit,number=Inf, coef = "cancervsNormal", adjust.method = "BH", sort.by = "none")
de_genes_3 <- de_genes[abs(de_genes$logFC) >= 1 & de_genes$adj.P.Val <= 0.05, ]

#Common Differentially expressed genes (cDEGs) identification
cDEGs_1=rownames(de_genes_1)
cDEGs_2<-rownames(de_genes_2)
cDEGs_3<-rownames(de_genes_3)
commonDEGs<- intersect(cDEGs_1,cDEGs_2)
cDEGs<- intersect(commonDEGs,cDEGs_3)
# For weighted gene co-expression network analysis (WGCNA)
phenoData<- read.csv("missing_free_phenoData_GSE42568.csv",row.names = 1)
data<- read.csv("missing_free_data_GSE42568.csv",row.names = 1)
data_frame<- data.frame(data)
data_mod <- data_frame[rownames(data_frame) %in% cDEGs, ] 
datExpr<- data_mod
tdatExpr<- t(data_mod)
datExpr <- as.data.frame(t(datExpr))
dim(datExpr)
datTraits<- phenoData
powers <- c(c(1:10), seq(from =10, to=30, by=2))
sft <- pickSoftThreshold(datExpr, powerVector=powers, verbose =5, networkType="signed")
sizeGrWindow(9,5)
par(mfrow= c(1,2))
cex1=0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab= substitute(paste(bold("Soft Threshold (power)"))), ylab=substitute(paste(bold("Scale Free Topology Model Fit, signed R^2"))),font=2, type= "n", main= paste("Scale Independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, cex=cex1,font=2, col="red")
abline(h=0.80, col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab= substitute(paste(bold("Soft Threshold (power)"))), ylab=substitute(paste(bold("Mean Connectivity"))),font=2, type="n", main = paste("Mean Connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,font=2, col="red")
enableWGCNAThreads()
softPower <- 16
adjacency <- adjacency(datExpr, power = softPower, type = "signed") 
head(adjacency)
TOM <- TOMsimilarity(adjacency, TOMType="signed")
dissTOM <- 1-TOM
geneTree <- flashClust(as.dist(dissTOM), method="average")
plot(geneTree, xlab="",font=2,ylab = substitute(paste(bold('Height'))), sub="", main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE, hang=0.04)

minModuleSize <- 20
dynamicMods = cutreeDynamic(dendro= geneTree, distM= dissTOM, deepSplit=2, pamRespectsDendro= FALSE, minClusterSize = minModuleSize)
dynamicColors= labels2colors(dynamicMods)
table(dynamicColors)
MEList= moduleEigengenes(datExpr, colors= dynamicColors,softPower = softPower)
MEs= MEList$eigengenes
MEDiss= 1-cor(MEs)
METree= flashClust(as.dist(MEDiss), method= "average")
save(dynamicMods, MEList, MEs, MEDiss, METree, file= "Network_allSamples_signed_RLDfiltered.RData")
plot(METree, main= "Clustering of module eigengenes",font=2, xlab= "", sub= "")
MEDissThres = 0.25
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight= MEDissThres, verbose =3)
mergedColors = merge$colors
mergedMEs = merge$newMEs

plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged Dynamic"), dendroLabels= FALSE,font=2,ylab = substitute(paste(bold('Height'))),cex.colorLabels=0.8, hang=0.03, addGuide= TRUE, guideHang=0.05)
moduleColors = mergedColors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use= "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
textMatrix= paste(signif(moduleTraitCor, 2), "\n(",
                  signif(moduleTraitPvalue, 1), ")", sep= "")
dim(textMatrix)= dim(moduleTraitCor)
par(mar= c(4.5, 9,1.5, 1))
labeledHeatmap(Matrix= moduleTraitCor,
               xLabels= names(datTraits),
               cex.lab = 1,
               font.lab.x = 2,
               font.lab.y = 2,
               font.text=2,
               font.axis=2,
               yLabels= names(MEs),
               ySymbols= names(MEs),
               colorLabels= FALSE,
               colors= blueWhiteRed(50),
               textMatrix= textMatrix,
               setStdMargins= FALSE,
               cex.text= 0.8,
               zlim= c(-0.5,.5),
               main= paste("Module-Trait Relationships"))
# For predictive accuracy analysis of hub genes (HubGs)
data <- read.csv('pred_accuracy.csv')
features <- c('TOP2A', 'CCNB1', 'CCNB2', 'CCNA2', 'BUB1B', 'CDCA8', 'CDK1', 'PBK', 'TTK', 'MELK')
X <- data[, features]
y <- as.factor(data$Label)  
levels(y) <- c("control", "case")  
preProc <- preProcess(X, method = c("center", "scale"))
X_scaled <- predict(preProc, X)
set.seed(42)
trainIndex <- createDataPartition(y, p = 0.7, list = FALSE)
X_train <- X_scaled[trainIndex, ]
X_test <- X_scaled[-trainIndex, ]
y_train <- y[trainIndex]
y_test <- y[-trainIndex]
calculate_metrics <- function(y_true, y_pred, y_prob = NULL) {
  cm <- confusionMatrix(y_pred, y_true, positive = "case")
  metrics <- list(
    Accuracy = cm$overall['Accuracy'],
    Precision = cm$byClass['Pos Pred Value'],
    Recall_Sensitivity = cm$byClass['Sensitivity'],
    Specificity = cm$byClass['Specificity'],
    F1_Score = cm$byClass['F1'],
    AUC = NA  
  )
  if (!is.null(y_prob) && all(!is.na(y_prob))) {
    tryCatch({
      roc_obj <- roc(response = as.numeric(y_true) - 1, predictor = y_prob)
      metrics$AUC <- roc_obj$auc
    }, error = function(e) {
      warning(paste("ROC calculation failed:", e$message))
    })
  }
  return(metrics)
}
models <- list(
  SVM = list(
    method = "svmRadial",
    tuneGrid = data.frame(sigma = 0.1, C = 1),
    trControl = trainControl(method = "none", classProbs = TRUE),
    preProcess = c("center", "scale")
  ),
  RandomForest = list(
    method = "rf",
    tuneGrid = data.frame(mtry = floor(sqrt(length(features)))),
    trControl = trainControl(method = "none", classProbs = TRUE),
    ntree = 300,
    preProcess = c("center", "scale")
  ),
  XGBoost = list(
    method = "xgbTree",
    tuneGrid = data.frame(
      nrounds = 300,
      max_depth = 6,
      eta = 0.3,
      gamma = 0,
      colsample_bytree = 1,
      min_child_weight = 1,
      subsample = 1
    ),
    trControl = trainControl(method = "none", classProbs = TRUE)
  )
)
results <- list()
for (name in names(models)) {
  model <- train(
    x = X_train,
    y = y_train,
    method = models[[name]]$method,
    tuneGrid = models[[name]]$tuneGrid,
    trControl = models[[name]]$trControl,
    preProcess = models[[name]]$preProcess,
    ntree = if (name == "RandomForest") models[[name]]$ntree else NULL
  )
  y_pred <- predict(model, newdata = X_test)
  y_prob <- tryCatch({
    prob <- predict(model, newdata = X_test, type = "prob")[, "case"]
    if (any(is.na(prob))) NULL else prob
  }, error = function(e) {
    warning(paste("Probability prediction failed for", name, ":", e$message))
    NULL
  })
  results[[name]] <- calculate_metrics(y_test, y_pred, y_prob)
}
print_detailed_metrics <- function(metrics, model_name) {
  cat("\n", model_name, "Performance Metrics:\n")
  cat("============================================================\n")
  cat("Accuracy:", sprintf("%.4f", metrics$Accuracy), "\n")
  cat("Precision:", sprintf("%.4f", metrics$Precision), "\n")
  cat("Recall/Sensitivity:", sprintf("%.4f", metrics$Recall_Sensitivity), "\n")
  cat("F1-Score:", sprintf("%.4f", metrics$F1_Score), "\n")
  if (!is.na(metrics$AUC)) {
    cat("AUC:", sprintf("%.4f", metrics$AUC), "\n")
  } else {
    cat("AUC: Not available\n")
  }
  
  cat("\n")
  cat("============================================================\n")
}
for (name in names(results)) {
  print_detailed_metrics(results[[name]], name)
}
