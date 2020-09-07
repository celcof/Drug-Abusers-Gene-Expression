library(GEOquery)
library(factoextra)
library(randomForest)
library(MASS)
library(genefilter)
library(caret)
library(rScudo)
library(parallel)
library(e1071)
library(doParallel)
library(illuminaHumanv4.db) #installation from BiocManager in case

gse <- getGEO("GSE54839")
gse <- gse[[1]]
ex <- exprs(gse)

pca <- prcomp(t(ex))

# for plotting tableau is used:
# i will use the output of prcomp for performing the
# visualization on the platform

ind <- rep(c(rep("control",3),rep("abuser",3)),10)

# the split control-abuser is made
# by looking at the data set description

write.csv(data.frame(comp1=pca$x[,1], comp2=pca$x[,2], ind=ind), "pca.csv")

# scree plot. again, tableau is used so i will
# just output the importance dataframe

write.csv(t(summary(pca)$importance), "imp.csv")

# list of subject: to look whether different genetic maps of different
# people are located one to each other

sList <- c()
for (i in 1:20) {
  sList <- c(sList, rep(paste("S",i), 3))
}
write.csv(data.frame(comp1=pca$x[,1], comp2=pca$x[,2], sL = sList), "subject.csv")

# k-means

k <- 2
kmeans_result <- kmeans(t(ex), k)
kmeans_result <- kmeans_result$cluster
write.csv(data.frame(comp1=pca$x[,1], comp2=pca$x[,2], cluster=kmeans_result), "km.csv")

# hclust

dist_matrix <- dist(t(ex))
hc_result <- hclust(dist_matrix, method = "ave")
fviz_dend(hc_result, cex = 0.5, k=2, show_labels = F, k_colors=c("brown", "red"), main=NULL)

# rf
train <- c(rep(c(rep(T,3),rep(F,3),rep(F,3),rep(T,3),rep(T,3)),4))
rf <- randomForest(x=t(ex)[, train], y=as.factor(ind), ntree=100)
pr <- as.character(predict(rf, t(ex[,])))

checkPerf <- function(pr, ind) {
  perf <- c()
  for (i in 1:60){
    if (train[i]) {
      perf <- c(perf, "train")
    }
    else {
      if (pr[i] == ind[i]) {
        if (pr[i] == "abuser") {
          perf <- c(perf, "tp")
        } else {
          perf <- c(perf, "tn")
        }
      } else {
        if (pr[i] == "abuser") {
          perf <- c(perf, "fp")
        } else {
          perf <- c(perf, "fn")
        }
      }
    }
  }
  return(perf)
}

perf <- checkPerf(pr, ind)
write.csv(data.frame(comp1=pca$x[,1], comp2=pca$x[,2], perf=perf), "rf.csv")
write.csv(data.frame(rf$importance), "imp.csv")

# lda

tt <- rowttests(ex,as.factor(ind))
keepers <- which(tt$p.value<0.15)
ext <- ex[keepers,]
ext <- data.frame(t(ext))
mod <- lda(ind ~ ., data=ext, prior = c(0.5,0.5), subset = train)
pr <- as.character(predict(mod, ext[!train,])$class)
table(pred=as.factor(pr), as.factor(ind[!train]))
posPredValue(as.factor(pr), as.factor(ind[!train]), positive="abuser")
sensitivity(as.factor(pr), as.factor(ind[!train]), positive="abuser")

write.csv(data.frame(comp1=pca$x[,1], comp2=pca$x[,2], perf=perf), "lda.csv")

# lasso

dat <- t(ex)
ctrl <- trainControl(method = "cv",
                     classProbs = TRUE,
                     summaryFunction = twoClassSummary)
lassoFit <- train(dat[train,],
                  ind[train],
                  method = "glmnet",
                  family = "binomial",
                  tuneGrid = expand.grid(alpha = 1,
                                         lambda = seq(0,1,by=0.05)),
                  trControl = ctrl,
                  metric = "ROC")
plot(lassoFit)
pred <- predict(lassoFit, dat[!train,], type="raw",
                s-lassoFit$finalModel$lambdaOpt)
table(pred,as.factor(ind[!train]))
posPredValue(as.factor(pred), as.factor(ind[!train]), positive="abuser")
sensitivity(as.factor(pred), as.factor(ind[!train]), positive="abuser")
write.csv(data.frame(lassoFit$results), "lambda.csv")
pr <- predict(lassoFit, dat, type="raw",
              s-lassoFit$finalModel$lambdaOpt)
perf <- checkPerf(pr,ind)
write.csv(data.frame(comp1=pca$x[,1], comp2=pca$x[,2], perf=perf), "lasso.csv")

# rank-based signatures

cl <- parallel::makePSOCKcluster(4)
doParallel::registerDoParallel(cl)
try <- c(5, 10, 15, 20, 25)
model <- scudoModel(nTop = try, nBottom = try, N = 0.25)
control <- caret::trainControl(method = "cv", number = 6, summaryFunction = caret::multiClassSummary)
cvRes <- caret::train(x = t(ex)[train,], y = ind[train], method = model, trControl = control)
parallel::stopCluster(cl)

write.csv(cvRes$results, "rbs-cv.csv")

classRes <- scudoClassify(ex[,train], ex[,!train], N = 0.25,
                          nTop = 10, nBottom = 25,
                          trainGroups = as.factor(ind[train]), alpha = 0.5)
confusionMatrix(classRes$predicted, as.factor(ind[!train]))

pr <- as.character(classRes$predicted)
perf <- as.character(train)
for (i in 1:24){
  j <- which(train == F)[i]
  if (pr[i] == ind[j]) {
    if (pr[i] == "abuser") {
      perf[j] <- "tp" 
    } else {
      perf[j] <- "tn"
    }
  } else {
    if (pr[i] == "abuser") {
      perf[j] <- "fp"
    } else {
      perf[j] <- "fn"
    }
  }
}

write.csv(data.frame(comp1=pca$x[,1], comp2=pca$x[,2], perf=perf), "rbs.csv")

# save output of random forest importance text file

fileConn<-file("imp.txt")
imp <- row.names(rf$importance)
impU <- imp[order(rf$importance,decreasing=T)][1:200]
writeLines(impU, fileConn)
close(fileConn)

# pathfindr
library(KEGGREST)
library(KEGGgraph)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(pathfindR)
library(readr)

impWP <- tt[imp,]
impWP["initial_alias"] <- row.names(impWP)
gProfiler_hsapiens_7_9_2020_15_00_45 <- read_csv("gProfiler_hsapiens_7-9-2020_15-00-45.csv")
impWP <- merge(gProfiler_hsapiens_7_9_2020_15_00_45, impWP, all.x=TRUE, by="initial_alias")
impWP <- impWP[impWP$name != "None",c(3,8)]
aggregate(p.value~name, impWP, max)
pathfinder <- run_pathfindR(impWP,
              gene_sets = "BioCarta",
              visualize_enriched_terms = FALSE) 
cluster_enriched_terms(pathfinder)
term_gene_graph(pathfinder)
