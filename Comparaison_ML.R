#Mouhamadou Mansour Lo
#01/03/2022

## Comparaison des differentes methodes de reduction de dimensioon sur dataset artificiel ####


source(file="Simulation_ML.R")

donnees_swissroll <- SwissRoll(N=2000,Plot = FALSE)
donnees_brokenswissroll <- Broken_SwissRoll(N=2000,Plot = FALSE)
donnees_helix <- helix(N=2000,Plot = FALSE)
donnees_twinpeaks <- Twinpeaks(N=2000,Plot = FALSE)


#Les données vont etre comparés grace a l'erreur de generalisation du 1-NN entrainé dans l'espace de faible dimension obtenu a la suite des methodes


#Nous calculons les nouveaux vecteyrs puis nous les entrainons sur un 1-NN


#------------------------------------------------------------------------------------------------------
#Commencons d'abord par lle

tmp=sample(1:2000,500)

library(lle)
library(class)
library(caret)

donnees_swissroll.lle <-  lle(donnees_swissroll$samples, m = 2, k = 16) 


donnees_brokenswissroll.lle <-  lle(donnees_brokenswissroll$samples, m = 2, k = 19) 


donnees_helix.lle <-  lle(donnees_helix$samples, m = 1, k = 19) 


donnees_twinpeaks.lle <-  lle(donnees_twinpeaks$samples, m = 1, k = 14)


swissroll_train_X <- donnees_swissroll.lle$Y[-tmp,]
swissroll_train_Y <- donnees_swissroll$labels[-tmp]
swissroll_test_X <- donnees_swissroll.lle$Y[tmp,]
swissroll_test_Y <- as.factor(donnees_swissroll$labels[tmp])
swissrollknn <- knn(swissroll_train_X,swissroll_test_X,swissroll_train_Y,k=1)
confusionMatrix(swissrollknn ,swissroll_test_Y)
#Accuracy : 0.952  


donnees_brokenswissroll_train_X <- donnees_brokenswissroll.lle$Y[-tmp,]
donnees_brokenswissroll_train_Y <- donnees_brokenswissroll$labels[-tmp]
donnees_brokenswissroll_test_X <- donnees_brokenswissroll.lle$Y[tmp,]
donnees_brokenswissroll_test_Y <- as.factor(donnees_brokenswissroll$labels[tmp])
donnees_brokenswissrollknn <- knn(donnees_brokenswissroll_train_X,donnees_brokenswissroll_test_X,donnees_brokenswissroll_train_Y,k=1)
confusionMatrix(donnees_brokenswissrollknn ,donnees_brokenswissroll_test_Y)
#Accuracy : 0.938


donnees_helix_train_X <- as.matrix(donnees_helix.lle$Y[-tmp])
donnees_helix_train_Y <- donnees_helix$labels[-tmp]
donnees_helix_test_X <- as.matrix(donnees_helix.lle$Y[tmp])
donnees_helix_test_Y <- as.factor(donnees_helix$labels[tmp])
donnees_helixknn <- knn(donnees_helix_train_X,donnees_helix_test_X,donnees_helix_train_Y,k=1)
confusionMatrix(donnees_helixknn ,donnees_helix_test_Y)
#Accuracy : 0.196 


donnees_twinpeaks_train_X <- donnees_twinpeaks.lle$Y[-tmp,]
donnees_twinpeaks_train_Y <- donnees_twinpeaks$labels[-tmp]
donnees_twinpeaks_test_X <- donnees_twinpeaks.lle$Y[tmp,]
donnees_twinpeaks_test_Y <- as.factor(donnees_twinpeaks$labels[tmp])
donnees_twinpeaksknn <- knn(donnees_twinpeaks_train_X,donnees_twinpeaks_test_X,donnees_twinpeaks_train_Y,k=1)
confusionMatrix(donnees_twinpeaksknn ,donnees_twinpeaks_test_Y)
#Accuracy : 0.968  

#------------------------------------------------------------------------------------------------------

calcul_accuracy <- function(data,labels,num_dim){
  
tmp=sample(1:2000,500)

if(num_dim >= 2){
  
  train_X <- data[-tmp,]
  test_X <- data[tmp,]
  
}
else{
  train_X <- as.matrix(data[-tmp])
  test_X <- as.matrix(data[tmp])
}

train_Y <- labels[-tmp]
test_Y <- as.factor(labels[tmp])

knn <- knn(train_X,test_X,train_Y,k=1)

return(confusionMatrix(knn ,test_Y)$overall[1])
  
}

#------------------------------------------------------------------------------------------------------
#Ensuite MDS


donnees_swissroll.mds <- cmdscale(dist(donnees_swissroll$samples), k = 2, eig = FALSE, add = FALSE, 
                                  x.ret = T)

donnees_brokenswissroll.mds <- cmdscale(dist(donnees_brokenswissroll$samples), k = 2, eig = FALSE, add = FALSE, 
                                        x.ret = T)

donnees_helix.mds <- cmdscale(dist(donnees_helix$samples), k = 1, eig = FALSE, add = FALSE, 
                              x.ret = T)

donnees_twinpeaks.mds <- cmdscale(dist(donnees_twinpeaks$samples), k = 2, eig = FALSE, add = FALSE, 
                                  x.ret = T)

accuracy_swissroll.mds <- calcul_accuracy(donnees_swissroll.mds$points,donnees_swissroll$labels,2)
#0.576

accuracy_brokenswissroll.mds <- calcul_accuracy(donnees_brokenswissroll.mds$points,donnees_brokenswissroll$labels,2)
#0.704

accuracy_twinpeaks.mds <- calcul_accuracy(donnees_twinpeaks.mds$points,donnees_twinpeaks$labels,2)
#0.744

accuracy_helix.mds <- calcul_accuracy(donnees_helix.mds$points,donnees_helix$labels,1)
#0.424

#------------------------------------------------------------------------------------------------------

#Ensuite kernel pca


library(kernlab)
donnees_swissroll.kpca <- kpca(as.matrix(donnees_swissroll$samples), kernel="rbfdot", kpar = 
                                 list(sigma = 0.1), features = 2, th = 1e-4, na.action = na.omit)


donnees_brokenswissroll.kpca <- kpca(as.matrix(donnees_brokenswissroll$samples), kernel="rbfdot", kpar = 
                                       list(sigma = 0.1), features = 2, th = 1e-4, na.action = na.omit)


donnees_helix.kpca <- kpca(as.matrix(donnees_helix$samples), kernel="rbfdot", kpar = 
                             list(sigma = 0.1), features = 2, th = 1e-4, na.action = na.omit)


donnees_helix.kpca <- kpca(as.matrix(donnees_helix$samples), kernel="rbfdot", kpar = 
                             list(sigma = 0.1), features = 1, th = 1e-4, na.action = na.omit)


donnees_twinpeaks.kpca <- kpca(as.matrix(donnees_twinpeaks$samples), kernel="rbfdot", kpar = 
                                 list(sigma = 0.1), features = 2, th = 1e-4, na.action = na.omit)

accuracy_swissroll.kpca <- calcul_accuracy(donnees_swissroll.kpca@pcv,donnees_swissroll$labels,2)
#0.76

accuracy_brokenswissroll.kpca <- calcul_accuracy(donnees_brokenswissroll.kpca@pcv,donnees_brokenswissroll$labels,2)
#0.796

accuracy_helix.kpca <- calcul_accuracy(donnees_helix.kpca@pcv,donnees_helix$labels,1)
#0.266

accuracy_twinpeaks.kpca <- calcul_accuracy(donnees_twinpeaks.kpca@pcv,donnees_twinpeaks$labels,2)
#0.802

#------------------------------------------------------------------------------------------------------

#ISOMAP

library(vegan)

donnees_swissroll.isomap <- isomap(dist(donnees_swissroll$samples), ndim=2,k=7) #5 #7 #10

donnees_brokenswissroll.isomap <- isomap(dist(donnees_brokenswissroll$samples), ndim=2,k=30) #15 #13 #14

donnees_helix.isomap <- isomap(dist(donnees_helix$samples), ndim=1,k=5) #1,3 #1,5 #2,5 #2,5

donnees_twinpeaks.isomap <- isomap(dist(donnees_twinpeaks$samples), ndim=2,k=10) #5 #7 #10

accuracy_swissroll.isomap <- calcul_accuracy(donnees_swissroll.isomap$points,donnees_swissroll$labels,2)
#0.982

accuracy_brokenswissroll.isomap <- calcul_accuracy(donnees_brokenswissroll.isomap$points,donnees_brokenswissroll$labels,2)
#0.928

accuracy_helix.isomap <- calcul_accuracy(donnees_helix.isomap$points,donnees_helix$labels,1)
#0.232

accuracy_twinpeaks.isomap <- calcul_accuracy(donnees_twinpeaks.isomap$points,donnees_twinpeaks$labels,2)
#0.964

#------------------------------------------------------------------------------------------------------

#Et ENFIN LAPLACIAN EIGENMAPS


library(dimRed)


donnees_swissroll.LEM <- embed(donnees_swissroll$samples, "LaplacianEigenmaps",ndim=2,knn=10)#10 #20 #5


donnees_brokenswissroll.LEM <- embed(donnees_brokenswissroll$samples, "LaplacianEigenmaps",ndim=2,knn=25)#12 #20 #5


donnees_helix.LEM <- embed(donnees_helix$samples, "LaplacianEigenmaps",ndim=1,knn=10)#2,5 #1,5 #1,10 2,10


donnees_twinpeaks.LEM <- embed(donnees_twinpeaks$samples, "LaplacianEigenmaps",ndim=2,knn=7)#10 #5


accuracy_swissroll.LEM <- calcul_accuracy(donnees_swissroll.LEM@data@data,donnees_swissroll$labels,2)
#0.932

accuracy_brokenswissroll.LEM <- calcul_accuracy(donnees_brokenswissroll.LEM@data@data,donnees_brokenswissroll$labels,2)
#0.854

accuracy_helix.LEM  <- calcul_accuracy(donnees_helix.LEM@data@data,donnees_helix$labels,1)
#0.16

accuracy_twinpeaks.LEM  <- calcul_accuracy(donnees_twinpeaks.LEM@data@data,donnees_twinpeaks$labels,2)
#0.892
