#Mouhamadou Mansour Lo
#01/03/2022

## Estimation(Methodes) Manifold Learning ####

setwd("D:/OneDrive/Bureau")
source(file="Simulation_ML.R")

donnees_swissroll <- SwissRoll(N=2000,Plot = FALSE)
donnees_brokenswissroll <- Broken_SwissRoll(N=2000,Plot = FALSE)
donnees_helix <- helix(N=2000,Plot = FALSE)
donnees_twinpeaks <- Twinpeaks(N=2000,Plot = FALSE)


#---------------------------------------------------------------------------------------------------

# LLE


library(lle)

#calk_swissroll <- calc_k(donnees_swissroll$samples, 2) 
donnees_swissroll.lle <-  lle(donnees_swissroll$samples, m = 2, k = 16) 

require(rgl)
plot3d(donnees_swissroll.lle$X, xlab="x", ylab="y", zlab="z",col=donnees_swissroll$labels)
plot(donnees_swissroll.lle$Y, xlab="ti", ylab="z",col = donnees_swissroll$labels)


#calk_brokenswissroll <- calc_k(donnees_brokenswissroll$samples, 2) 
donnees_brokenswissroll.lle <-  lle(donnees_brokenswissroll$samples, m = 2, k = 19) 

require(rgl)
plot3d(donnees_brokenswissroll.lle$X, xlab="x", ylab="y", zlab="z",col=donnees_brokenswissroll$labels)
plot(donnees_brokenswissroll.lle$Y, xlab="ti", ylab="z",col = donnees_brokenswissroll$labels)


#calk_helix <- calc_k(donnees_helix$samples, 1) 
donnees_helix.lle <-  lle(donnees_helix$samples, m = 1, k = 19) 

require(rgl)
plot3d(donnees_helix.lle$X, xlab="x", ylab="y", zlab="z",col=donnees_helix$labels)
plot(donnees_helix.lle$Y, xlab="ti", ylab="z",col = donnees_helix$labels)


#calk_twinpeaks <- calc_k(donnees_twinpeaks$samples, 2) 
donnees_twinpeaks.lle <-  lle(donnees_twinpeaks$samples, m = 2, k = 14) 

require(rgl)
plot3d(donnees_twinpeaks.lle$X, xlab="x", ylab="y", zlab="z",col=donnees_twinpeaks$labels)
plot(donnees_twinpeaks.lle$Y, xlab="ti", ylab="z",col = donnees_twinpeaks$labels)

#---------------------------------------------------------------------------------------------------

# MDS

donnees_swissroll.mds <- cmdscale(dist(donnees_swissroll$samples), k = 2, eig = FALSE, add = FALSE, 
                     x.ret = T)
plot(donnees_swissroll.mds$points, xlab="ti", ylab="z",col = donnees_swissroll$labels)



donnees_brokenswissroll.mds <- cmdscale(dist(donnees_brokenswissroll$samples), k = 2, eig = FALSE, add = FALSE, 
                                  x.ret = T)
plot(donnees_brokenswissroll.mds$points, xlab="ti", ylab="z",col = donnees_brokenswissroll$labels)



donnees_helix.mds <- cmdscale(dist(donnees_helix$samples), k = 1, eig = FALSE, add = FALSE, 
                                        x.ret = T)
plot(donnees_helix.mds$points, xlab="ti", ylab="z",col = donnees_helix$labels)



donnees_twinpeaks.mds <- cmdscale(dist(donnees_twinpeaks$samples), k = 2, eig = FALSE, add = FALSE, 
                              x.ret = T)
plot(donnees_twinpeaks.mds$points, xlab="ti", ylab="z",col = donnees_twinpeaks$labels)


#---------------------------------------------------------------------------------------------------

# Kernel PCA


library(kernlab)
donnees_swissroll.kpca <- kpca(as.matrix(donnees_swissroll$samples), kernel="rbfdot", kpar = 
                    list(sigma = 0.1), features = 2, th = 1e-4, na.action = na.omit)

require(rgl)
plot3d(donnees_swissroll.kpca@xmatrix, xlab="x", ylab="y", zlab="z",col=donnees_swissroll$labels)
plot(rotated(donnees_swissroll.kpca), xlab="ti", ylab="z",col = donnees_swissroll$labels)


donnees_brokenswissroll.kpca <- kpca(as.matrix(donnees_brokenswissroll$samples), kernel="rbfdot", kpar = 
                                 list(sigma = 0.1), features = 2, th = 1e-4, na.action = na.omit)

require(rgl)
plot3d(donnees_brokenswissroll.kpca@xmatrix, xlab="x", ylab="y", zlab="z",col=donnees_brokenswissroll$labels)
plot(rotated(donnees_brokenswissroll.kpca), xlab="ti", ylab="z",col = donnees_brokenswissroll$labels)


donnees_helix.kpca <- kpca(as.matrix(donnees_helix$samples), kernel="rbfdot", kpar = 
                                       list(sigma = 0.1), features = 2, th = 1e-4, na.action = na.omit)

require(rgl)
plot3d(donnees_helix.kpca@xmatrix, xlab="x", ylab="y", zlab="z",col=donnees_helix$labels)
plot(rotated(donnees_helix.kpca), xlab="ti", ylab="z",col = donnees_helix$labels)


donnees_helix.kpca <- kpca(as.matrix(donnees_helix$samples), kernel="rbfdot", kpar = 
                             list(sigma = 0.1), features = 1, th = 1e-4, na.action = na.omit)

require(rgl)
plot3d(donnees_helix.kpca@xmatrix, xlab="x", ylab="y", zlab="z",col=donnees_helix$labels)
plot(rotated(donnees_helix.kpca), xlab="ti", ylab="z",col = donnees_helix$labels)


donnees_twinpeaks.kpca <- kpca(as.matrix(donnees_twinpeaks$samples), kernel="rbfdot", kpar = 
                             list(sigma = 0.1), features = 2, th = 1e-4, na.action = na.omit)

require(rgl)
plot3d(donnees_twinpeaks.kpca@xmatrix, xlab="x", ylab="y", zlab="z",col=donnees_twinpeaks$labels)
plot(rotated(donnees_twinpeaks.kpca), xlab="ti", ylab="z",col = donnees_twinpeaks$labels)


#---------------------------------------------------------------------------------------------------

# ISOMAP

library(vegan)

donnees_swissroll.isomap <- isomap(dist(donnees_swissroll$samples), ndim=2,k=7) #5 #7 #10
plot(donnees_swissroll.isomap$points, xlab="ti", ylab="z",col = rainbow(7,start = 0.2,end = 1)[donnees_swissroll$labels])


donnees_brokenswissroll.isomap <- isomap(dist(donnees_brokenswissroll$samples), ndim=2,k=30) #15 #13 #14
plot(donnees_brokenswissroll.isomap$points, xlab="ti", ylab="z",col = rainbow(7,start = 0.2,end = 1)[donnees_brokenswissroll$labels])

donnees_helix.isomap <- isomap(dist(donnees_helix$samples), ndim=1,k=5) #1,3 #1,5 #2,5 #2,5
plot(donnees_helix.isomap$points, xlab="ti", ylab="z",col = rainbow(7,start = 0.2,end = 1)[donnees_helix$labels])

donnees_twinpeaks.isomap <- isomap(dist(donnees_twinpeaks$samples), ndim=2,k=10) #5 #7 #10
plot(donnees_twinpeaks.isomap$points, xlab="ti", ylab="z",col = rainbow(7,start = 0.2,end = 1)[donnees_twinpeaks$labels])

#---------------------------------------------------------------------------------------------------

# LAPLACIAN EIGENMAPS

library(dimRed)


donnees_swissroll.LEM <- embed(donnees_swissroll$samples, "LaplacianEigenmaps",ndim=2,knn=10)#10 #20 #5
plot(donnees_swissroll.LEM@data@data, xlab="ti", ylab="z",col = donnees_swissroll$labels)

donnees_brokenswissroll.LEM <- embed(donnees_brokenswissroll$samples, "LaplacianEigenmaps",ndim=2,knn=25)#12 #20 #5
plot(donnees_brokenswissroll.LEM@data@data, xlab="ti", ylab="z",col = donnees_brokenswissroll$labels)

donnees_helix.LEM <- embed(donnees_helix$samples, "LaplacianEigenmaps",ndim=1,knn=10)#2,5 #1,5 #1,10 2,10
plot(donnees_helix.LEM@data@data, xlab="ti", ylab="z",col = donnees_helix$labels)


donnees_twinpeaks.LEM <- embed(donnees_twinpeaks$samples, "LaplacianEigenmaps",ndim=2,knn=7)#10 #5
plot(donnees_twinpeaks.LEM@data@data, xlab="ti", ylab="z",col = donnees_twinpeaks$labels)