#Mouhamadou Mansour Lo
#27/02/2022

## Estimation(Calibration) Manifold Learning ####

setwd("D:/OneDrive/Bureau")
source(file="Simulation_ML.R")

donnees_swissroll <- SwissRoll(N=1000,Plot = FALSE)
donnees_brokenswissroll <- Broken_SwissRoll(N=1000,Plot = FALSE)
donnees_helix <- helix(N=1000,Plot = FALSE)
donnees_twinpeaks <- Twinpeaks(N=1000,Plot = FALSE)

# 1 ere partie : Calibration par CorrDim

library(pdist)
# fonction qui d?termine la corr?lation dimension pour plusieurs Epsilon
corrDim <- function(data, epsilon = 10^seq(-2, 1, length.out = 100)){

  n <- nrow(data)
  lEps <- length(epsilon)

  C <- rep(0,lEps)
  distances <- as.matrix(pdist(data,data))
  for(k in 1:lEps){
    C[k] <- sum(distances<epsilon[k]) /n / (n-1)
  }


  return.df <- data.frame(C, epsilon)
  return.df$epsilon <- epsilon
  return.df$C <- C
  return(return.df)

}
# Fonction pour la d?riv?e
derivate <- function(x, y) {
  ll     <- length(y)
  deltax <- x[2] - x[1] # assumes equally spaced grid
  deltaf <- y[3:ll] - y[1:(ll - 2)]
  return(c(NA, deltaf / 2 / deltax, NA))
}

corrDim_swissroll <- corrDim(donnees_swissroll$samples)
corrDim_brokenswissroll <- corrDim(donnees_brokenswissroll$samples)
corrDim_helix <- corrDim(donnees_helix$samples)
corrDim_twinpeaks <- corrDim(donnees_twinpeaks$samples)



plot(log10(corrDim_swissroll$epsilon),
     derivate(log10(corrDim_swissroll$epsilon), log10(corrDim_swissroll$C)),
     type = 'l', xlim = c(-2,1))

abline(h = 1.93, untf = FALSE, col=2)

plot(log10(corrDim_brokenswissroll$epsilon),
     derivate(log10(corrDim_brokenswissroll$epsilon), log10(corrDim_brokenswissroll$C)),
     type = 'l', xlim = c(-2,1))

abline(h = 1.92, untf = FALSE, col=2)

plot(log10(corrDim_helix$epsilon),
     derivate(log10(corrDim_helix$epsilon), log10(corrDim_helix$C)),
     type = 'l', xlim = c(-2,1))

abline(h = 1, untf = FALSE, col=2)

plot(log10(corrDim_twinpeaks$epsilon),
     derivate(log10(corrDim_twinpeaks$epsilon), log10(corrDim_twinpeaks$C)),
     type = 'l', xlim = c(-2,1))

abline(h = 1.93, untf = FALSE, col=2)

#--------------------------------------------------------------------------------------------------------------

# 2 eme partie : Calibration par ACP Norm?e

donnees_swissroll.pca <- prcomp(donnees_swissroll$samples, center=T, scale. = T)
plot(donnees_swissroll.pca$sdev,type = 'b')

donnees_brokenswissroll.pca <- prcomp(donnees_brokenswissroll$samples, center=T, scale. = T)
plot(donnees_brokenswissroll.pca$sdev,type = 'b')

#Pour les deux premiers datasets nous ne remarquons pas vraiment de "coude"

donnees_helix.pca <- prcomp(donnees_helix$samples, center=T, scale. = T)
plot(donnees_helix.pca$sdev,type = 'b')

#Pour ce dataset , la variance entre les composantes principales reste constante

donnees_twinpeaks.pca <- prcomp(donnees_twinpeaks$samples, center=T, scale. = T)
plot(donnees_twinpeaks.pca$sdev,type = 'b')

#Ici nous avons un "coude" a 2

#--------------------------------------------------------------------------------------------------------------

# 3 eme partie : Calibration par kernelpca

#Pour cela , on utilise le package kernlab


library(kernlab)

donnees_swissroll.kpca <- kpca(as.matrix(donnees_swissroll$samples), kernel="rbfdot", kpar =
                    list(sigma = 0.1), features = 0, th = 1e-4, na.action = na.omit)
plot(eig(donnees_swissroll.kpca),  xlim = c(0,10),type='b')


donnees_brokenswissroll.kpca <- kpca(as.matrix(donnees_brokenswissroll$samples), kernel="rbfdot", kpar =
                                 list(sigma = 0.1), features = 0, th = 1e-4, na.action = na.omit)
plot(eig(donnees_brokenswissroll.kpca),  xlim = c(0,10),type='b')

donnees_helix.kpca <- kpca(as.matrix(donnees_helix$samples), kernel="rbfdot", kpar =
                                       list(sigma = 0.1), features = 0, th = 1e-4, na.action = na.omit)
plot(eig(donnees_helix.kpca),  xlim = c(0,10),type='b')

donnees_twinpeaks.kpca <- kpca(as.matrix(donnees_twinpeaks$samples), kernel="rbfdot", kpar =
                             list(sigma = 0.1), features = 0, th = 1e-4, na.action = na.omit)
plot(eig(donnees_twinpeaks.kpca),  xlim = c(0,10),type='b')


#--------------------------------------------------------------------------------------------------------------

# 4 eme partie : Calibration par Maximum Likelihood Estimation of Intrinsic Dimension

#Pour cela , on utilise le package intrinsicDimension

library(intrinsicDimension)


MLE_ID_swissroll <- maxLikGlobalDimEst(donnees_swissroll$samples, k = 10)

MLE_ID_brokenswissroll <- maxLikGlobalDimEst(donnees_brokenswissroll$samples, k = 10)

MLE_ID_helix <- maxLikGlobalDimEst(donnees_helix$samples, k = 5,unbiased = TRUE)

MLE_ID_twinpeaks <- maxLikGlobalDimEst(donnees_twinpeaks$samples, k = 10)

