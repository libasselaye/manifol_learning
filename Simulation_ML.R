#Mouhamadou Mansour Lo
#26/02/2022

## Simulation Manifold Learning ####

make_label <- function(X){
  
  maxi <- max(X)
  val <- floor(maxi/6)
  labels1 <- c()
  labels2 <- c()
  
  for(i in 1:length(X)){
    
    labels1 <- c(labels1,floor(X[i]/val))
  }
  
  unique <- sort(unique(labels1),decreasing = FALSE)
  
  tabassign <- cbind(unique,c(1:length(unique)))
  
  for(i in 1:length(X)){
    
    row <- which(grepl(labels1[i], tabassign[,1]))
    
    labels2 <- c(labels2,tabassign[row,2])

  }
    
  return(labels2)
}

SwissRoll <- function(N=5000,Plot=TRUE){
  
  q_i = runif(N, 0 , 1)
  p_i = runif(N, 0 , 1)
  ti = (3 * pi / 2) * (1 + 2 * p_i)
  
  x = ti * cos(ti)
  y = ti * sin(ti)
  z = 30 * q_i
  
  samples = cbind(x, y, z)
  
  labels <- make_label(ti)
  
  ## plot and return samples
  if(Plot){
    require(rgl)
    plot3d(samples, xlab="x", ylab="y", zlab="z",col = rainbow(7,start = 0.2,end = 1)[labels]);
  }
  
  return(list(samples=samples,labels=labels,t=cbind(ti,z)))
  
}

Broken_SwissRoll <- function(N=5000,Plot=TRUE){
  
  p_i_1 = runif(ceiling(N/2), 0 , 1)
  p_i_2 = runif(floor(N/2), 0 , 1)
  q_i_1 = runif(ceiling(N/2), 0 , 1)
  q_i_2 = runif(floor(N/2), 0 , 1)
  
  ti_1 = (3 * pi / 2) * (1 + 2 * p_i_1 * 0.4)
  ti_2 = (3 * pi / 2) * (1 + 2 * (p_i_2 * 0.4 + 0.6))
  
  ti = c(ti_1,ti_2)
  
    
  x = c(ti_1 * cos(ti_1),ti_2 * cos(ti_2))
  y = c(ti_1 * sin(ti_1),ti_2 * sin(ti_2))
  z = c(30 * q_i_1,30 * q_i_2)
  
  samples = cbind(x, y, z)
  
  labels <- make_label(ti)
  
  ## plot and return samples
  if(Plot){
    require(rgl)
    plot3d(samples, xlab="x", ylab="y", zlab="z",col = rainbow(7,start = 0.2,end = 1)[labels]);
  }
  
  return(list(samples=samples,labels=labels,t=cbind(ti,z)))
  
}


helix <- function(N=5000,Plot=TRUE){
  
  p_i = runif(N, 0 , 1)
  ti <- c(1:N)/N
  ti <- ti*2*pi
  
  x = (2+cos(8*ti))*cos(ti)
  y = (2+cos(8*ti))*sin(ti)
  z = sin(8*ti)
  
  samples = cbind(x, y, z)
  
  labels <- make_label(ti) 
  
  ## plot and return samples
  if(Plot){
    require(rgl)
    plot3d(samples, xlab="x", ylab="y", zlab="z",col = rainbow(7,start = 0.2,end = 1)[labels]);
  }
  
  return(list(samples=samples,labels=labels,t=cbind(ti,z)))
  
}


Twinpeaks <- function(N=5000,Plot=TRUE){
  
  q_i = runif(N, 0 , 1)
  p_i = runif(N, 0 , 1)
  
  
  x = 1 - 2*p_i
  y = 1 - 2*q_i
  
  z = sin(pi*x)* tanh(3*y)
  
  samples = cbind(x,y,z)
  
  labels <- make_label(x*30 +50)
  
  ## plot and return samples
  if(Plot){
    require(rgl)
    plot3d(samples, xlab="x", ylab="y", zlab="z",col=rainbow(7,start = 0.2,end = 1)[labels]);
  }

  return(list(samples=samples,labels=labels,t=cbind(x,y)))
  
}

donnees_swissroll <- SwissRoll(N=5000)
plot(donnees_swissroll$t, xlab="ti", ylab="z",col = rainbow(7,start = 0.2,end = 1)[donnees_swissroll$labels])
donnees_brokenswissroll <- Broken_SwissRoll(N=5000)
plot(donnees_brokenswissroll$t, xlab="ti", ylab="z",col = rainbow(7,start = 0.2,end = 1)[donnees_brokenswissroll$labels])
donnees_helix <- helix(N=5000)
plot(donnees_helix$t, xlab="ti", ylab="z",col = rainbow(7,start = 0.2,end = 1)[donnees_helix$labels])
donnees_twinpeaks <- Twinpeaks(N=5000)
plot(donnees_twinpeaks$t, xlab="ti", ylab="z",col = rainbow(7,start = 0.2,end = 1)[donnees_twinpeaks$labels])

