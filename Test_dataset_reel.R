#Mouhamadou Mansour Lo
#25/02/2022

#Chargement du jeu de donn√©e MNIST

Train_image <- file("train-images.idx3-ubyte", "rb")
Train_Label <- file("train-labels.idx1-ubyte", "rb")

Test_image <- file("t10k-images.idx3-ubyte", "rb")
Test_Label <- file("t10k-labels.idx1-ubyte", "rb")

# helper function for visualization
show_digit = function(arr784, col = gray(12:1 / 12), ...) {
  image(matrix(as.matrix(arr784[-785]), nrow = 28)[, 28:1], col = col, ...)
}

#Code de chargement source
#https://gist.github.com/daviddalpiaz/ae62ae5ccd0bada4b9acd6dbc9008706
# load image files
load_image_file = function(f) {
  readBin(f, 'integer', n = 1, size = 4, endian = 'big')
  n    = readBin(f, 'integer', n = 1, size = 4, endian = 'big')
  nrow = readBin(f, 'integer', n = 1, size = 4, endian = 'big')
  ncol = readBin(f, 'integer', n = 1, size = 4, endian = 'big')
  x = readBin(f, 'integer', n = n * nrow * ncol, size = 1, signed = FALSE)
  close(f)
  data.frame(matrix(x, ncol = nrow * ncol, byrow = TRUE))
}

# load label files
load_label_file = function(f) {
  readBin(f, 'integer', n = 1, size = 4, endian = 'big')
  n = readBin(f, 'integer', n = 1, size = 4, endian = 'big')
  y = readBin(f, 'integer', n = n, size = 1, signed = FALSE)
  close(f)
  y
}

# load images
train = load_image_file(Train_image)
test  = load_image_file(Test_image)

# load labels
#train$y = as.factor(load_label_file(Train_Label))
#test$y  = as.factor(load_label_file(Test_Label))

# view test image
show_digit(train[1, ])



#-----------------------------------------------------------------------------

#Choix du nombre de dimension intresinque par Maximum Likelihood Estimation of Intrinsic Dimension


#Pour cela , on utilise le package intrinsicDimension

ind_sample <- sample(nrow(train), 2000)

sample_Mnist <- as.matrix(train[ind_sample,])

labels <- load_label_file(Train_Label)
label_sampl <- labels[ind_sample]

library(intrinsicDimension)


MLE_ID_MNIST <- maxLikGlobalDimEst(sample_Mnist, k = 5)

#-----------------------------------------------------------------------------

#Nous obtenons ~ 13 dimensions comme resultat


#Nous allons reduire la dimension avec PCA puis isomap

#Commencons d'abord par une PCA norm√©e

sample_Mnist.pca <- prcomp(sample_Mnist)
plot(sample_Mnist.pca$sdev,type = 'b')

library(factoextra)
fviz_eig(sample_Mnist.pca,ncp=25)

pc.use <- 13
trunc_MNIST <- sample_Mnist.pca$x[,1:pc.use] 

#-----------------------------------------------------------------------------

#Ensuite reduisons la dimension avec isomap

library(vegan)

sample_Mnist.isomap <- isomap(dist(sample_Mnist), ndim=13,k=7)

#-----------------------------------------------------------------------------

#Testons ces nouveaux jeux de donnÈes avec  1-NN , RandomForest et un package developpÈ dans le cadre dans projet de Model Based Learning

calcul_accuracy <- function(data,labels,num_dim){
  
  tmp=sample(1:2000,700)
  
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

#testons d'abord avec isomap
library(class)
library(caret)

accuracy_MNIST.isomap <- calcul_accuracy(sample_Mnist.isomap$points,label_sampl,13)
#0.881

accuracy_MNIST.pca <- calcul_accuracy(trunc_MNIST,label_sampl,13)
#0.834

#Testons avec Random forest

mnist_isomap = sample_Mnist.isomap$points

tmp=sample(1:2000,700)

train_X <- as.matrix(mnist_isomap[-tmp,])
test_X <- as.matrix(mnist_isomap[tmp,])

train_Y <- as.factor(label_sampl[-tmp])
test_Y <- as.factor(label_sampl[tmp])

fit = randomForest::randomForest(train_Y ~ ., data = train_X)
fit$confusion
test_pred = predict(fit, test_X)
mean(test_pred == test_Y)
table(predicted = test_pred, actual = test_Y)

#Accuracy 0.8885714

library(EMpackage)

mytest <- Classif_MM(train_X,train_Y,test_X,3)

t <- table(mytest$prediction,test_Y)
confusionMatrix(t)
#Accuracy : 0.8814 