###### ------------------------------------------------------------
###### ------------------------------------------------------------
###### NAME: Spectral Degree Reduction based on Juan's publication
###### DATE: July 2018
###### Version: 1
###### ------------------------------------------------------------
###### ------------------------------------------------------------
###### COMMENTS:
###### The function SpecEmb from the LOE package works, both in finding the nearest partners
###### and doing the reduction on the small example
###### I cannot get it to work on the spehere example, it needs some more work to figure out
###### what is wrong with the function
###### I will continue working on this
###### but also, there is another point: what does the orthogonality of the base 
###### with which the lower dimention manifold is created have to do with the independence 
###### or conditional independence of the variables or dimensions?
###### There might be more information on the original article:
###### http://web.cse.ohio-state.edu/~belkin.8/papers/LEM_NC_03.pdf
###### Here are the slides from the presentation:
###### https://juanitorduz.github.io/documents/orduz_pydata2018.pdf
###### and here is the OP
###### https://juanitorduz.github.io/laplacian_eigenmaps_dim_red.html
if(!require(dimRed)) install.packages("dimRed",repos = "http://cran.us.r-project.org")
library("dimRed")
if(!require(diffusionMap)) install.packages("diffusionMap",repos = "http://cran.us.r-project.org")
library("diffusionMap")

if(!require(Graph)) install.packages("igraph",repos = "http://cran.us.r-project.org")
library("igraph")

if(!require(loe)) install.packages("loe",repos = "http://cran.us.r-project.org")
library("loe")

if(!require(ggplot2)) install.packages("ggplot2",repos = "http://cran.us.r-project.org")
library("ggplot2")

if(!require(scatterplot3d)) install.packages("scatterplot3d",repos = "http://cran.us.r-project.org")
library("scatterplot3d")


if(!require(FNN)) install.packages("FNN",repos = "http://cran.us.r-project.org")
library("FNN")

if(!require(reshape)) install.packages("reshape",repos = "http://cran.us.r-project.org")
library("reshape")


adj_mat <- matrix(c(0,1,1,1,1,0,0,0,1,0,0,0,1,0,0,0), nrow = 4, ncol = 4)
deg_mat <- Diagonal(rowsum(adj_mat,group = c(1,1,1,1)), n = 4)
laplacian <- as.matrix(deg_mat-adj_mat)

#GARI(adj_mat, deg_mat)
#
y <- spec.emb(adj_mat,1)

lambda_1 <- 1

a1 <- laplacian %*% y
a2 <- lambda_1*(deg_mat %*% y)
result <- a1 - a2



#make sphere in R

number_points <- 10
r <- runif(number_points, min = 0, max = 2*pi)
rho <- runif(number_points, min = 0, max = 2*pi)

x <- sin(rho)*cos(r)
y <- sin(rho)*sin(r)
z <- cos(rho)
data <- as.data.frame(matrix(c(x, y, z), ncol =3))
colnames(data) <- c("x_ax","y_ax", "z_ax")


x <- data$x_ax[!(abs(data$z_ax)>0.95)]
y <- data$y_ax[!(abs(data$z_ax)>0.95)]
z <- data$z_ax[!(abs(data$z_ax)>0.95)]

data_nopo <- as.data.frame(matrix(c(x, y, z), ncol =3))
colnames(data_nopo) <- c("x_ax","y_ax", "z_ax")
x1 <- rnorm(length(data_nopo$x_ax))

s_colors <- rgb(abs(data_nopo$x_ax), abs(data_nopo$y_ax), abs(data_nopo$z_ax), maxColorValue = 1)
scatterplot3d(data_nopo, pch = 16, main="Sphere",
              xlab = "X",
              ylab = "Y",
              zlab = "Z",
              box = FALSE,
              color = s_colors)

ggplot2::ggplot(data = data_nopo, aes(x = x_ax, y = y_ax)) + geom_point( color = s_colors)

n_components <- 2
prin <- prcomp(x = data_nopo, scale. = T, rank. = n_components)
#prin$center
#biplot(prin)
prinx <- prin$x[,1]
priny <- prin$x[,2]
prin_comp <- as.data.frame(matrix(c(prinx, priny), ncol = 2))
colnames(prin_comp) <- c("first", "second")

ggplot2::ggplot(data = prin_comp, aes(x = first, y = second)) + geom_point( color = s_colors)

n_neighbors <- 2

# Fit the object.
nn <- get.knn(data = data_nopo, n_neighbors) 
nn <- as.data.frame(nn)
nn$id <- rownames(nn)
nn_melt <- melt(nn, id.vars = c("id"))
nn_melt <- melt(nn)
nn_melt <- nn_melt[c("id", "value")]
adj_mat <- as.data.frame(as.matrix(get.adjacency(graph.edgelist(as.matrix(nn_melt), directed=TRUE))))
colnames(adj_mat)
head(as.matrix(adj_mat))
View((adj_mat))
t <- sort(as.integer(colnames(adj_mat)))
adj_mat <- adj_mat[order(colnames(adj_mat))]
adj_mat <- sort()
#adj_mat2 <- make.distmat(data_nopo)
#
#

devtools::install_github(c("duncantl/XMLRPC", "duncantl/RWordPress"))
library(RWordPress)
library(knitr)


# If need be, set your working directory to the location where you stored the Rmd file. 
setwd("/Users/sebastianmartinez/Dropbox/0. UoG/Projects/spectraldimreduc")




spec_emb <- spec.emb(A = as.matrix(adj_mat), p = n_components)
spec1 <- spec_emb[,1]
spec2 <- spec_emb[,2]
spec <- as.data.frame(matrix(c(spec1, spec2), ncol = 2))
colnames(spec) <- c("first", "second")


ggplot2::ggplot(data = spec, aes(x = first, y = second)) + geom_point( color = s_colors)
plot(se[,1], se[,2])





lpvs <- matrix(rnorm(25), 5, 5)
lpvs <- apply(lpvs, 2, function(x) { return (abs(x)/sqrt(sum(x^2))) })
lpvs <- apply(adj_mat, 2, function(x) { return (abs(x)/sqrt(sum(x^2))) })
RDP <- sample_dot_product(adj_mat)
RDP <- graph_from_adjacency_matrix(adj_mat)
embed <- embed_adjacency_matrix(RDP, 1)




dat <- loadDataSet("3D S Curve")
## use the S4 Class directly:
diffmap <- DiffusionMaps()
emb <- diffmap@fun(dat, diffmap@stdpars)
## simpler, use embed():
emb2 <- embed(dat, "DiffusionMaps")
plot(emb, type = "2vars")
samp <- sample(floor(nrow(dat) / 10))
embsamp <- diffmap@fun(dat[samp], diffmap@stdpars)
embother <- embsamp@apply(dat[-samp])
plot(embsamp, type = "2vars")
points(embother)



