## library(foreach)
## library(doParallel)
data <- read.csv("../model_data.csv", header=TRUE)
library(deldir)
library(rgl)
library(png)
## data <- read.csv("model_data.csv", header=TRUE)
## cores=detectCores()
## cl <- makeCluster(cores[1]-1) #not to overload your computer
## registerDoParallel(cl)
Nmax <- dim(data)[1]/22500
mat1 <- matrix(c(.7, -.7, 0, 0,
               .3, .3, .9, 0,
               -.7, -.7, .4, 0,
               0, 0, 0, 1), nrow=4, ncol=4, byrow=TRUE)
mat2 <- matrix(c(1, 0, 0, 0,
               0, .4, .9, 0,
               0, -.9, .4, 0,
               0, 0, 0, 1), nrow=4, ncol=4, byrow=TRUE)
mat3 <- matrix(c(-1, 0, 0, 0,
               0, -.4, .9, 0,
               0, .9, .4, 0,
               0, 0, 0, 1), nrow=4, ncol=4, byrow=TRUE)
## out <- foreach(i=1:Nmax) %dopar% {
for(i in 1:Nmax){v
    library(plot3D)
    library(neuralnet)
    subdata <- data[((i-1)*22500 + 1):(i*22500),]
    subdata2 <- subdata[subdata$x1>0,-1:-3]
    ## nn <- neuralnet(inv~x1+y1, subdata2, hidden=c(36),
    ## algorithm="slr", stepmax = 1e+05, act.fct="tanh")
    ## nn.predict <- compute(nn, subdata[,c("x", "y")])$net.result
    name <- paste("../figs/plot", i, ".jpg", sep = "")
    ## name <- paste("figs/plot", i, ".jpg", sep = "")
    jpeg(file = name, width=3*480, height=1*480)
    par(mfrow=c(1, 3))
    scatter3D(subdata$x, subdata$y, subdata$v, theta=-45, zlim=c(0,300))
    plotdev(theta=0)
    plotdev(theta=-180)

    par3d(windowRect = 50 + c(0, 0, 640, 640))
    plot3d(deldir(subdata2$x1, subdata2$y1, z = subdata2$inv), zlim=c(0,4))
    rgl.viewpoint(scale=c(1, 1, 18/4), userMatrix=mat1)
    name <- paste("../figs/inv_plot_1_", i, ".png", sep = "")
    snapshot3d(name)
    rgl.viewpoint(scale=c(1, 1, 18/4), userMatrix=mat2)
    name <- paste("../figs/inv_plot_2_", i, ".png", sep = "")
    snapshot3d(name)
    rgl.viewpoint(scale=c(1, 1, 18/4), userMatrix=mat3)
    name <- paste("../figs/inv_plot_3_", i, ".png", sep = "")
    snapshot3d(name)

    par3d(windowRect = 50 + c(0, 0, 640, 640))
    plot3d(deldir(subdata2$x1, subdata2$y1, z = subdata2$p), zlim=c(0,4))
    rgl.viewpoint(scale=c(1, 1, 18/4), userMatrix=mat1)
    name <- paste("../figs/price_plot_1_", i, ".png", sep = "")
    snapshot3d(name)
    rgl.viewpoint(scale=c(1, 1, 18/4), userMatrix=mat2)
    name <- paste("../figs/price_plot_2_", i, ".png", sep = "")
    snapshot3d(name)
    rgl.viewpoint(scale=c(1, 1, 18/4), userMatrix=mat3)
    name <- paste("../figs/price_plot_3_", i, ".png", sep = "")
    snapshot3d(name)

    ## scatter3D(subdata$x, subdata$y, nn.predict, theta=-45, zlim=c(0, 4))
    ## plotdev(theta=0)
    ## plotdev(theta=-180)
    ## scatter3D(x, y, p, theta=-45)
    ## plotdev(theta=0)
    ## plotdev(theta=-180)
    ## scatter3D(subdata2$x1, subdata2$y1, subdata2$inv, theta=-45, zlim=c(0, 4))
    ## plotdev(theta=0)
    ## plotdev(theta=-180)
    dev.off()
}
## stopCluster(cl)
