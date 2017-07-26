library(foreach)
library(doParallel)
data <- read.csv("model_data.csv", header=TRUE)
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)
Nmax <- dim(data)[1]/22500
## for(i in 1:Nmax){
out <- foreach(i=1:Nmax) %dopar% {
    library(plot3D)
    subdata <- data[((i-1)*22500 + 1):(i*22500),]
    attach(subdata)
    name <- paste("figs/plot", i, ".jpg", sep = "")
    jpeg(file = name, width=3*480, height=2*480)
    par(mfrow=c(3, 3))
    scatter3D(x, y, v, theta=-45)
    plotdev(theta=0)
    plotdev(theta=-180)
    scatter3D(x, y, p, theta=-45)
    plotdev(theta=0)
    plotdev(theta=-180)
    scatter3D(x, y, inv, theta=-45)
    plotdev(theta=0)
    plotdev(theta=-180)
    dev.off()
    detach(subdata)
    Sys.sleep(.5)
}
stopCluster(cl)
