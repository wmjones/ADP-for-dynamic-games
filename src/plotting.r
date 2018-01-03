library(foreach)
library(doParallel)

path <- ""
data <- read.csv(paste(path, "model_data.csv", sep=""), header=TRUE)

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

col <-  rainbow(100)
color_pallet <- function(x){
    for(i in 1:99){
        if(y[[i+1]][1]>x & x>=y[[i]][1]){
            return(col[i])
        }
    }
    return(col[100])
}

tmp <- data[data$p!=-999,]
v.split <- split(data$v[order(data$v)], ceiling(seq_along(data$v)/(length(data$v)/100)))
y <- v.split
v.colors <- unlist(sapply(data$v, color_pallet))
p.split <- split(tmp$p[order(tmp$p)], ceiling(seq_along(tmp$p)/(length(tmp$p)/100)))
inv.split <- split(tmp$inv[order(tmp$inv)], ceiling(seq_along(tmp$inv)/(length(tmp$inv)/100)))
data <- data.frame(data, v.colors)
v.min <- min(data$v)
v.max <- max(data$v)
inv.min <- min(tmp$inv)-.01
inv.max <- max(tmp$inv)
p.min <- min(tmp$p)-.0001
p.max <- max(tmp$p)+.0001

cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

out <- foreach(i=1:Nmax) %dopar% {
## for(i in 1:Nmax){
    ## library(plot3D)
    library(rgl)
    library(deldir)
    library(png)
    subdata <- data[((i-1)*22500 + 1):(i*22500),]
    subdata2 <- subdata[subdata$x1>0,-1:-3]
    y <- p.split
    p.colors <- unlist(sapply(subdata2$p, color_pallet))
    y <- inv.split
    inv.colors <- unlist(sapply(subdata2$inv, color_pallet))

    name <- paste(path, "figs/", i, "_plot.jpeg", sep = "")
    jpeg(file = name, width=1440, height=1440)
    par(mfrow=c(3, 3), mai=c(0.1,0.1,0.1,0.1))

    name <- paste(path, "figs/", i, "_1_val_plot.png", sep = "")
    par3d(windowRect = c(0, 0, 480, 480))
    plot3d(subdata$x, subdata$y, subdata$v, type="p", col=subdata$v.colors,
           zlim=c(v.min, v.max), box=FALSE, xlab="", ylab="", zlab="", lit=FALSE)
    rgl.viewpoint(scale=c(1, 1, 17/(v.max-v.min)), userMatrix=mat1)
    snapshot3d(name)
    ## rgl.postscript(paste(path, "figs/val_plot_1_", i, ".svg", sep = ""), fmt="svg")
    ## img <- rsvg(paste(path, "figs/val_plot_1_", i, ".svg", sep = ""), width=200, height=200)
    img <- readPNG(name)
    frame()
    plot.window(0:1, 0:1)
    usr<-par("usr")
    rasterImage(img, usr[1], usr[3], usr[2], usr[4])
    ## rasterImage(img, 0, 1, 0, 1)

    name <- paste(path, "figs/", i, "_2_val_plot.png", sep = "")
    rgl.viewpoint(scale=c(1, 1, 17/(v.max-v.min)), userMatrix=mat2)
    ## plot3d(subdata$x, subdata$y, subdata$v, type="p", col=subdata$v.colors,
    ##        xlab="S_i", ylab="S_{-i}", zlab="Value", box=FALSE)
    snapshot3d(name)
    img <- readPNG(name)
    frame()
    plot.window(0:1, 0:1)
    usr<-par("usr")
    rasterImage(img, usr[1], usr[3], usr[2], usr[4])
    ## rasterImage(img, 0, 1, 0, 1)

    name <- paste(path, "figs/", i, "_3_val_plot.png", sep = "")
    rgl.viewpoint(scale=c(1, 1, 17/(v.max-v.min)), userMatrix=mat3)
    ## plot3d(subdata$x, subdata$y, subdata$v, type="p", col=subdata$v.colors,
    ##        xlab="S_i", ylab="S_{-i}", zlab="Value", box=FALSE)
    snapshot3d(name)
    img <- readPNG(name)
    frame()
    plot.window(0:1, 0:1)
    usr<-par("usr")
    rasterImage(img, usr[1], usr[3], usr[2], usr[4])
    ## rasterImage(img, 0, 1, 0, 1)
    clear3d()

    ## name <- paste(path, "figs/", i, "_1_val_plot.png", sep = "")
    ## par3d(windowRect = c(0, 0, 480, 480))
    ## plot3d(subdata$x, subdata$y, sin((subdata$x-min(subdata$x))/(max(subdata$x)-min(subdata$x))*2*pi)+sin((subdata$y-min(subdata$y))/(max(subdata$y)-min(subdata$y))*2*pi), type="p", col=subdata$v.colors,
    ##        zlim=c(v.min, v.max), box=FALSE, xlab="", ylab="", zlab="", lit=FALSE)
    ## rgl.viewpoint(scale=c(1, 1, 17/(v.max-v.min)), userMatrix=mat1)
    ## snapshot3d(name)
    ## ## rgl.postscript(paste(path, "figs/val_plot_1_", i, ".svg", sep = ""), fmt="svg")
    ## ## img <- rsvg(paste(path, "figs/val_plot_1_", i, ".svg", sep = ""), width=200, height=200)
    ## img <- readPNG(name)
    ## frame()
    ## plot.window(0:1, 0:1)
    ## usr<-par("usr")
    ## rasterImage(img, usr[1], usr[3], usr[2], usr[4])
    ## ## rasterImage(img, 0, 1, 0, 1)

    ## name <- paste(path, "figs/", i, "_2_val_plot.png", sep = "")
    ## rgl.viewpoint(scale=c(1, 1, 17/(v.max-v.min)), userMatrix=mat2)
    ## ## plot3d(subdata$x, subdata$y, subdata$v, type="p", col=subdata$v.colors,
    ## ##        xlab="S_i", ylab="S_{-i}", zlab="Value", box=FALSE)
    ## snapshot3d(name)
    ## img <- readPNG(name)
    ## frame()
    ## plot.window(0:1, 0:1)
    ## usr<-par("usr")
    ## rasterImage(img, usr[1], usr[3], usr[2], usr[4])
    ## ## rasterImage(img, 0, 1, 0, 1)

    ## name <- paste(path, "figs/", i, "_3_val_plot.png", sep = "")
    ## rgl.viewpoint(scale=c(1, 1, 17/(v.max-v.min)), userMatrix=mat3)
    ## ## plot3d(subdata$x, subdata$y, subdata$v, type="p", col=subdata$v.colors,
    ## ##        xlab="S_i", ylab="S_{-i}", zlab="Value", box=FALSE)
    ## snapshot3d(name)
    ## img <- readPNG(name)
    ## frame()
    ## plot.window(0:1, 0:1)
    ## usr<-par("usr")
    ## rasterImage(img, usr[1], usr[3], usr[2], usr[4])
    ## ## rasterImage(img, 0, 1, 0, 1)
    ## clear3d()

    ## scatter3D(subdata2$x, subdata2$y, subdata2$inv, theta=-45, zlim=c(inv.min, inv.max))
    ## plotdev(theta=0)
    ## plotdev(theta=-180)

    par3d(windowRect = c(0, 0, 480, 480))
    persp3d(as.numeric(levels(as.factor(subdata2$x1))),
            as.numeric(levels(as.factor(subdata2$y1))),
            matrix(subdata2$inv, nrow=length(levels(as.factor(subdata2$x1)))),
            facets=FALSE, col=inv.colors, lit=FALSE, box=FALSE,
            xlab="", ylab="", zlab="", zlim=c(inv.min, inv.max))
    ## plot3d(deldir(subdata2$x1, subdata2$y1, z = subdata2$inv), col=inv.colors, zlab="", xlab="", ylab="", zlim=c(inv.min, inv.max), box=FALSE, lit=FALSE)
    ## col <- cm.colors(20)[1 + round(19*(subdata2$inv - min(subdata2$inv))/diff(range(subdata2$inv)))]
    ## persp3d(deldir(subdata2$x1, subdata2$y1, z = subdata2$inv), col=inv.colors, lit=FALSE)
    rgl.viewpoint(scale=c(1, 1, 17/(inv.max-inv.min)), userMatrix=mat1)
    name <- paste(path, "figs/", i, "_1_inv_plot.png", sep = "")
    snapshot3d(name)
    img <- readPNG(name)
    frame()
    plot.window(0:1, 0:1)
    usr<-par("usr")
    rasterImage(img, usr[1], usr[3], usr[2], usr[4])

    rgl.viewpoint(scale=c(1, 1, 18/4), userMatrix=mat2)
    name <- paste(path, "figs/", i, "_2_inv_plot.png", sep = "")
    snapshot3d(name)
    img <- readPNG(name)
    frame()
    plot.window(0:1, 0:1)
    usr<-par("usr")
    rasterImage(img, usr[1], usr[3], usr[2], usr[4])

    rgl.viewpoint(scale=c(1, 1, 17/(inv.max-inv.min)), userMatrix=mat3)
    name <- paste(path, "figs/", i, "_3_inv_plot.png", sep = "")
    snapshot3d(name)
    img <- readPNG(name)
    frame()
    plot.window(0:1, 0:1)
    usr<-par("usr")
    rasterImage(img, usr[1], usr[3], usr[2], usr[4])
    clear3d()

    par3d(windowRect = c(0, 0, 480, 480))
    persp3d(as.numeric(levels(as.factor(subdata2$x1))),
            as.numeric(levels(as.factor(subdata2$y1))),
            matrix(subdata2$p, nrow=length(levels(as.factor(subdata2$x1)))),
            facets=FALSE, col=p.colors, lit=FALSE, box=FALSE,
            xlab="", ylab="", zlab="", zlim=c(p.min, p.max))
    ## plot3d(deldir(subdata2$x1, subdata2$y1, z = subdata2$p), col=p.colors, zlab="Price", box=FALSE, zlim=c(p.min, p.max), lit=FALSE, xlab="", ylab="", zlab="")
    rgl.viewpoint(scale=c(1, 1, 17/(p.max-p.min)), userMatrix=mat1)
    name <- paste(path, "figs/", i, "_1_price_plot.png", sep = "")
    snapshot3d(name)
    img <- readPNG(name)
    frame()
    plot.window(0:1, 0:1)
    usr<-par("usr")
    rasterImage(img, usr[1], usr[3], usr[2], usr[4])

    rgl.viewpoint(scale=c(1, 1, 17/(p.max-p.min)), userMatrix=mat2)
    name <- paste(path, "figs/", i, "_2_price_plot.png", sep = "")
    snapshot3d(name)
    img <- readPNG(name)
    frame()
    plot.window(0:1, 0:1)
    usr<-par("usr")
    rasterImage(img, usr[1], usr[3], usr[2], usr[4])

    rgl.viewpoint(scale=c(1, 1, 17/(p.max-p.min)), userMatrix=mat3)
    name <- paste(path, "figs/", i, "_3_price_plot.png", sep = "")
    snapshot3d(name)
    img <- readPNG(name)
    frame()
    plot.window(0:1, 0:1)
    usr<-par("usr")
    rasterImage(img, usr[1], usr[3], usr[2], usr[4])
    clear3d()

    dev.off()
}
## stopCluster(cl)

    ## library(neuralnet)
    ## nn <- neuralnet(inv~x1+y1, subdata2, hidden=c(36),
    ## algorithm="slr", stepmax = 1e+05, act.fct="tanh")
    ## nn.predict <- compute(nn, subdata[,c("x", "y")])$net.result
    ## scatter3D(subdata$x, subdata$y, nn.predict, theta=-45, zlim=c(0, 4))
    ## plotdev(theta=0)
    ## plotdev(theta=-180)
    ## scatter3D(x, y, p, theta=-45)
    ## plotdev(theta=0)
    ## plotdev(theta=-180)
    ## scatter3D(subdata2$x1, subdata2$y1, subdata2$inv, theta=-45, zlim=c(0, 4))
    ## plotdev(theta=0)
    ## plotdev(theta=-180)
    ## surface3d(as.numeric(levels(as.factor(subdata$x1))),
    ##           as.numeric(levels(as.factor(subdata$y1))),
    ##           t(matrix(subdata2$p, nrow=5, byrow=TRUE)))
