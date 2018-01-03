library(neuralnet)
library(ggplot2)
m <- 10
x1 <- 0:m/m
y1 <- x1
y1[x1>1/2] <- 1
y1[x1<=1/2] <- 0
data1 <- data.frame(x1,y1)
x <- 0:99/99
y <- x
y[x>1/2] <- 1
y[x<=1/2] <- 0
z <- x*0
data <- data.frame(x,y,z,z,z,z,z)
for (i in 1:5){
    nn <- neuralnet(y1~x1, data=data1, hidden=c(m), act.fct="logistic", linear.output=TRUE, threshold=1e-4)
    data[,i+2] <- compute(nn, x)$net.result
}
plot <- ggplot(data=data) +
    geom_line(aes(x=x, y=y, colour = "(x>1/2)?1:0")) +
    geom_line(aes(x=x, y=data$z, colour = "ANN 1")) +
    geom_line(aes(x=x, y=data$z.1, colour = "ANN 2")) +
    geom_line(aes(x=x, y=data$z.2, colour = "ANN 3")) +
    geom_line(aes(x=x, y=data$z.3, colour = "ANN 4")) +
    geom_line(aes(x=x, y=data$z.4, colour = "ANN 5")) +
    geom_point(data=data.frame(x1,y1), aes(x=x1,y=y1)) +
    ylab("f(x)") + ylim(-.5,1.5)
ggsave("plot_ANN_Ind_Func.png", width=6, height=3)




## library(ggplot2)
## data1 <- read.csv("model_data.csv")
## jpeg(file = "sin_plot.jpg", width=4*240, height=3*240)
## x <- rep(0:9/9, 10)
## y <- rep(0:9/9, each=10)
## z <- sin(x*2*pi) + sin(y*2*pi)
## d1 <- as.matrix(cbind(data1, 1), 100)
## d2 <- as.matrix(cbind(x,y,2))
## data2 <- data.frame(rbind(d1, d2))
## names(data2) <- c("x", "y", "id")
## ggplot(data2, aes(x, y, colour=id)) + geom_point() + theme(legend.position="none")
## dev.off()

## XOR <- c(0,1,1,0)
## xor.data <- data.frame(expand.grid(c(0,1), c(0,1)), XOR)
## print(net.xor <- neuralnet( XOR~Var1+Var2, xor.data, hidden=2, rep=5))
## compute(net.xor,xor.data)
## plot(net.xor, rep="best")
## plot(nn, rep="best")

## data(infert, package="datasets")
## print(net.infert <- neuralnet(case~parity+induced+spontaneous, infert,
## err.fct="ce", linear.output=FALSE, likelihood=TRUE))
