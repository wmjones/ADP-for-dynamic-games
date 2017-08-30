library(ggplot2)
x <- 0:500/500*17+1
g1 <- function(x){
    if(x<=5)
        return(3*x-4)
    else
        return(12+log(2-exp(16-3*x)))
}
Discrete_g <- sapply(x,g1)
Discrete_g[is.nan(Discrete_g)] <- 0
g2 <- function(x){
    14.7/(1+exp(4-1.2*x))-2
}
Continuous_g <- sapply(x,g2)
data <- data.frame(x,Discrete_g,Continuous_g)
plot <- ggplot(aes(x=x), data=data) + geom_line(aes(y=Discrete_g, colour="Discrete_g")) + geom_line(aes(y=Continuous_g, colour = "Continuous_g")) + xlab(expression(omega)) + ylab(expression(g(omega)))
ggsave("plot_g.png", width=6,height=4)
