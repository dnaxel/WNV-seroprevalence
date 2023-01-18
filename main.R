graphics.off()
rm(list = ls())
library(ggplot2)
library(RColorBrewer)
library(deSolve)
library(mvtnorm)

#color palettes
colback <- brewer.pal(6,"Pastel2")
collin <- brewer.pal(6,"Set2")
colconf <- brewer.pal(6,"Dark2")
palettes <- data.frame(brewer.pal(6,"BuGn"), brewer.pal(6,"YlOrRd"), brewer.pal(6,"GnBu"), brewer.pal(6,"PuRd"), brewer.pal(6,"YlGn"), brewer.pal(6,"YlOrRd"))
colcorr <- brewer.pal(11, "RdYlBu")

#parameters
C <- -1;  #log frac total infected 
m <- 33;  #peak time (weeks)
sigma <- 3;
gam <- 0.44;
rho <- 0.1; #% dead birds due to WNV

params <- c(C, m, sigma, gam, rho)  #!order is important!
x0 <- c(1,0,0) #0 initial cases
times  <- seq(20, 47, by=0.05) 
parstot1 <- rep(params,6)

######################################
##### Variance covariance matrix #####
######################################

params <- c(C, m, sigma, gam, rho)  #!order is important!
x0 <- c(1,0,0) #0 initial cases
times  <- seq(20, 47, by=0.05) 
parstot1 <- rep(params,6)

#perform an exploratory sample of the posterior to estimate a varcov matrix for the parameters
prior.var <- diag(c(rep(c(1,5^0.5,108^0.5),6),1))
res.expl <- MCMC.varcov(parstot, bird.data, lgth = 40000, varcov = prior.var, step_size = 0.05)
#compute the varcov matrix
tmp <- t(res.expl[2:20, 30001:40000]) #!check dimensions!
covmat <- cov(tmp, tmp)

##############################
##### POSTERIOR SAMPLING #####
##############################

res_multi <- MCMC.multistart(n.chain = 10, N = 10000, burn.in = 4000, covmat = as.matrix(covmat), step_size = 0.5, var.a = 1, gamma = 0.441)
#plot density estimations for the llk to keep chains only at convergence
llks <- matrix(0, 60000, 2)
llks[,2] <- res_multi$output[,21]
for(j in 1:10){
  llks[(1+6000*(j-1)):(6000*j),1] <- j
}
llks <- as.data.frame(llks)
llks[,1] <- as.factor(llks[,1])
colnames(llks) <- c("Chain", "logpost")
p1<-ggplot(llks, aes(x=logpost, fill=Chain)) + geom_density(alpha=0.3, size=1) +
  theme(panel.grid = element_line(colour = "gray70"), panel.background = element_rect(fill = "gray98"), plot.caption = element_text(size = 12), plot.title=element_text(size=15, face = "bold"))+
  labs(title='Chain-dependent MCMC convergence')
ggsave( "multi_chain_density.png", p1, dpi=1200, height = 7, width = 11)
