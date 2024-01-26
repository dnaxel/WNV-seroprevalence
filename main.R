graphics.off()
rm(list = ls())
library(ggplot2)
library(cowplot)
library(deSolve)
library(mvtnorm)
library(openxlsx)
set.seed(5041)

### Color palettes for plotting
library(RColorBrewer)
year.colors <- c(brewer.pal(7,"Dark2"), "darkred", "purple4", "deeppink4")
subregion.colors <- year.colors[1:3]
background.colors <- brewer.pal(3, "Pastel2")


###############################################
### DATA LOADING & PARAMETER INIZIALIZATION ###
###############################################

# Data of captured birds, tested for WNV positivity; clustered by subregion and 2-weeks time periods 
corvid.data.by.species <- read.xlsx("corvid_data.xlsx")
corvid.data <- aggregate(corvid.data.by.species[,c("TOTAL", "WNV")],
                       by = corvid.data.by.species[,c("WEEK", "YEAR", "SUBREGION")],
                       sum)
corvid.data <- cbind(corvid.data, INCIDENCE = NA, CI1 = NA, CI2 = NA)
corvid.data[,c('INCIDENCE', 'CI1', 'CI2')] <- t(apply(corvid.data[,c('WNV', 'TOTAL')], 1, 
                                               function(x){tmp <- prop.test(x[1], x[2])
                                               c(tmp$estimate, tmp$conf.int)}))

# Data of birds for WNV antibodies (subregion 2 only)
serology.data <- read.xlsx("serology_data.xlsx")
serology.data$TIME <- serology.data$WEEK/52 + serology.data$YEAR
serology.data[,c('SEROPREVALENCE', 'CI1', 'CI2')] <- t(apply(serology.data[,c('IMMUNE', 'TOTAL')], 1, 
                                                         function(x){tmp <- prop.test(x[1], x[2])
                                                         return(c(tmp$estimate, tmp$conf.int))}))

### Parameters
C <- -2        # C*: log fraction total infected in a year (to estimate) 
m <- 35        # peak time (weeks; to estimate) 
sigma <- 3     # length parameter of the infectious period (to estimate)
gam <- 7/15.8  # Recovery rate (weeks^-1; fixed)
rho <- 0.2     # fraction of fatal WNV infections (fixed)
q <- 0.418     # decay of seroprevalence after spring due to newborns (fixed)
var.b <- 1     # prior standard deviation of log(b) (sample bias parameter; fixed)
params <- c(C, m, sigma, gam, rho)  #parameters of the ODE model

# Other variables
x0 <- c(1,0,0) # initial condition (Y=2013, w=20)
times  <- seq(20, 46, by = 0.05)    # simulation times
times2 <- seq(46.1, 72, 0.1) 
n.times <- length(times)
n.times2 <- length(times2)
n.times12 <- n.times+n.times2
years <- 2013:2022
n.years <- length(years)
subregions <- 1:3
n.subregions <- length(subregions)

############################################
### STATISTICAL & MATHEMATICAL MODELLING ###
############################################

# Obtain predicted observed incidence (y*(t)) from predicted incidence (y(t)); y: WNV incidence; b: sampling bias
apply.bias <- function(y, b)  b/(1/y + (b - 1))

# Log-likelihood
# input params: parameter vector of length 151 (50 params/subregion + b*); C*=log(C), b*=log(b)
log.likelihood <- function(params, dat, q = 0.418){ 
  
  # Boundary conditions: sigma>=0, 20<=m<=46 
  if(sum(params[seq(3, 150, 5)]<0) + sum(params[seq(2, 149, 5)]>46) + sum(params[seq(2, 149, 5)]<20)) return(-Inf)
  b <- exp(params[151])
  llk <- 0  #initialize output
  
  for(s in subregions-1){
    .x0 <- x0
    for(y in 1:n.years-1){
      # extract parameter set and data of subregion s, year y
      .params <- c(exp(params[50*s+5*y+1]), params[(2:5)+50*s+5*y])
      .dat <- dat[dat$YEAR==(2013+y) & dat$SUBREGION==(s+1),]
      
      # Simulate epidemiological model:
      sim <- as.data.frame(ode(.x0, times, SIR, .params[1:5]))
      
      # Boundary condition: x(t)>=0
      if(sum(sim$`1`<0)) return(-Inf)
      
      # Extract y*(t) at time points of interest:
      y.star <- apply.bias(sim[sim[,1]%in%(.dat$WEEK),3], b)
      y.star[y.star<0] <- 0 # some negative values are caused by deSolve approximation
      
      # Binomial log-likelihood
      llk <- llk + sum(dbinom(.dat$WNV, .dat$TOTAL, y.star, log = TRUE))
      
      # z(final_time) is required to compute next year's x0
      z.f <- sim[nrow(sim), 4]   
      .x0 <- c(1-z.f*q, 0, z.f*q)
    }
  }
  return(llk)
}

# Function employed to provide estimates and credible intervals:
quantile95 <- function(x){return(as.vector(quantile(x, c(0.5, 0.025, 0.975))))}

# ODE Model (deSolve package)
SIR <- function(t, X, param)
{
  x <- X[1];
  y <- X[2];
  z <- X[3];
  C <- param[1];
  m <- param[2];
  sigma <- param[3];
  gam <- param[4];
  rho <- param[5];
  
  dx <- -C*dnorm(t, mean=m, sd=sigma)+gam*rho*x*y;
  dy <- C*dnorm(t, mean=m, sd=sigma)-gam*y*(1-rho*y);
  dz <- gam*(1-rho)*y+gam*rho*y*z;
  
  return(list(c(dx,dy,dz)))
}

# Simulation of ODE model for all years:
sim.complete <- function(param, x0 = c(1,0,0), q = 0.418, Gam = 0.443, Rho = 0.20){
  sim.tot <- data.frame()
  for(year in years){
    sim <- as.data.frame(ode(x0, times, SIR, 
                             c(exp(param[1+3*(year-2013)]), param[2:3 + 3*(year-2013)], Gam, Rho)))
    z.f <- sim[nrow(sim),4]
    x0 <- c(1-q*z.f, 0, q*z.f)
    sim.tot <- rbind(sim.tot, sim)
  }
  return(sim.tot)
}

# Simulation of z(t) for a subregion (assuming linear decay during the breeding period)
# param is a 30-length vector of (C*, m, sigma) parameters
z.t <- function(param, Gam = 0.443, Rho = 0.20, q = 0.418, x.0 = c(1, 0, 0)){
  res <- c()
  for(year in years){
    sim <- as.data.frame(ode(x.0, times, SIR, 
                             c(exp(param[1+(year-2013)*3]), param[2:3 + (year-2013)*3], Gam, Rho)))
    x.0 <- c(1-q*sim[n.times, 4], 
             0, q*sim[n.times, 4])
    res <- c(res, sim[,4], sim[n.times, 4]*(1-(times2-46)*(1-q)/26))
  }
  return(res)
}


##########################
### POSTERIOR SAMPLING ###
##########################

# MCMC Metropolis-Hastings implementation
# params: 151-length vector of parameters (C*, m, sigma, gamma, rho) for each year and subregion + parameter b*
# dat: full data clustered by subregion and biweekly periods, spanning years 2013-2022
# varcov: 91x91 variance covariance matrix
# var.b: prior variance of b*
MCMC.MH <- function(params, dat, lgth = 50000, varcov = diag(rep(1, 91)), var.b = 1, step.size = 0.05, q = 0.418)  
{
  par.to.opt<-c(rep(seq(1, 5*n.years*n.subregions, 5), rep(3, n.years*n.subregions))+0:2, 5*n.years*n.subregions+1) #indexes of the parameters to estimate
  varcov <- varcov*step.size
  
  # compute log posterior = log prior(b,C) + log likelihood
  lpost <- dnorm(params[length(params)], 0, var.b, log = TRUE) + sum(params[seq(1, 5*n.years*n.subregions, 5)]) + 
    log.likelihood(params, dat, q = q)
  k <- 1  # accepted steps counter
  N <- 1  # total steps counter
  
  # initialize output and store I.C.
  res <- matrix(0, lgth, 3*n.years*n.subregions+3)
  res[1,] <- c(N, params[par.to.opt], lpost)
  while (N<lgth){
    N <- N+1
    
    # generate a new candidate inside the support of the posterior distribution
    lpost1 <- NA
    while(is.na(lpost1) | is.null(lpost1)){
      par1 <- params
      
      # Sample candidate
      par1[par.to.opt] <- rmvnorm(1, params[par.to.opt], varcov)   
      
      # Compute log(posterior density) of the candidate
      lpost1 <- dnorm(par1[length(par1)], 0, var.b, log = TRUE) + sum(par1[seq(1, 5*n.years*n.subregions, 5)]) + 
        log.likelihood(par1, dat, q = q)
    }
    LR <- lpost1-lpost        #log("a posteriori" density Ratio)
    if (LR > -rexp(1)){       #acceptance condition; -rexp = log(runif(0,1))
      params <- par1
      lpost <- lpost1
      k <- k+1
      show(paste(k, N, round(k/N*100), round(lpost1, 1)))
    }
    
    # Store step N:
    res[N,] <- c(N, params[par.to.opt], lpost)        
  }
  return(res)
}

# Sample the first state of the Markov chain (multi-start MCMC algorithm)
sample.initial.state <- function(data, Q = 0.418, Var.b = 1, Gam = gam, Rho = rho){
  llk <- -Inf
  # Generate initial parameters 
  while(llk == -Inf | is.null(llk)){
    #beta distribution for C improves run-time
    params <- c(log(rbeta(n.years*n.subregions, 1, 10)), 
              runif(n.years*n.subregions, 20, 46), 
              rexp(n.years*n.subregions, 0.2), 
              rep(Gam, n.years*n.subregions), 
              rep(Rho, n.years*n.subregions))
    params <- c(params[rep(1:(n.years*n.subregions), rep(5, n.years*n.subregions)) + 
                       seq(0, n.years*n.subregions*4, n.years*n.subregions)], #rearrange vector
              rnorm(1, 0, Var.b)) #sample b*
    llk <- log.likelihood(params, data, Q)
  }
  return(params)
}


### First MCMC run: used to estimate a sample variance-covariance matrix

# For this run, an arbitrary 91x91 diagonal variance-covariance matrix is employed:
prior.var <- diag(c(rep(c(.5, 5, 5), n.years*n.subregions), var.b))

# params.full is the starting point of the Markov chain
params.full <- sample.initial.state(corvid.data)
exploratory.chain <- MCMC.MH(params.full, corvid.data, lgth = 100000, varcov = prior.var, step.size = 0.001, q = q)  

# Traceplot of the log posterior density (log-posterior at convergence should be ~ -380):
# Burn-in is 50% (red vertical line)
ggplot(as.data.frame(exploratory.chain[,c(1,93)]), aes(x=V1, y=V93)) + geom_line(size=0.5)+
  theme(panel.grid = element_line(colour = "gray70"), panel.background = element_rect(fill = "gray98"), plot.caption = element_text(size = 12), plot.title=element_text(size=15, face = "bold"))+
  labs(title='log-posterior: traceplot') + xlab("Cycle") + ylab("Log-posterior density") +
  geom_vline(xintercept = 50000, color="red")

# Estimate the variance-covariance matrix
covmat <- as.matrix(cov(exploratory.chain[50001:100000, 2:92], exploratory.chain[50001:100000, 2:92]))

# MCMC run with the sample var-cov matrix:
params.full <- sample.initial.state(corvid.data)
mcmc.sample <- MCMC.MH(params.full, corvid.data, lgth = 100000, varcov = covmat, step.size = 0.033, q = q) 


################
### Figure 1 ###
################

library(ggmap)
library(osmdata)
library(raster)
library(viridis)
#library(ggnewscale)

### Get poligons for all municipalities
italy.towns <- getData("GADM",country="ITA",level=3)
emilia.towns <- italy.towns[which(italy.towns$NAME_1=="Emilia-Romagna"),]
emilia.towns$NAME_3 <- toupper(emilia.towns$NAME_3)
# towns that were recently transferred from Marche to Emilia-Romagna must be manually added:
marche.towns <- italy.towns[which(italy.towns$NAME_1=="Marche"),]
new.towns <- c("Casteldelci", "Maiolo", "Pennabilli", "Novafeltria", "San Leo", "Sant' Agata Feltria", "Talamello")
emilia.data <- rbind(fortify(emilia.towns), fortify(marche.towns[which(marche.towns$NAME_3 %in% new.towns),]))

### Get bird data clustered by municipality
corvids.by.town <- read.xlsx("corvids_towns.xlsx")
emilia.borders <- merge(corvids.by.town[,c("id", "TOTAL", "INCIDENCE", "SUBREGION")], cbind(n = 1:nrow(emilia.data), emilia.data), all=T)
emilia.borders <- emilia.borders[order(emilia.borders$n),]
emilia.borders$n_specimens[is.na(emilia.borders$TOTAL)] <- 0

### Get Map
# an API key to download Stadia maps can be obtained at https://client.stadiamaps.com/signup/

your.API.key <- "06cb6ca4-36bf-4d5a-b2be-8cc8abbdc3a5"
register_stadiamaps(key = your.API.key)
map.emilia <- get_map(location = getbb("emilia-romagna"), maptype = "stamen_terrain_background", source = "stadia")

### Figure 1A (Bubble Plot)
p1A <- ggmap(map.emilia) +
  geom_polygon(aes(x = long, y = lat, group = group, fill = SUBREGION), data = emilia.borders,
               alpha = 0.5, linewidth = 0.8, color = rgb(0,0,0,0.15)) +
  scale_fill_manual(values = subregion.colors, name = "Subregion") +
  geom_point(data = corvids.by.town, aes(x = long, y = lat, size = TOTAL, color = INCIDENCE)) +
  scale_size_continuous(range=c(1, 5), trans = "log10", name = "Sample size", breaks = c(1, 10, 100, 1000, 10000)) +
  scale_color_viridis(name = "WNV incidence", option = "plasma", trans = "sqrt") +
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank())
p1A
ggsave("Fig1A.jpeg", p1A, dpi = 1200, width = 9, height = 6)

### Figure 1B (Italy map)
italy.poligon <- fortify(getData("GADM", country="ITA", level=0))
p1B <- ggplot() +
  geom_polygon(aes(x = long, y = lat, group = group), data = italy.poligon,
               size = 0.5, color = NA, fill = "grey80", alpha = 1) +
  geom_polygon(aes(x = long, y = lat, group = group, fill = SUBREGION), data = emilia.borders,
               alpha = 1, color = NA) +
  scale_fill_manual(values = subregion.colors) +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = "white"), 
        axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(), axis.line = element_blank(),
        legend.position = "none")
p1B
ggsave("Fig1B.jpeg", p1B, dpi = 1200, width = 8, height = 11)


###############
### Table 2 ###
###############

### Corvid count
(T1 <- aggregate(corvid.data[,c("TOTAL", "WNV")], by = corvid.data[,c("YEAR", "SUBREGION")], sum))

### Count by species
aggregate(corvid.data.by.species[,c("TOTAL", "WNV")],
          by = list(corvid.data.by.species[,"LATIN_NAME"]),
          sum)


###########################
### Figure 3 / Table S1 ###
###########################

### Estimate of b
round(exp(quantile95(mcmc.sample[50000:70000,92])), 2)

### Estimate parameters from the MCMC output
get.estimated.params <- function(chain, burn.in = 50000){
  res <- data.frame(t(apply(chain[burn.in:nrow(chain), 2:91], 2, quantile95)))
  res[1+3*0:29,] <- exp(res[1+3*0:29,])
  colnames(res) <- c("ESTIMATE", "CI1", "CI2")
  res <- round(res, 2)
  res$SUBREGION <- rep(1:3, rep(30, 3))
  res$YEAR <- as.factor(rep(rep(years, rep(3, 10)), 3))
  res$PARAMETER <- rep(c("C", "m", "sigma"), 30)
  return(res)
}

est.params <- get.estimated.params(mcmc.sample)

### Table 3
cbind(round(est.params[,1:3], 2), est.params[,4:5])

### Point-range plots
# y-Axes limits
C.lim <- c(0, ceiling(max(est.params[est.params$PARAMETER=="C",1:3])/0.4)*0.4)
m.lim <- c(floor(min(est.params[est.params$PARAMETER=="m",1:3])/5)*5, 
           ceiling(max(est.params[est.params$PARAMETER=="m",1:3])/5)*5)
sigma.lim <- c(floor(min(est.params[est.params$PARAMETER=="sigma",1:3])/10)*10, 
               ceiling(max(est.params[est.params$PARAMETER=="sigma",1:3])/10)*10)

# plotting all subregions and years for each parameter
est.params$YEAR <- est.params$YEAR+0.15*(as.numeric(est.params$SUBREGION)-2)

.pC <- ggplot(data = est.params[est.params$PARAMETER=="C",]) + 
  theme(panel.grid = element_line(colour = "gray70"), panel.background = element_rect(fill = "gray99"), 
        axis.text = element_text(size = 12), axis.title=element_text(size=15),
        panel.grid.minor.x = element_blank(), panel.grid.major.x = element_line(linetype = "dotted"), panel.grid.minor.y = element_line(linetype = "dashed")) +
  geom_pointrange(aes(x=YEAR, y=ESTIMATE, ymin=CI1, ymax=CI2, color = SUBREGION), size=1.5, fatten=2, linewidth = 1) +
  geom_point(aes(x=YEAR, y=ESTIMATE), fill="black", size = 1) + scale_color_manual(values=subregion.colors) +
  labs(title="", caption="") + ylab("C") + xlab("") +
  scale_y_continuous(breaks = seq(C.lim[1], C.lim[2], 0.4), limits = C.lim) + scale_x_continuous(breaks = seq(2013, 2022), limits = c(2012.85, 2022.15)) + 
  theme(legend.position = "none")

.pM <- ggplot(data = est.params[est.params$PARAMETER=="m",]) + 
  theme(panel.grid = element_line(colour = "gray70"), panel.background = element_rect(fill = "gray99"), 
        axis.text = element_text(size = 12), axis.title=element_text(size=15),
        panel.grid.minor.x = element_blank(), panel.grid.major.x = element_line(linetype = "dotted"), panel.grid.minor.y = element_line(linetype = "dashed")) +
  geom_pointrange(aes(x=YEAR, y=ESTIMATE, ymin=CI1, ymax=CI2, color = SUBREGION), size=1.5, fatten=2, linewidth = 1) +
  geom_point(aes(x=YEAR, y=ESTIMATE), fill="black", size = 1) + scale_color_manual(values=subregion.colors) +
  labs(title="", caption="") + ylab("m") + xlab("") +
  scale_y_continuous(breaks = seq(m.lim[1], m.lim[2], 5), limits = m.lim) + scale_x_continuous(breaks = seq(2013, 2022), limits = c(2012.85, 2022.15)) + 
  theme(legend.position = "none")

.pS <- ggplot(data = est.params[est.params$PARAMETER=="sigma",]) + 
  theme(panel.grid = element_line(colour = "gray70"), panel.background = element_rect(fill = "gray99"), 
        axis.text = element_text(size = 12), axis.title=element_text(size=15),
        panel.grid.minor.x = element_blank(), panel.grid.major.x = element_line(linetype = "dotted"), panel.grid.minor.y = element_line(linetype = "dashed")) +
  geom_pointrange(aes(x=YEAR, y=ESTIMATE, ymin=CI1, ymax=CI2, color = SUBREGION), size=1.5, fatten=2, linewidth = 1) +
  geom_point(aes(x=YEAR, y=ESTIMATE), fill="black", size = 1) + scale_color_manual(values=subregion.colors) +
  xlab("Year") + ylab("\u03c3") + labs(title="", caption="") + 
  scale_y_continuous(breaks = seq(sigma.lim[1], sigma.lim[2], 10), limits = sigma.lim) + scale_x_continuous(breaks = seq(2013, 2022), limits = c(2012.85, 2022.15)) + 
  theme(legend.position = "none")

pF <- plot_grid(.pC,
                .pM,
                .pS,
                nrow = 3, ncol = 1, scale = 0.96, labels = c("A", "B", "C"))
ggsave("Fig3.jpeg", pF, height = 19.2, width = 9.6)


#####################
### Figures 2 & 4 ###
#####################

#library(ggpubr)
#library(gridExtra)

# obtain incidence + seroprevalence plots by undersampling from the MCMC output 
# res: MCMC output; burn.in: discard first n steps of MCMC; n.samples: size of the sample used to obtain the estimate
generate.outputs <- function(res, burn.in = 0, n.samples = 200, q  = 0.418, to.plot = T){
  N <- nrow(res)
  # MC estimate for all parameters
  b <- exp(median(res[(burn.in+1):N, 2+n.years*n.subregions*3]))
  par_f <- lapply(subregions, function(c){
    matrix(apply(res[1+(c-1)*3*n.years+1:(n.years*3)], 2, median), 3, n.years)
  })
  # initialize outputs
  inc.plots <- list()
  sero.plots <- list()
  fit.performance <- data.frame(YEAR = rep(years, n.subregions), 
                                SUBREGION = rep(subregions, rep(n.years, n.subregions)), 
                                N_SUCCESSES = NA, 
                                TOTAL = NA)
  zf <- list()
  y.t <- list()
  z.t <- list()
  
  opt.params <- 1:3
  idx <- sample((burn.in+1):N, n.samples, replace = FALSE)   # perform undersampling from the MCMC output
  
  for(s in subregions){
    show(paste("Computing: Subregion", s))
    corvid.data.tmp <- corvid.data[corvid.data$SUBREGION == s,]
    chain <- res[idx, c(1, 1:(n.years*n.subregions) + 1 + (s-1)*n.years*n.subregions, 92:93)]
    
    
    # dataframes containing n.samples simulations of y(t) and z(t):
    n.sim.y <- apply(chain[,2:(ncol(chain)-2)], 1, 
                     function(p){
                       sim.complete(p, q = q)[,3]
                     })
    n.sim.z <- apply(chain[,2:(ncol(chain)-2)], 1, 
                     function(p){
                       z.t(p, q = q)
                     })
    # apply bias to the observed incidence
    n.sim.y.b <- apply.bias(n.sim.y, b)
    
    # compute z(tf)
    tmp.zf <- data.frame(t(apply(n.sim.z[seq(n.times, n.times12*n.years, n.times12),], 1, quantile95)))
    colnames(tmp.zf) <- c("zf", "min", "max")
    zf[[s]] <- tmp.zf
    
    # plotting
    y.plot <- data.frame(cbind(rep(times, n.years),
                               apply(n.sim.y, 1, median),
                               t(apply(n.sim.y.b, 1, quantile95))))
    colnames(y.plot) <- c("time", "Y", "Y.b", "min", "max")
    y.t[[s]] <- y.plot
    z.plot <- data.frame(cbind(rep(c(times, times2), n.years)/52 + rep(years, rep(n.times12, n.years)),
                               t(apply(n.sim.z, 1, quantile95))))
    colnames(z.plot) <- c("time", "Z", "min", "max")
    z.t[[s]] <- z.plot
    if(to.plot){
      inc.plots[[s]] <- lapply(years, function(year){
        ggplot(data = y.plot[1:n.times + (year-2013)*n.times,]) + 
          geom_line(aes(x=time, y=Y), size=.5, color="black", linetype="dotted") + 
          geom_pointrange(data=corvid.data.tmp[corvid.data.tmp$YEAR==year,], 
                          aes(x=WEEK, y=INCIDENCE, ymin=CI1, ymax=CI2), 
                          color=subregion.colors[s], size=1, fatten=2) + 
          theme(panel.grid = element_line(colour = background.colors[s], size = 0.5), 
                panel.background = element_rect(fill = "gray100"), 
                panel.grid.minor = element_blank(),
                panel.grid.major.x = element_blank(),
                plot.caption = element_text(size = 12), plot.title=element_text(size = 15, face = "bold"),
                axis.ticks = element_blank(),
                axis.line = element_line(size = 0.75))+
          coord_cartesian(xlim = c(20,46), ylim=c(0, 0.5)) + ylab("") + xlab("") +
          geom_line(aes(x=time, y=Y.b), size=1, color="black") +
          geom_ribbon(aes(x=time, ymin=min, ymax=max), alpha=0.2, fill = subregion.colors[s])
      })
      
      p.sero <- ggplot(data = z.plot) + 
        theme(panel.grid = element_line(colour = "gray70"), panel.background = element_rect(fill = "gray99"), 
              plot.caption = element_text(size = 12), plot.title=element_text(size=15, face = "bold"),
              panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_line(linetype = "dotted")) +
        ylab("Seroprevalence") + xlab("") + geom_ribbon(aes(x=time, ymin=min, ymax=max), alpha=0.3, fill = subregion.colors[s]) +
        labs(title="", caption="") + scale_x_continuous(breaks = years+0.5, limits = c(2013.35, 2022.9), labels = paste("1st", "July", years, sep=" ")) +
        scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1))
      for(year in 2013:2023){
        p.sero <- p.sero + geom_line(data = z.plot[1:n.times+n.times12*(year-2013),], aes(x=time, y=Z), size=1, color="black") +
          geom_line(data = z.plot[(1+n.times):n.times12+n.times12*(year-2013),], aes(x=time, y=Z), size=1, color="black", linetype="dotted")
      }
      if(s == 2) p.sero <- p.sero + geom_pointrange(data=serology.data, 
                                            aes(x = TIME, y = SEROPREVALENCE, ymin = CI1, ymax = CI2), 
                                            color = "black", size = 0.5, fatten = 1)
      sero.plots[[s]] <- p.sero
    }
    
    # compute the number of predictions falling inside data's CIs
    fit.performance[fit.performance$SUBREGION == s, 3:4] <- t(sapply(years,
                           function(year){
                             weeks <- corvid.data.tmp[corvid.data.tmp$YEAR == year, "WEEK"]
                             tmpY <- y.plot[1:n.times + (year-2013)*n.times, c("time", "Y.b")] 
                             selY <- tmpY[tmpY$time %in% weeks, "Y.b"]
                             return(c(sum(corvid.data.tmp$CI1[corvid.data.tmp$YEAR == year] < selY & 
                                            selY < corvid.data.tmp$CI2[corvid.data.tmp$YEAR == year]), 
                                      length(selY)))
                           }))
  }
  return(list(z.f = zf, par.f = par_f, bias = b, indexes = idx, seroprevalence.plots = sero.plots, inc.pl = inc.plots, fit.performance = fit.performance, yt = y.t, zt = z.t))
}

outputs <- generate.outputs(as.data.frame(mcmc.sample), burn.in = 50000, n.samples = 200)

# fit performance
apply(outputs$fit.performance[,c("N_SUCCESSES", "TOTAL")], 2, sum)

# seroprevalence values at final times (z(tf))
outputs$z.f

# 10x3 incidence plots
fig2 <- plot_grid(plot_grid(outputs$inc.pl[[c(1,1)]] + ylab("Incidence (2013)"), outputs$inc.pl[[c(1,2)]] + ylab("Incidence (2014)"), outputs$inc.pl[[c(1,3)]] + ylab("Incidence (2015)"), outputs$inc.pl[[c(1,4)]] + ylab("Incidence (2016)"), outputs$inc.pl[[c(1,5)]] + ylab("Incidence (2017)"),
                          outputs$inc.pl[[c(1,6)]] + ylab("Incidence (2018)"), outputs$inc.pl[[c(1,7)]] + ylab("Incidence (2019)"), outputs$inc.pl[[c(1,8)]] + ylab("Incidence (2020)"), outputs$inc.pl[[c(1,9)]] + ylab("Incidence (2021)"), outputs$inc.pl[[c(1,10)]] + ylab("Incidence (2022)") + xlab("Time (weeks)"),
                          ncol = 1, nrow = 10),
                plot_grid(outputs$inc.pl[[c(2,1)]], outputs$inc.pl[[c(2,2)]], outputs$inc.pl[[c(2,3)]], outputs$inc.pl[[c(2,4)]], outputs$inc.pl[[c(2,5)]],
                          outputs$inc.pl[[c(2,6)]], outputs$inc.pl[[c(2,7)]], outputs$inc.pl[[c(2,8)]], outputs$inc.pl[[c(2,9)]], outputs$inc.pl[[c(2,10)]] + xlab("Time (weeks)"),
                          ncol = 1, nrow = 10),
                plot_grid(outputs$inc.pl[[c(3,1)]], outputs$inc.pl[[c(3,2)]], outputs$inc.pl[[c(3,3)]], outputs$inc.pl[[c(3,4)]], outputs$inc.pl[[c(3,5)]],
                          outputs$inc.pl[[c(3,6)]], outputs$inc.pl[[c(3,7)]], outputs$inc.pl[[c(3,8)]], outputs$inc.pl[[c(3,9)]], outputs$inc.pl[[c(3,10)]] + xlab("Time (weeks)"),
                          ncol = 1, nrow = 10),
                nrow = 1, ncol = 3, scale = 0.97)
ggsave("Fig2_2.jpeg", fig2, height = 19.2, width = 13.4)

fig4 <- plot_grid(outputs$seroprevalence.plots[[1]], 
                outputs$seroprevalence.plots[[2]], 
                outputs$seroprevalence.plots[[3]] + xlab("Time (years)"),
                nrow = 3)
ggsave("Fig4.jpeg", fig4, height = 3500, width = 3500, units = "px")

############################
### Sensitivity Analysis ###   (Figures S3, 5, S4, S5; tables S2-3-4-5-6-7)
############################

### Perturbed parameters:
gam.p <- 7*5/4/15.8 # Ti -20%
gam.m <- 7*5/6/15.8 # Ti +20%
rho.m <- 0.1 # rho -50%
rho.p <- 0.3 # rho +50%
q.m <- 0.418*3/4 # q -25%
q.p <- 0.418*5/4 # q +25%

### Run MCMC for each perturbed parameter
# Gamma_plus
params.full <- sample.initial.state(corvid.data, Gam = gam.p)
chain_gam.p <- MCMC.MH(params.full, corvid.data, lgth = 100000, varcov = covmat, step.size = 0.03, q = q)
# Gamma_minus
params.full <- sample.initial.state(corvid.data, Gam = gam.m)
chain_gam.m <- MCMC.MH(params.full, corvid.data, lgth = 100000, varcov = covmat, step.size = 0.02, q = q)
# Rho_plus
params.full <- sample.initial.state(corvid.data, Rho = rho.p)
chain_rho.p <- MCMC.MH(params.full, corvid.data, lgth = 100000, varcov = covmat, step.size = 0.033, q = q)
# Rho_minus
params.full <- sample.initial.state(corvid.data, Rho = rho.m)
chain_rho.m <- MCMC.MH(params.full, corvid.data, lgth = 100000, varcov = covmat, step.size = 0.027, q = q)
# q_plus
params.full <- sample.initial.state(corvid.data, Q = q.p)
chain_q.p <- MCMC.MH(params.full, corvid.data, lgth = 100000, varcov = covmat, step.size = 0.033, q = q.p)
# q_minus
params.full <- sample.initial.state(corvid.data, Q = q.m)
chain_q.m <- MCMC.MH(params.full, corvid.data, lgth = 100000, varcov = covmat, step.size = 0.033, q = q.m)

### This function produces a sensitivity plot of z(t) given a subregion and a single perturbed parameter 
SA.plots <- function(chain, chain.p, chain.m, gam.p = gam, gam.m = gam, rho.p = rho, rho.m = rho, q.p = q, q.m = q, n.samples = 200, subregion = 2, burn.in = 40000, colors, par.label = ""){
  
  # Estimated parameters
  par_f <- list(def = apply(chain[,2:92], 2, quantile95), 
                plus = apply(chain.p[,2:92], 2, quantile95), 
                min = apply(chain.m[,2:92], 2, quantile95))
  
  # undersampling from the MCMC output
  N <- min(nrow(chain), nrow(chain.p), nrow(chain.m))
  idx <- sample((burn.in+1):N, n.samples, replace = FALSE)
  
  # compute z(t) for all sampled parameter sets
  immune <- list(def = data.frame(apply(chain[idx, 2:31+(subregion-1)*30], 1, function(x){z.t(as.numeric(x))})),  # default parameter
                 plus = data.frame(apply(chain.p[idx, 2:31+(subregion-1)*30], 1, function(x){z.t(as.numeric(x), Gam = gam.p, Rho = rho.p, q = q.p)})),    # "plus" parameter
                 minus = data.frame(apply(chain.m[idx, 2:31+(subregion-1)*30], 1, function(x){z.t(as.numeric(x), Gam = gam.m, Rho = rho.m, q = q.m)})))   # "minus" parameter
  
  # time values for the final plot
  time.tot <- rep(c(times, times2), n.years) + 
    rep(52*0:(n.years-1), rep(n.times12, n.years))
  
  # compute z(t) median and credible intervals for all 3 conditions
  immune.tot <- as.data.frame(cbind(2013 + time.tot/52,
                                    t(apply(immune$def, 1, quantile95)),
                                    t(apply(immune$plus, 1, quantile95)),
                                    t(apply(immune$minus, 1, quantile95))))
  colnames(immune.tot) <- c("Time", "def", "min.def", "max.def", "plus", "min.plus", "max.plus", "minus", "min.minus", "max.minus")
  # Plot (3 solid/dotted lines + ribbons representing 95% credible intervals)
  SA.plot <- ggplot(immune.tot) +
    geom_ribbon(aes(x = Time, ymin = min.def, ymax = max.def), alpha = .2) +
    geom_ribbon(aes(x=Time, ymin=min.plus, ymax=max.plus), alpha = .2, fill = "red") +
    geom_ribbon(aes(x=Time, ymin=min.minus, ymax=max.minus), alpha = .2, fill = "blue") +
    ylab("Seroprevalence") + xlab("") + scale_color_manual(name = par.label, values = colors) +
    scale_x_continuous(breaks = years+0.5, limits = c(2013, 2023), labels = paste("1st July", years, sep=" ")) +
    scale_y_continuous(breaks = 0:5/5) +
    theme(
      panel.grid = element_line(colour = "gray70"),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_line(linetype = "dotted"),
      panel.background = element_rect(fill = "gray98"),
    )
  for(y in 2013:2023){
    SA.plot <- SA.plot + geom_line(data = immune.tot[1:n.times+n.times12*(y-2013),], aes(x=Time, y=def, color=labels(colors)[2]), size=1) +
      geom_line(data = immune.tot[(1+n.times):n.times12+n.times12*(y-2013),], aes(x=Time, y=def, color=labels(colors)[2]), size=1, linetype="dotted") +
      geom_line(data = immune.tot[1:n.times+n.times12*(y-2013),], aes(x=Time, y=plus, color=labels(colors)[3]), size=1) +
      geom_line(data = immune.tot[(1+n.times):n.times12+n.times12*(y-2013),], aes(x=Time, y=plus, color=labels(colors)[3]), size=1, linetype="dotted") +
      geom_line(data = immune.tot[1:n.times+n.times12*(y-2013),], aes(x=Time, y=minus, color=labels(colors)[1]), size=1) +
      geom_line(data = immune.tot[(1+n.times):n.times12+n.times12*(y-2013),], aes(x=Time, y=minus, color=labels(colors)[1]), size=1, linetype="dotted")
  }
  return(SA.plot)
}

### Generate plots for each subregion by combining each parameter's sensitivity plot
# Color associations for the legend
colors.gam <- c("0.369" = "blue", "0.443" = "black", "0.554" = "red")
colors.rho <- c("0.10" = "blue", "0.20" = "black", "0.30" = "red")
colors.q <- c("0.313" = "blue", "0.418" = "black", "0.522" = "red")

# Subregions 1, 2, 3 (figure C, 4, D)
fig_names <- c("FigS3.jpeg", "Fig5.jpeg", "FigS4.jpeg")
for(i in subregions){
  SAplot_gam <- SA.plots(chain = mcmc.sample, chain.p = chain_gam.p, chain.m = chain_gam.m, 
                          gam.p = gam.p, gam.m = gam.m, subregion = i, colors = colors.gam, par.label = "\u03B3")
  SAplot_rho <- SA.plots(chain = mcmc.sample, chain.p = chain_rho.p, chain.m = chain_rho.m, 
                         rho.p = rho.p, rho.m = rho.m, subregion = i, colors = colors.rho, par.label = "\u03C1")
  SAplot_q <- SA.plots(chain = mcmc.sample, chain.p = chain_q.p, chain.m = chain_q.m, 
                       q.p = q.p, q.m = q.m, subregion = i, colors = colors.q, par.label = "q")
  ggsave(fig_names[i],
         plot_grid(SAplot_gam, 
                   SAplot_rho, 
                   SAplot_q + xlab("Time (years)"),
                   nrow = 3),
         height = 4900, width = 3500, units = "px")
}


### Estimate mean shift of z(t0) and z(tf) in S2 after the perturbation of q
z.t_default <- c(z.t(param = as.numeric(apply(mcmc.sample[50000:nrow(mcmc.sample), 32:61], 2, median))))
z.t_q.p <- c(z.t(param = as.numeric(apply(chain_q.p[50000:nrow(chain_q.p), 32:61], 2, median)), q = q.p))
z.t_q.m <- c(z.t(param = as.numeric(apply(chain_q.m[50000:nrow(chain_q.m), 32:61], 2, median)), q = q.m))

mean_shift_t0 <- mean(c(rep(0, 2),
                        abs(z.t_q.p[1+1:9*n.times12]/z.t_default[1+1:9*n.times12])-1,
                        abs(z.t_q.m[1+1:9*n.times12]/z.t_default[1+1:9*n.times12])-1))
mean_shift_tf <- mean(c(abs(z.t_q.p[n.times+0:9*n.times12]/z.t_default[n.times+0:9*n.times12])-1,
                        abs(z.t_q.m[n.times+0:9*n.times12]/z.t_default[n.times+0:9*n.times12])-1))


### Figure E (sensitivity on x_0)
library(latex2exp)
# alternative initial condition (computed from serology.data)
x0.alt <- c(1-sum(serology.data$IMMUNE[1:4])/sum(serology.data$TOTAL[1:4]), 
            0, sum(serology.data$IMMUNE[1:4])/sum(serology.data$TOTAL[1:4]))

# sample the posterior distribution considering the alternative initial condition
x0 <- x0.alt
params.full <- sample.initial.state(corvid.data)
chain_x0 <- MCMC.MH(params.full, corvid.data, lgth = 100000, varcov = covmat, step.size = 0.03, q = q)
x0 <- c(1,0,0)

# this function produces the plot of z(t) given the standard and alternative initial conditions
SA.plot.x0 <- function(chain, chain.x0, x0.alt, n.samples = 200, cl = 2) {
  par_f <- list(def = apply(chain[, 2:92], 2, quantile95),
                alt = apply(chain.x0[, 2:92], 2, quantile95))
  N <- min(nrow(chain), nrow(chain.x0))
  par <- params
  
  # undersampling from the MCMC output(s)
  idx <- sample(1:N, n.samples, replace = FALSE)
  mcmc.data <- rbind(chain[idx, ], chain.x0[idx, ])[, 2:31 + (cl - 1) * 30]
  
  # prepare time values for the final plot
  time.tot <- rep(c(times, times2), n.years) +
    rep(52 * 0:(n.years - 1), rep(n.times12, n.years))
  immune <- list(def = data.frame(apply(mcmc.data[1:n.samples, ], 1,
                                        function(x){
                                          tmp.params <- as.numeric(x)
                                          z.t(tmp.params)
                                        })),
                 alt = data.frame(apply(mcmc.data[1:n.samples + n.samples, ], 1,
                                        function(x){
                                          tmp.params.alt <- as.numeric(x)
                                          z.t(tmp.params.alt, x.0 = x0.alt)
                                        })))
  
  
  immune.tot <- as.data.frame(cbind(2013 + time.tot/52,
                                    t(apply(immune$def, 1, quantile95)),
                                    t(apply(immune$alt, 1, quantile95))))
  colnames(immune.tot) <- c("Time", "def", "min.def", "max.def", "alt", "min.alt", "max.alt")
  colors <- c("1" = "black", "0.55" = "blue")
  
  #plot the total seroprevalence prediction over 10 years
  SAplot <- ggplot(immune.tot) +
    theme(
      panel.grid = element_line(colour = "gray70"),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_line(linetype = "dotted"),
      panel.background = element_rect(fill = "gray98"),
    ) +
    geom_ribbon(aes(x = Time, ymin = min.def, ymax = max.def), alpha = 0.2) +
    geom_ribbon(aes(x = Time, ymin = min.alt, ymax = max.alt), alpha = 0.2, fill = "blue") +
    scale_x_continuous(
      breaks = years + 0.5,
      limits = c(2013, 2023),
      labels = paste("1st", "July", years, sep = " ")) +
    labs(x = "Time", y = "Seroprevalence") +
    scale_color_manual(values = colors) +
    scale_y_continuous(breaks = 0:5/5)
  for (y in 2013:2023)
    SAplot <-
    SAplot + geom_line(data = immune.tot[1:n.times + n.times12 * (y - 2013), ],
                   aes(x = Time, y = alt, color = "0.55"), size = 0.8) +
    geom_line(data = immune.tot[(1 + n.times):n.times12 + n.times12 * (y - 2013), ],
      aes(x = Time, y = alt, color = "0.55"),
      size = 0.8, linetype = "dashed") +
    geom_line(data = immune.tot[1:n.times + n.times12 * (y - 2013), ],
              aes(x = Time, y = def, color = "1"), size = 0.8) +
    geom_line(data = immune.tot[(1 + n.times):n.times12 + n.times12 * (y - 2013), ],
      aes(x = Time, y = def, color = "1"), size = 0.8, linetype = "dashed")
  
  return(SAplot)
}

SAplot_x0 <- SA.plot.x0(chain = mcmc.sample, chain.x0 = chain_x0, x0.alt = x0.alt)
ggsave("FigE.jpeg", SAplot_x0 + labs(color = TeX(r'($x_{2013,S2}(20)$)')), height = 6, width = 15)


### Tables of estimated parameters after parameter perturbation

S1 <- get.estimated.params(chain_gam.m)
S2 <- get.estimated.params(chain_gam.p)
S3 <- get.estimated.params(chain_rho.m)
S4 <- get.estimated.params(chain_rho.p)
S5 <- get.estimated.params(chain_q.m)
S6 <- get.estimated.params(chain_q.p)
S7 <- get.estimated.params(chain_x0)


### Table 3: Estimated sampling biases

bSA <- rbind(exp(quantile95(chain_gam.m[50000:nrow(chain_gam.m),92])),
             exp(quantile95(chain_gam.p[50000:nrow(chain_gam.p),92])),
             exp(quantile95(chain_rho.m[50000:nrow(chain_rho.m),92])),
             exp(quantile95(chain_rho.p[50000:nrow(chain_rho.p),92])),
             exp(quantile95(chain_q.m[50000:nrow(chain_q.m),92])),
             exp(quantile95(chain_q.p[50000:nrow(chain_q.p),92])),
             exp(quantile95(chain_x0[50000:nrow(chain_x0),92])))


#########################
### CORRELATION PLOTS ### 
#########################

### Figure A: correlation between early infected corvids & human cases in Emilia-Romagna

# This function produces a linear regression plot given corvid incidence (x) and human cases (y) data 
plotLR <- function(x, y, names = c("x", "y")){
  Data <- data.frame(YEAR = as.factor(years), X = x, Y = y)
  corr <- lm(Y ~ X, data = Data)
  x.lin <- seq(0, 1, 0.01)*max(x)*1.1
  regr.line <- data.frame(X = x.lin, Y = corr$coefficients[1]+x.lin*corr$coefficients[2])
  
  p1 <- ggplot(Data) + geom_point(aes(x=X, y=Y, colour=YEAR), size=2.5) +
    geom_line(data = regr.line, aes(X, Y), linetype = "dotted", size = 0.7) + 
    theme(panel.grid = element_line(colour = "gray80"), panel.background = element_rect(fill = "gray98"))+
    scale_color_discrete(type = year.colors) + 
    xlab(names[1]) + ylab(names[2])
  return(p1)
}

est.params.corr <- est.params[,1]
dim(est.params.corr) <- c(3, 30)

# Early infected corvids (up to week 30)
# computed as integral from t=20 to t=30 of lambda(t)*x(t)
early.inf <- apply(est.params.corr, 2, function(x){
  x <- as.numeric(x)
  return(x[1]*(pnorm(30, x[2], x[3])-pnorm(20, x[2], x[3])))})
early.inf <- early.inf[1:10] + early.inf[11:20] + early.inf[21:30]

# Human WNV cases in Emilia-Romagna for each year 
# (reference: https://www.epicentro.iss.it/westnile/bollettino)
humanCasesER <- c(16, 7, 17, 20, 10, 100, 4, 4, 18, 69)

pA <- plotLR(early.inf,
             humanCasesER,
             names = c("Early infected corvids", "Human cases"))
ggsave("FigA.jpeg", pA, dpi=1200, width = 12, height = 7)


### Figure B: correlation between observed WNV incidence in mosquitoes and corvids

# load mosquito data 
# (reference: Marini G et al (2020), "A quantitative comparison of West Nile virus incidence from 2013 to 2018 in Emilia-Romagna, Italy")
mosquito.data <- read.xlsx("mosquito_data.xlsx")
mosquito.data$INCIDENCE <- mosquito.data$WNV/mosquito.data$TOTAL #add observed incidence

# get cumulative corvid data (without subregions)
corvid.data.full <- aggregate(corvid.data[,c(4,5)], by = corvid.data[,c("WEEK", "YEAR")], sum)
corvid.data.full$INCIDENCE <- corvid.data.full$WNV/corvid.data.full$TOTAL

# merge data
lm.data <- merge(mosquito.data,
                      corvid.data.full,
                      by=c("YEAR", "WEEK"))
lm.data$YEAR <- as.factor(lm.data$YEAR)

# linear regression
mosquito.lm <- lm(INCIDENCE.y ~ INCIDENCE.x, data = lm.data)
x <- seq(0, 0.3, 0.01)
regr.line <- data.frame(x = x, y = mosquito.lm$coefficients[1]+x*mosquito.lm$coefficients[2])

# regression plot
pB <- ggplot(lm.data) + geom_point(aes(x=INCIDENCE.x, y=INCIDENCE.y, colour=YEAR), size=1.5) +
  geom_line(data = regr.line, aes(x, y), linetype = "dotted", size = 0.7) + 
  theme(panel.grid = element_line(colour = "gray80"), panel.background = element_rect(fill = "gray98"), plot.caption = element_text(size = 12), plot.title=element_text(size=15, face = "bold"))+
  scale_color_discrete(type = year.colors) + 
  ylab("Corvid incidence") + xlab("Mosquito incidence")
ggsave("FigB.jpeg", pB, dpi=1200, width = 12, height = 7)
