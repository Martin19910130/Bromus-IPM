##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##    Lambda calculation, Bootstrap
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())
gc()
options(stringsAsFactors = F)

require(dplyr)
require(ggplot2)
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###         Load functions
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  1. function, size range
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Transforms all values below/above limits in min/max size
#this is to deal with unintentional eviction
x_range <- function(x, pars){
  pmin(pmax(x, pars$L), pars$U)
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  2. inverse logit function
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Create the function to apply the inverse logit
inv_logit <- function(x){ exp(x)/(1+exp(x)) } 

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  3. Growth function
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Growth (transition) from size x to size y
gxy <- function(y,x,pars){
  xb <- x_range(x, pars)
  # returns a *probability density distribution* for each x value
  return(dnorm(y, 
               mean = pars$grow_b0 + pars$grow_b1*xb, 
               sd   = pars$grow_sd) )
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  4. Survival at size x function
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sx<-function(x,pars){
  xb <- x_range(x, pars)
  # survival prob. of each x size class 
  return( inv_logit(pars$surv_b0 + pars$surv_b1 * xb) )
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  5. transition function
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# transition: Survival * growth
pxy<-function(y,x,pars){
  xb <- x_range(x, pars)
  return( sx(xb,pars) * gxy(y,xb,pars) )
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 6. Production of flowers
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fx1<-function(x,pars){
  # probability of flowering of each x size class 
  return( inv_logit(pars$flowprob_b0 + 
                      pars$flowprob_b1 * x) )
}

fx2<-function(x,pars){
  #number of flowers per flowering plant based on size class
  return( exp( pars$seednum_b0 + 
                 pars$seednum_b1 * x) )
}

#put them together- probability of flowering * flowers per flowering plant
fx3<-function(x,pars){
  xb <- x_range(x, pars)
  return( fx1(xb,pars) * fx2(xb,pars) )
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 7. Size distribution of new plants
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
newpl <- function(y, pars){
  dnorm( x    = y,
         mean = pars$new_size_mean,
         sd   = pars$new_size_sd )
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Unintentional eviction using Merow's MEE article suggestion -----------
#             Kernel function
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ker_mee <- function(pars){
  
  # set up necessary variables
  h   <- (pars$U-pars$L)/pars$n             #Bin size
  b   <- pars$L+c(0:pars$n)*h               #Lower boundaries of bins 
  y   <- 0.5*(b[1:pars$n]+b[2:(pars$n+1)])  #Bins' midpoints
  n   <- pars$n
  
  # Fertility matrix
  Fmat            <- matrix(0,(n+1),(n+1))
  
  # April seedlings go in top row
  Fmat[1,2:(n+1)] <- fx3(y,pars) * pars$aprilSL_per_fl
  
  
  # November seedlings go into 50 by 50 plant matrix
  Fmat[2:(n+1),2:(n+1)] <- fx3(y,pars) * pars$novSL_per_fl * pars$SL_surv * newpl(y,pars) * h 
  
  # Growth/survival transition matrix
  Smat            <- matrix(0,(n+1),(n+1))
  Gmat            <- matrix(0,(n+1),(n+1))
  Tmat            <- matrix(0,(n+1),(n+1))
  
  # survival matrix
  Smat[2:(n+1),1] <- pars$SL_surv
  for(i in 2:(n+1)){ 
    Smat[2:(n+1),i] <- sx(y[i-1], pars)
  }
  
  # Growth transitions among cts sizes
  Gmat[2:(n+1),
       2:(n+1)]   <- outer(y,y,gxy,pars) * h
  
  # seedlings stage distribution
  Gmat[2:(n+1),1] <- newpl(y,pars) * h  
  
  # Growth/survival transitions among cts sizes
  Tmat            <- Gmat
  
  # fix eviction for seedligns
  Gmat[2,1]  <- Gmat[2,1]+1-sum(Gmat[2:(n+1),1])
  Tmat[,1]   <- Gmat[,1]*Smat[,1]
  
  # fix eviction of large adults
  for(i in 2:(n+1)){
    Gmat[n,i]  <- Gmat[n,i]+1-sum(Gmat[2:(n+1),i])
    Tmat[,i]   <- Gmat[,i]*Smat[,i]
  }
  
  # Full Kernel is simply a summation of fertility and transition matrix
  k_yx            <- Fmat + Tmat     
  
  return(list(k_yx    = k_yx,
              Fmat    = Fmat,
              Tmat    = Tmat,
              Gmat    = Gmat,
              meshpts = y ))
  
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## bootstrap function
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
boot_lam <- function(ii)
{
  # set the random number generator
  set.seed(ii)
  
  # sample the data
  data <- sample_n(demo_dat, nrow(demo_dat), replace = T)
  data_s <- sample_n(seed_dat, nrow(seed_dat), replace = T)

  # calculate the mean seeds 
  mean_seed <- mean(data_s$number_of_seeds, na.rm = T)
  
  # add the mean seeds to the demography data
  data$mean_seed <- ifelse(data$sample_year == 2018, mean_seed, NA)
  data$seed_per_ind <- data$mean_seed * data$number_of_flowers
  
  demo_dat$mean_seed <- ifelse(demo_dat$sample_year == 2018, mean_seed, NA)
  demo_dat$seed_per_ind <- demo_dat$mean_seed * demo_dat$number_of_flowers
  
  #--------------------------------------
  # survival
  #--------------------------------------
  sr_mod <- glm(survival ~ logsizet0, data = data, family = binomial())
  
  #--------------------------------------
  # growth
  #--------------------------------------
  gr_mod <- lm(logsizet1 ~ logsizet0, data = data)
  
  #--------------------------------------
  # flower probability
  #--------------------------------------
  fl_mod <- glm(flower ~ logsizet0, data = data, family = binomial())
  
  #--------------------------------------
  # number of seeds per flowering plant
  #--------------------------------------
  data2 <- subset(data, flower == "1")
  sefl_mod <- glm(round(seed_per_ind) ~ logsizet0, data = data2, family='poisson')
  
  #--------------------------------------
  # seedling per seed
  #--------------------------------------
  data3 <- subset(demo_dat, subplot <= 3)
  
  # get the seedling count per subplot
  Sl_per_subplot_apr18 <- unique(data3[, c("plot", "subplot", "apr2018SL")])
  Sl_per_subplot_nov18 <- unique(data3[, c("plot", "subplot", "nov18SL")])
  Sl_per_subplot_apr19 <- unique(data3[, c("plot", "subplot", "apr19SL")])
  
  # calculate the sum of seedling per plot
  Sl_per_plot_apr18 <- aggregate(Sl_per_subplot_apr18$apr2018SL, 
                                 by = list(Sl_per_subplot_apr18$plot), FUN = sum, na.rm = T)
  Sl_per_plot_nov18 <- aggregate(Sl_per_subplot_nov18$nov18SL, 
                                 by = list(Sl_per_subplot_nov18$plot), FUN = sum, na.rm = T)
  Sl_per_plot_apr19 <- aggregate(Sl_per_subplot_apr19$apr19SL, 
                                 by = list(Sl_per_subplot_apr19$plot), FUN = sum, na.rm = T)
  
  # calculate the sum of seeds per plot
  seedling_per_seed <- aggregate(data3$seed_per_ind, by = list(data3$plot), FUN = sum, na.rm = T)
  
  # calculate the turnover of seeds to seedlings
  seedling_per_seed$seed_SL_nov18 <- Sl_per_plot_nov18$x / seedling_per_seed$x
  seedling_per_seed$seed_SL_apr19 <- Sl_per_plot_apr19$x / seedling_per_seed$x
  seedling_per_seed[seedling_per_seed == Inf] <- NA
  
  #--------------------------------
  # Seedling survival
  #--------------------------------
  sl_surv <- sum(data3$new_plant, na.rm = T) / (sum(Sl_per_plot_apr18$x) + sum(Sl_per_plot_nov18$x))
  
  #--------------------------------
  # new Individual size
  #--------------------------------
  data_new_ind <- subset(data, new_plant == 1)
  new_size_mean <- mean(data_new_ind$logsizet1, na.rm = T)
  new_size_sd <- sd(data_new_ind$logsizet1, na.rm = T)
  
  # store the parameters of the models
  pars <- list(  surv_b0       = coef(sr_mod)[1],
                 surv_b1       = coef(sr_mod)[2],
                 grow_b0       = coef(gr_mod)[1],
                 grow_b1       = coef(gr_mod)[2],
                 grow_sd       = summary(gr_mod)$sigma,
                 seednum_b0    = coef(sefl_mod)[1],
                 seednum_b1    = coef(sefl_mod)[2],
                 flowprob_b0   = coef(fl_mod)[1],
                 flowprob_b1   = coef(fl_mod)[2],
                 new_size_mean  = new_size_mean,
                 new_size_sd    = new_size_sd,
                 novSL_per_fl   = mean(seedling_per_seed$seed_SL_nov18, na.rm = T),
                 aprilSL_per_fl = mean(seedling_per_seed$seed_SL_apr19, na.rm = T),
                 SL_surv        = sl_surv,
                 L = min( c(data$logsizet0,
                            data$logsizet1),
                          na.rm=T),
                 U = max( c(data$logsizet0,
                            data$logsizet1),
                          na.rm=T),
                 L_obs = min( c(data$logsizet0,
                                data$logsizet1),
                              na.rm=T),
                 U_obs = max( c(data$logsizet0,
                                data$logsizet1),
                              na.rm=T),
                 n = 200
  )
  
  # calculate Lambda
  Re(eigen(ker_mee(pars)$k_yx)$value[1])
}

####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####           Read data
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("C://Users/ma22buky/Documents/Julia_Paper/Done/")
dat_d <- read.csv("Bro_Demography.csv")
dat_s <- read.csv("Bro_Seed_data.csv")

## Delete the new plants which grow bigger then 30
dat_d <- dat_d[-(which(dat_d$sizet1 >= 30 & dat_d$new_plant == 1)),]

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##            Ambient mowing
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
demo_dat <- subset(dat_d, treatment == "ambient_mowing")
seed_dat <- subset(dat_s, treatment == "ambient_mowing")

lamb_amb_mow <- sapply(1:10, FUN = boot_lam)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##            Ambient grazing
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
demo_dat <- subset(dat_d, treatment == "ambient_grazing")
seed_dat <- subset(dat_s, treatment == "ambient_grazing")

lamb_amb_gra <- sapply(1:10, FUN = boot_lam)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##          Future mowing
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
demo_dat <- subset(dat_d, treatment == "future_mowing")
seed_dat <- subset(dat_s, treatment == "future_mowing")

lamb_fut_mow <- sapply(1:10, FUN = boot_lam)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##          Future grazing
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
demo_dat <- subset(dat_d, treatment == "future_grazing")
seed_dat <- subset(dat_s, treatment == "future_grazing")

lamb_fut_gra <- sapply(1:10, FUN = boot_lam)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##          Plot mean Lambdas
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Create data frame with mean and Confidence intervalls
Results <- data.frame(lambda_mean = c(mean(lamb_amb_mow), mean(lamb_amb_gra), 
                                      mean(lamb_fut_mow), mean(lamb_fut_gra)),
                      lambda_U = c(quantile(lamb_amb_mow, probs = c(0.025, 0.975))[2],
                                   quantile(lamb_amb_gra, probs = c(0.025, 0.975))[2],
                                   quantile(lamb_fut_mow, probs = c(0.025, 0.975))[2],
                                   quantile(lamb_fut_gra, probs = c(0.025, 0.975))[2]),
                      lambda_L = c(quantile(lamb_amb_mow, probs = c(0.025, 0.975))[1],
                                   quantile(lamb_amb_gra, probs = c(0.025, 0.975))[1],
                                   quantile(lamb_fut_mow, probs = c(0.025, 0.975))[1],
                                   quantile(lamb_fut_gra, probs = c(0.025, 0.975))[1]),
                      treatment = c("ambient_mowing", "ambient_grazing", "future_mowing", "future_grazing"),
                      landuse = c("mowing", "grazing", "mowing", "grazing"),
                      climate = c("ambient", "ambient", "future", "future"))

## Choose color palette
myPalette <- c("#0072B2", "#D55E00")

ggplot(Results, aes(x = climate, y = lambda_mean, color = climate, shape = landuse)) + theme_classic() +
  geom_line(aes(group = landuse ), col = "black") + scale_shape_manual(values=c(16, 17)) +
  scale_colour_manual(values=myPalette) +
  geom_point( size = 5, position = position_dodge(width = 0.08)) +
  geom_errorbar(aes(ymax = lambda_U, ymin = lambda_L), width = .17, size = 1.2, position = position_dodge(0.08)) +
  theme(text = element_text(size = 20), 
        legend.position = c(0.87, 0.89), 
        legend.background = element_rect(fill="transparent"),
        legend.title = element_blank())  +
  ylab("Asymptotic population growth rate (λ)") + xlab("Climate treatment") + ggtitle("Population growth rates") + 
  theme(legend.position = "none")