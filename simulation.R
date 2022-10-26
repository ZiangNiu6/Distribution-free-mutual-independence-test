
########Type-I error comparison################
#################################################



########Gaussian Setting################

# compute threshold
n <- 500
set.seed(1)
emprical <- gensamdistrjdcov(n, dim_list = rep(3, 3), niter=1000)

# specify the group of variables
group <- c(1,1,1,2,2,2,3,3,3)
type1_rdhdcov_G <- rep(0, 500)
type1_hdcov_G <- rep(0, 500)
type1_matteson <- rep(0, 500)
type1_dhsic <- rep(0, 500)

# 1000 replications
for (k in 1:500) {
  X <- matrix(0, nrow=n, ncol = 9)
  a <- c(1, 4, 7)
  b <- c(2,3,5,6,8,9)
  #  ber <- rbinom(rep(1, 1500), rep(1, 1500), prob = rep(0.05, 1500))
  #  X[,1] <- ber[1:500]*rmvc(500, mu=rep(0,1), S=1) + (1-ber[1:500])*rmvnorm(500, mean = rep(0,1), sigma = diag(1))
  #  X[,4] <- ber[501:1000]*rmvc(500, mu=rep(0,1), S=1) + (1-ber[501:1000])*rmvnorm(500, mean = rep(0,1), sigma = diag(1))
  #  X[,7] <- ber[1001:1500]*rmvc(500, mu=rep(0,1), S=1) + (1-ber[1001:1500])*rmvnorm(500, mean = rep(0,1), sigma = diag(1))
  X[,a] <- rmvnorm(n, mean = rep(0,3), sigma = diag(3))
  X[,b] <- rmvnorm(n, mean = rep(0,6), sigma = diag(6))
  X_list <- list(X[,1:3], X[, 4:6], X[, 7:9])
  statis <- computestatisticjdcov(X, dim_list = rep(3, 3))
  type1_rdhdcov_G[k] <- length(which(emprical >= statis)) / (length(emprical) + 1)
  print(type1_rdhdcov_G[k])
  
  # JdCov
  statis <- jdcov.test(X_list, stat.type = "V", B = 500, alpha = 0.05)
  type1_hdcov_G[k] <- statis$p.value
  print(type1_hdcov_G[k])
  
  
  ## matteson's test
  type1_matteson[k] <- boot_matteson(X, B =500, group = group)
  print(type1_matteson[k])
  
  ## dHSIC test
  type1_dhsic[k] <- dhsic.test(X_list, B = 500)$p.value
  print(type1_dhsic[k])
  
  print(k)
}

Gaussian_type1 <- data.frame(rdhdcov = type1_rdhdcov_G, hdcov = type1_hdcov_G, matteson = type1_matteson, dHSIC = type1_dhsic)

write.csv(Gaussian_type1, "C:/Users/ziangniu/Desktop/Simulation/result/Gaussian_type1.csv")

################Copula Setting######################
n <- 500
set.seed(1)
group <- c(1,1,1,2,2,2,3,3,3)
emprical <- gensamdistrjdcov(n, dim_list = rep(3, 3), niter=1000)
type1_rdhdcov_CG <- rep(0, 500)
type1_hdcov_CG <- rep(0, 500)
type1_matteson_CG <- rep(0, 500)
type1_dhsic_CG <- rep(0, 500)
for (k in 1:500) {
  X <- matrix(0, nrow=500, ncol = 9)
  a <- c(1, 4, 7)
  b <- c(2,3,5,6,8,9)
  #  ber <- rbinom(rep(1, 1500), rep(1, 1500), prob = rep(0.05, 1500))
  #  X[,1] <- ber[1:500]*rmvc(500, mu=rep(0,1), S=1) + (1-ber[1:500])*rmvnorm(500, mean = rep(0,1), sigma = diag(1))
  #  X[,4] <- ber[501:1000]*rmvc(500, mu=rep(0,1), S=1) + (1-ber[501:1000])*rmvnorm(500, mean = rep(0,1), sigma = diag(1))
  #  X[,7] <- ber[1001:1500]*rmvc(500, mu=rep(0,1), S=1) + (1-ber[1001:1500])*rmvnorm(500, mean = rep(0,1), sigma = diag(1))
  X[,a] <- (rmvnorm(500, mean = rep(0,3), sigma = diag(3)))**3
  X[,b] <- (rmvnorm(500, mean = rep(0,6), sigma = diag(6)))**3
  X_list <- list(X[,1:3], X[, 4:6], X[, 7:9])
  statis <- computestatisticjdcov(X, dim_list = rep(3, 3))
  type1_rdhdcov_CG[k] <- length(which(emprical >= statis)) / (length(emprical) + 1)
  print(type1_rdhdcov_CG[k])
  
  # JdCov
  statis <- jdcov.test(X_list, stat.type = "V", B = 500, alpha = 0.05)
  type1_hdcov_CG[k] <- statis$p.value
  print(type1_hdcov_CG[k])
  
  
  ## matteson's test
  type1_matteson_CG[k] <- boot_matteson(X, B =500, group = group)
  print(type1_matteson_CG[k])
  
  ## dHSIC test
  type1_dhsic_CG[k] <- dhsic.test(X_list, B = 500)$p.value
  print(type1_dhsic_CG[k])
  
  print(k)
}


CopulaGaussian_type1 <- data.frame(rdhdcov = type1_rdhdcov_CG, hdcov = type1_hdcov_CG, matteson = type1_matteson_CG, dHSIC = type1_dhsic_CG)

write.csv(CopulaGaussian_type1, "C:/Users/ziangniu/Desktop/Simulation/result/CopulaGaussian_type1.csv")




########Cauchy Setting################
n <- 500
set.seed(1)
group <- c(1,1,1,2,2,2,3,3,3)
emprical <- gensamdistrjdcov(n, dim_list = rep(3, 3), niter=1000)
type1_rdhdcov_C <- rep(0, 500)
type1_hdcov_C <- rep(0, 500)
type1_matteson_C <- rep(0, 500)
type1_dhsic_C <- rep(0, 500)
for (k in 1:500) {
  X <- matrix(rcauchy(9*500, 0, 1), nrow=n, ncol=9)
  #  X[,b] <- rmvnorm(500, mean = rep(0,6), sigma = diag(6))
  X_list <- list(X[,1:3], X[, 4:6], X[, 7:9])
  
  # RJdCov
  statis <- computestatisticjdcov(X, dim_list = rep(3, 3))
  type1_rdhdcov_C[k] <- length(which(emprical >= statis)) / (length(emprical) + 1)
  print(type1_rdhdcov_C[k])
  
  # JdCov
  statis <- jdcov.test(X_list, stat.type = "V", B = 500, alpha = 0.05)
  type1_hdcov_C[k] <- statis$p.value
  print(type1_hdcov_C[k])
  
  
  ## matteson's test
  type1_matteson_C[k] <- boot_matteson(X, B =500, group = group)
  print(type1_matteson_C[k])
  
  ## dHSIC test
  type1_dhsic_C[k] <- dhsic.test(X_list, B = 500)$p.value
  print(type1_dhsic_C[k])
  
  print(k)
}

Cauchy_type1 <- data.frame(rdhdcov = type1_rdhdcov_C, hdcov = type1_hdcov_C, matteson = type1_matteson_C, dHSIC = type1_dhsic_C)

write.csv(Cauchy_type1, "C:/Users/ziangniu/Desktop/Simulation/result/Cauchy_type1.csv")




########Power Analysis################
#######################################


###matrix###
sig_T <- function(rho, d){
  sig_T <- matrix(0, nrow=d, ncol=d)
  for (i in 1:d) {
    for (j in 1:d) {
      sig_T[i,j] <- rho**(abs(i-j))
    }
  }
  sig_T
}


###Gaussian### vary rho within \{0.005,0.01,0.015,0.02,0.025\}

n <- 500
set.seed(1)
emprical <- gensamdistrjdcov(n, dim_list = rep(3, 3), niter=1000)
group <- c(1,1,1,2,2,2,3,3,3)
power_rdhdcov_GT2 <- matrix(0, 500, 3)
power_hdcov_GT2 <- matrix(0, 500, 3)
power_matteson_GT2 <- matrix(0, 500, 3)
power_dhsic_GT2 <- matrix(0, 500, 3)
rho <- seq(0.15,0.25,0.05)
for (i in 1:length(rho)) {
  for (k in 1:500) {
    X <- as.matrix(rmvnorm(n, mean = rep(0,9), sigma = sig_T(rho[i], 9)))
    X_list <- list(X[,1:3], X[, 4:6], X[, 7:9])
    # RJdCov
    statis <- computestatisticjdcov(X, dim_list = rep(3, 3))
    power_rdhdcov_GT2[k, i] <- length(which(emprical >= statis)) / (length(emprical) + 1)
    print(power_rdhdcov_GT2[k, i])
    
    # JdCov
    statis <- jdcov.test(X_list, stat.type = "V", B = 500, alpha = 0.05)
    power_hdcov_GT2[k, i] <- statis$p.value
    print(power_hdcov_GT2[k, i])
    
    
    ## matteson's test
    power_matteson_GT2[k, i] <- boot_matteson(X, B =500, group = group)
    print(power_matteson_GT2[k, i])
    
    ## dHSIC test
    power_dhsic_GT2[k, i] <- dhsic.test(X_list, B = 500)$p.value
    print(power_dhsic_GT2[k, i])
    
    print(k)
  }
  print("finish one round")
}

GT2_power <- data.frame(rdhdcov = power_rdhdcov_GT2, hdcov = power_hdcov_GT2, matteson = power_matteson_GT2, dHSIC = power_dhsic_GT2)

write.csv(GT2_power, "C:/Users/ziangniu/Desktop/Simulation/result/GT2_power.csv")






#######banded covariance matrix setting########

sig_B <- function(rho, d){
  sig_B <- matrix(0, nrow=d, ncol=d)
  for (i in 1:d) {
    for (j in 1:d) {
      if(abs(i-j)<=2){
        if( i != j){
          sig_B[i,j] <- rho
        }else{
          sig_B[i, j] <- 1
        }
      }
    }
  }
  sig_B
}



n <- 500
set.seed(1)
emprical <- gensamdistrjdcov(n, dim_list = rep(3, 3), niter=1000)
group <- c(1,1,1,2,2,2,3,3,3)
power_rdhdcov_GB2 <- matrix(0, 500, 3)
power_hdcov_GB2 <- matrix(0, 500, 3)
power_matteson_GB2 <- matrix(0, 500, 3)
power_dhsic_GB2 <- matrix(0, 500, 3)
rho <- seq(0.15,0.25,0.05)
for (i in 1:length(rho)) {
  for (k in 1:500) {
    X <- as.matrix(rmvnorm(n, mean = rep(0,9), sigma = sig_B(rho[i], 9)))
    X_list <- list(X[,1:3], X[, 4:6], X[, 7:9])
    # RJdCov
    statis <- computestatisticjdcov(X, dim_list = rep(3, 3))
    power_rdhdcov_GB2[k, i] <- length(which(emprical >= statis)) / (length(emprical) + 1)
    print(power_rdhdcov_GB2[k, i])
    
    # JdCov
    statis <- jdcov.test(X_list, stat.type = "V", B = 500, alpha = 0.05)
    power_hdcov_GB2[k, i] <- statis$p.value
    print(power_hdcov_GB2[k, i])
    
    
    ## matteson's test
    power_matteson_GB2[k, i] <- boot_matteson(X, B =500, group = group)
    print(power_matteson_GB2[k, i])
    
    ## dHSIC test
    power_dhsic_GB2[k, i] <- dhsic.test(X_list, B = 500)$p.value
    print(power_dhsic_GB2[k, i])
    
    print(k)
  }
  print("finish one round")
}

GB2_power <- data.frame(rdhdcov = power_rdhdcov_GB2, hdcov = power_hdcov_GB2, matteson = power_matteson_GB2, dHSIC = power_dhsic_GB2)

write.csv(GB2_power, "C:/Users/ziangniu/Desktop/Simulation/result/GB2_power.csv")





#######heavy-tailed regression#########


n <- 500
set.seed(1)
emprical <- gensamdistrjdcov(n, dim_list = rep(3, 3), niter=1000)
group <- c(1,1,1,2,2,2,3,3,3)
power_rdhdcov_H2 <- matrix(0, 500, 3)
power_hdcov_H2 <- matrix(0, 500, 3)
power_matteson_H2 <- matrix(0, 500, 3)
power_dhsic_H2 <- matrix(0, 500, 3)
rho <- seq(0.06, 0.1, length.out = 3)
for (i in 1:length(rho)) {
  for (k in 1:500) {
    X <- matrix(rcauchy(500*9, rep(0, 500*9), rep(1, 500*9)), nrow = 500, ncol = 9) + rho[i]*rcauchy(500)
    X_list <- list(X[,1:3], X[, 4:6], X[, 7:9])
    
    # RJdCov
    statis <- computestatisticjdcov(X, dim_list = rep(3, 3))
    power_rdhdcov_H2[k, i] <- length(which(emprical >= statis)) / (length(emprical) + 1)
    print(power_rdhdcov_H2[k, i])
    
    # JdCov
    statis <- jdcov.test(X_list, stat.type = "V", B = 500, alpha = 0.05)
    power_hdcov_H2[k, i] <- statis$p.value
    print(power_hdcov_H2[k, i])
    
    
    ## matteson's test
    power_matteson_H2[k, i] <- boot_matteson(X, B =500, group = group)
    print(power_matteson_H2[k, i])
    
    ## dHSIC test
    power_dhsic_H2[k, i] <- dhsic.test(X_list, B = 500)$p.value
    print(power_dhsic_H2[k, i])
    
    print(k)
  }
  print("finish one round")
}

H2_power <- data.frame(rdhdcov = power_rdhdcov_H2, hdcov = power_hdcov_H2, matteson = power_matteson_H2, dHSIC = power_dhsic_H2)

write.csv(H2_power, "C:/Users/ziangniu/Desktop/Simulation/result/H2_power.csv")



#########sinusoid distribution###############


# rejection sampling

rej_sam <- function(l, n, d, dim){
  data <- matrix(0, nrow = n, ncol = d*dim)
  for (i in 1:dim) {
    from_ind <- 1
    while (from_ind < n) {
      X <- matrix(runif(n*d, -pi, pi), nrow = n, ncol = d)
      un_den <- 1 + apply(sin(l*X), 1, prod)
      u <- runif(n)
      accept_id <- which(u < un_den / 2)
      # specify the end index
      end_ind <- min(n, from_ind + length(accept_id) - 1)
      no_sam <- end_ind - from_ind + 1
      data[from_ind:end_ind, seq(i, d*dim, by = d)] <- X[accept_id[1:no_sam], ]
      from_ind <- end_ind + 1
    }
  }
  data
}

###### vary l

n <- 500
emprical <- gensamdistrjdcov(n, dim_list = rep(3, 3), niter=1000, fixgrid = matrix(rnorm(n*9), n, 9), c = 10)
set.seed(1)
group <- c(1,1,1,2,2,2,3,3,3)
power_rdhdcov_sin <- matrix(0, 500, 5)
power_hdcov_sin <- matrix(0, 500, 5)
power_matteson_sin <- matrix(0, 500, 5)
power_dhsic_sin <- matrix(0, 500, 5)
l <- c(0.1, 0.5, 1, 2, 5)
for (i in 1:length(l)) {
  for (k in 1:500) {
    X <- rej_sam(l[i], n, 3, 3)
    X_list <- list(X[,1:3], X[, 4:6], X[, 7:9])
    
    # RJdCov
    statis <- computestatisticjdcov(X, dim_list = rep(3, 3), gridch = matrix(rnorm(n*9), n, 9),  c = 10)
    power_rdhdcov_sin[k, i] <- length(which(emprical >= statis)) / (length(emprical) + 1)
    print(power_rdhdcov_sin[k, i])
    
    # JdCov
    statis <- jdcov.test(X_list, stat.type = "V", B = 500, alpha = 0.05)
    power_hdcov_sin[k, i] <- statis$p.value
    print(power_hdcov_sin[k, i])
    
    
    ## matteson's test
    power_matteson_sin[k, i] <- boot_matteson(X, B =500, group = group)
    print(power_matteson_sin[k, i])
    
    ## dHSIC test
    power_dhsic_sin[k, i] <- dhsic.test(X_list, B = 500)$p.value
    print(power_dhsic_sin[k, i])
    
    print(k)
  }
  print("finish one round")
}

H_power <- data.frame(rdhdcov = power_rdhdcov_H, hdcov = power_hdcov_H, matteson = power_matteson_H, dHSIC = power_dhsic_H)

write.csv(H_power, "C:/Users/ziangniu/Desktop/Simulation/result/H_power.csv")


### additive sin dependence#####


n <- 500
set.seed(1)
emprical <- gensamdistrjdcov(n, dim_list = rep(3, 3), niter=1000)
group <- c(1,1,1,2,2,2,3,3,3)
power_rdhdcov_s2 <- matrix(0, 500, 3)
power_hdcov_s2 <- matrix(0, 500, 3)
power_matteson_s2<- matrix(0, 500, 3)
power_dhsic_s2 <- matrix(0, 500, 3)
l <- seq(0.3, 0.5, length.out = 3)
for (i in 1:length(l)) {
  for (k in 1:500) {
    X <- matrix(rcauchy(500*9, rep(0, 500*9), rep(1, 500*9)), nrow = 500, ncol = 9) + sin(l[i]*rcauchy(500))
    X_list <- list(X[,1:3], X[, 4:6], X[, 7:9])
    
    # RJdCov
    statis <- computestatisticjdcov(X, dim_list = rep(3, 3))
    power_rdhdcov_s2[k, i] <- length(which(emprical >= statis)) / (length(emprical) + 1)
    print(power_rdhdcov_s2[k, i])
    
    # JdCov
    statis <- jdcov.test(X_list, stat.type = "V", B = 500, alpha = 0.05)
    power_hdcov_s2[k, i] <- statis$p.value
    print(power_hdcov_s2[k, i])
    
    
    ## matteson's test
    power_matteson_s2[k, i] <- boot_matteson(X, B =500, group = group)
    print(power_matteson_s2[k, i])
    
    ## dHSIC test
    power_dhsic_s2[k, i] <- dhsic.test(X_list, B = 500)$p.value
    print(power_dhsic_s2[k, i])
    
    print(k)
  }
  print("finish one round")
}

s2_power <- data.frame(rdhdcov = power_rdhdcov_s2, hdcov = power_hdcov_s2, matteson = power_matteson_s2, dHSIC = power_dhsic_s2)

write.csv(s2_power, "C:/Users/ziangniu/Desktop/Simulation/result/s2_power.csv")


########higher-order dependence########

###generate data######


## test the pairwise independence
n <- 500
set.seed(1)
emprical <- gensamdistrhodcov(500, dim_list = rep(3, 2), niter=1000)
power_rdhdcov_Ghod <- rep(0, 500)
power_hdcov_Ghod <- rep(0, 500)
for (k in 1:500) {
  X <- matrix(0, nrow=500, ncol = 9)
  a <- c(1, 4, 7)
  b <- c(2,3,5,6,8,9)
  X <- matrix(rcauchy(500*9, 0, 1), nrow = 500, ncol = 9)
  sign_list <- X[, 1]*X[, 4]*X[, 7]
  plus_index <- which(sign_list<=0)
  neg_index <- which(sign_list>0)
  X[, 7][plus_index] <- -X[, 7][plus_index]
  X[, 7][neg_index] <- X[, 7][neg_index]
  X_list <- list(X[,4:6], X[, 7:9])
  
  ## rdhod
  statis <- computestatisticjdcov(X[,4:9], dim_list = rep(3, 2))
  power_rdhdcov_Ghod[k] <- length(which(emprical >= statis)) / (length(emprical) + 1)
  
  ## jdhod
  statis <- jdcov.test(X_list, stat.type = "V", alpha = 0.05)
  power_hdcov_Ghod[k] <- statis$p.value
  
  
  print(power_rdhdcov_Ghod[k])
  print(power_hdcov_Ghod[k])
  print(k)
}

mix_2_power <- data.frame(rdhdcov = power_rdhdcov_Ghod, hdcov = power_hdcov_Ghod)

write.csv(mix_2_power, "C:/Users/ziangniu/Desktop/Simulation/result/mix_2_power.csv")


## test the third order dependence

n <- 500
set.seed(1)
emprical <- gensamdistrhodcov(500, dim_list = rep(3, 3), niter=1000)
power_rdhdcov_Ghod_3 <- rep(0, 500)
power_hdcov_Ghod_3 <- rep(0, 500)
for (k in 1:500) {
  X <- matrix(0, nrow=500, ncol = 9)
  a <- c(1, 4, 7)
  b <- c(2,3,5,6,8,9)
  X <- matrix(rcauchy(500*9, 0, 1), nrow = 500, ncol = 9)
  sign_list <- X[, 1]*X[, 4]*X[, 7]
  plus_index <- which(sign_list<=0)
  neg_index <- which(sign_list>0)
  X[, 7][plus_index] <- -X[, 7][plus_index]
  X[, 7][neg_index] <- X[, 7][neg_index]
  X_list <- list(X[,1:3], X[,4:6], X[,7:9])
  
  ## rdhod
  statis <- computestatisticrdcov(X, dim_list = rep(3, 3))
  power_rdhdcov_Ghod_3[k] <- length(which(emprical >= statis)) / (length(emprical) + 1)
  
  ## jdhod
  statis <- boot_hodcov(X_list, 500)
  power_hdcov_Ghod_3[k] <- statis
  
  
  print(power_rdhdcov_Ghod_3[k])
  print(power_hdcov_Ghod_3[k])
  print(k)
}

mix_3_power <- data.frame(rdhdcov = power_rdhdcov_Ghod_3, hdcov = power_hdcov_Ghod_3)

write.csv(mix_3_power, "C:/Users/ziangniu/Desktop/Simulation/result/mix_3_power.csv")


## test the joint independence
n <- 500
set.seed(1)
emprical <- gensamdistrjdcov(500, dim_list = rep(3, 3), niter=1000)
power_rdhdcov_P <- rep(0, 500)
power_hdcov_P <- rep(0, 500)
for (k in 1:500) {
  X <- matrix(0, nrow=500, ncol = 9)
  a <- c(1, 4, 7)
  b <- c(2,3,5,6,8,9)
  X <- matrix(rcauchy(500*9, 0, 1), nrow = 500, ncol = 9)
  sign_list <- X[, 1]*X[, 4]*X[, 7]
  plus_index <- which(sign_list<=0)
  neg_index <- which(sign_list>0)
  X[, 7][plus_index] <- -X[, 7][plus_index]
  X[, 7][neg_index] <- X[, 7][neg_index]
  X_list <- list(X[,1:3], X[, 4:6], X[, 7:9])
  
  ## rdhod
  statis <- computestatisticjdcov(X, dim_list = rep(3, 3))
  power_rdhdcov_P[k] <- length(which(emprical >= statis)) / (length(emprical) + 1)

  statis <- jdcov.test(X_list, stat.type = "V", alpha = 0.05)
  power_hdcov_P[k] <- statis$p.value
  
  
  print(power_rdhdcov_P[k])
  print(power_hdcov_P[k])
  print(k)
}

mix_P_power <- data.frame(rdhdcov = power_rdhdcov_P, hdcov = power_hdcov_P)

write.csv(mix_P_power, "C:/Users/ziangniu/Desktop/Simulation/result/mix_P_power.csv")


###Gaussian case######

## test the pairwise independence
n <- 500
set.seed(1)
emprical <- gensamdistrhodcov(500, dim_list = rep(3, 2), niter=1000)
power_rdhdcov_Ghod <- rep(0, 500)
power_hdcov_Ghod <- rep(0, 500)
for (k in 1:500) {
  X <- matrix(0, nrow=500, ncol = 9)
  a <- c(1, 4, 7)
  b <- c(2,3,5,6,8,9)
  X <- matrix(rnorm(500*9, 0, 1), nrow = 500, ncol = 9)
  sign_list <- X[, 1]*X[, 4]*X[, 7]
  plus_index <- which(sign_list<=0)
  neg_index <- which(sign_list>0)
  X[, 7][plus_index] <- -X[, 7][plus_index]
  X[, 7][neg_index] <- X[, 7][neg_index]
  X_list <- list(X[,4:6], X[, 7:9])
  
  ## rdhod
  statis <- computestatisticjdcov(X[,4:9], dim_list = rep(3, 2))
  power_rdhdcov_Ghod[k] <- length(which(emprical >= statis)) / (length(emprical) + 1)
  
  ## jdhod
  statis <- jdcov.test(X_list, stat.type = "V", alpha = 0.05)
  power_hdcov_Ghod[k] <- statis$p.value
  
  
  print(power_rdhdcov_Ghod[k])
  print(power_hdcov_Ghod[k])
  print(k)
}

mix_2_power <- data.frame(rdhdcov = power_rdhdcov_Ghod, hdcov = power_hdcov_Ghod)

write.csv(mix_2_power, "C:/Users/ziangniu/Desktop/Simulation/result/mix_G2_power.csv")


## test the third order dependence

n <- 500
set.seed(1)
emprical <- gensamdistrhodcov(500, dim_list = rep(3, 3), niter=1000)
power_rdhdcov_Ghod_3 <- rep(0, 500)
power_hdcov_Ghod_3 <- rep(0, 500)
for (k in 1:500) {
  X <- matrix(0, nrow=500, ncol = 9)
  a <- c(1, 4, 7)
  b <- c(2,3,5,6,8,9)
  X <- matrix(rnorm(500*9, 0, 1), nrow = 500, ncol = 9)
  sign_list <- X[, 1]*X[, 4]*X[, 7]
  plus_index <- which(sign_list<=0)
  neg_index <- which(sign_list>0)
  X[, 7][plus_index] <- -X[, 7][plus_index]
  X[, 7][neg_index] <- X[, 7][neg_index]
  X_list <- list(X[,1:3], X[,4:6], X[,7:9])
  
  ## rdhod
  statis <- computestatisticrdcov(X, dim_list = rep(3, 3))
  power_rdhdcov_Ghod_3[k] <- length(which(emprical >= statis)) / (length(emprical) + 1)
  
  ## jdhod
  statis <- boot_hodcov(X_list, 500)
  power_hdcov_Ghod_3[k] <- statis
  
  
  print(power_rdhdcov_Ghod_3[k])
  print(power_hdcov_Ghod_3[k])
  print(k)
}

mix_3_power <- data.frame(rdhdcov = power_rdhdcov_Ghod_3, hdcov = power_hdcov_Ghod_3)

write.csv(mix_3_power, "C:/Users/ziangniu/Desktop/Simulation/result/mix_G3_power.csv")


## test the joint independence
n <- 500
set.seed(1)
emprical <- gensamdistrjdcov(500, dim_list = rep(3, 3), niter=1000)
power_rdhdcov_P <- rep(0, 500)
power_hdcov_P <- rep(0, 500)
for (k in 1:500) {
  X <- matrix(0, nrow=500, ncol = 9)
  a <- c(1, 4, 7)
  b <- c(2,3,5,6,8,9)
  X <- matrix(rnorm(500*9, 0, 1), nrow = 500, ncol = 9)
  sign_list <- X[, 1]*X[, 4]*X[, 7]
  plus_index <- which(sign_list<=0)
  neg_index <- which(sign_list>0)
  X[, 7][plus_index] <- -X[, 7][plus_index]
  X[, 7][neg_index] <- X[, 7][neg_index]
  X_list <- list(X[,1:3], X[, 4:6], X[, 7:9])
  
  ## rdhod
  statis <- computestatisticjdcov(X, dim_list = rep(3, 3))
  power_rdhdcov_P[k] <- length(which(emprical >= statis)) / (length(emprical) + 1)
  
  statis <- jdcov.test(X_list, stat.type = "V", alpha = 0.05)
  power_hdcov_P[k] <- statis$p.value
  
  
  print(power_rdhdcov_P[k])
  print(power_hdcov_P[k])
  print(k)
}

mix_P_power <- data.frame(rdhdcov = power_rdhdcov_P, hdcov = power_hdcov_P)

write.csv(mix_P_power, "C:/Users/ziangniu/Desktop/Simulation/result/mix_GP_power.csv")


###student case with degree 3######

## test the pairwise independence
n <- 500
set.seed(1)
emprical <- gensamdistrhodcov(500, dim_list = rep(3, 2), niter=1000)
power_rdhdcov_Ghod <- rep(0, 500)
power_hdcov_Ghod <- rep(0, 500)
for (k in 1:500) {
  X <- matrix(0, nrow=500, ncol = 9)
  a <- c(1, 4, 7)
  b <- c(2,3,5,6,8,9)
  X <- matrix(rt(500*9, 3), nrow = 500, ncol = 9)
  sign_list <- X[, 1]*X[, 4]*X[, 7]
  plus_index <- which(sign_list<=0)
  neg_index <- which(sign_list>0)
  X[, 7][plus_index] <- -X[, 7][plus_index]
  X[, 7][neg_index] <- X[, 7][neg_index]
  X_list <- list(X[,4:6], X[, 7:9])
  
  ## rdhod
  statis <- computestatisticjdcov(X[,4:9], dim_list = rep(3, 2))
  power_rdhdcov_Ghod[k] <- length(which(emprical >= statis)) / (length(emprical) + 1)
  
  ## jdhod
  statis <- jdcov.test(X_list, stat.type = "V", alpha = 0.05)
  power_hdcov_Ghod[k] <- statis$p.value
  
  
  print(power_rdhdcov_Ghod[k])
  print(power_hdcov_Ghod[k])
  print(k)
}

mix_2_power <- data.frame(rdhdcov = power_rdhdcov_Ghod, hdcov = power_hdcov_Ghod)

write.csv(mix_2_power, "C:/Users/ziangniu/Desktop/Simulation/result/mix_T2_power.csv")


## test the third order dependence

n <- 500
set.seed(1)
emprical <- gensamdistrhodcov(500, dim_list = rep(3, 3), niter=1000)
power_rdhdcov_Ghod_3 <- rep(0, 500)
power_hdcov_Ghod_3 <- rep(0, 500)
for (k in 1:500) {
  X <- matrix(0, nrow=500, ncol = 9)
  a <- c(1, 4, 7)
  b <- c(2,3,5,6,8,9)
  X <- matrix(rt(500*9, 3), nrow = 500, ncol = 9)
  sign_list <- X[, 1]*X[, 4]*X[, 7]
  plus_index <- which(sign_list<=0)
  neg_index <- which(sign_list>0)
  X[, 7][plus_index] <- -X[, 7][plus_index]
  X[, 7][neg_index] <- X[, 7][neg_index]
  X_list <- list(X[,1:3], X[,4:6], X[,7:9])
  
  ## rdhod
  statis <- computestatisticrdcov(X, dim_list = rep(3, 3))
  power_rdhdcov_Ghod_3[k] <- length(which(emprical >= statis)) / (length(emprical) + 1)
  
  ## jdhod
  statis <- boot_hodcov(X_list, 500)
  power_hdcov_Ghod_3[k] <- statis
  
  
  print(power_rdhdcov_Ghod_3[k])
  print(power_hdcov_Ghod_3[k])
  print(k)
}

mix_3_power <- data.frame(rdhdcov = power_rdhdcov_Ghod_3, hdcov = power_hdcov_Ghod_3)

write.csv(mix_3_power, "C:/Users/ziangniu/Desktop/Simulation/result/mix_T3_power.csv")


## test the joint independence
n <- 500
set.seed(1)
emprical <- gensamdistrjdcov(500, dim_list = rep(3, 3), niter=1000)
power_rdhdcov_P <- rep(0, 500)
power_hdcov_P <- rep(0, 500)
for (k in 1:500) {
  X <- matrix(0, nrow=500, ncol = 9)
  a <- c(1, 4, 7)
  b <- c(2,3,5,6,8,9)
  X <- matrix(rt(500*9, 3), nrow = 500, ncol = 9)
  sign_list <- X[, 1]*X[, 4]*X[, 7]
  plus_index <- which(sign_list<=0)
  neg_index <- which(sign_list>0)
  X[, 7][plus_index] <- -X[, 7][plus_index]
  X[, 7][neg_index] <- X[, 7][neg_index]
  X_list <- list(X[,1:3], X[, 4:6], X[, 7:9])
  
  ## rdhod
  statis <- computestatisticjdcov(X, dim_list = rep(3, 3))
  power_rdhdcov_P[k] <- length(which(emprical >= statis)) / (length(emprical) + 1)
  
  statis <- jdcov.test(X_list, stat.type = "V", alpha = 0.05)
  power_hdcov_P[k] <- statis$p.value
  
  
  print(power_rdhdcov_P[k])
  print(power_hdcov_P[k])
  print(k)
}

mix_P_power <- data.frame(rdhdcov = power_rdhdcov_P, hdcov = power_hdcov_P)

write.csv(mix_P_power, "C:/Users/ziangniu/Desktop/Simulation/result/mix_TP_power.csv")



###student case with degree 2######

## test the pairwise independence
n <- 500
set.seed(1)
emprical <- gensamdistrhodcov(500, dim_list = rep(3, 2), niter=1000)
power_rdhdcov_Ghod <- rep(0, 500)
power_hdcov_Ghod <- rep(0, 500)
for (k in 1:500) {
  X <- matrix(0, nrow=500, ncol = 9)
  a <- c(1, 4, 7)
  b <- c(2,3,5,6,8,9)
  X <- matrix(rt(500*9, 2), nrow = 500, ncol = 9)
  sign_list <- X[, 1]*X[, 4]*X[, 7]
  plus_index <- which(sign_list<=0)
  neg_index <- which(sign_list>0)
  X[, 7][plus_index] <- -X[, 7][plus_index]
  X[, 7][neg_index] <- X[, 7][neg_index]
  X_list <- list(X[,4:6], X[, 7:9])
  
  ## rdhod
  statis <- computestatisticjdcov(X[,4:9], dim_list = rep(3, 2))
  power_rdhdcov_Ghod[k] <- length(which(emprical >= statis)) / (length(emprical) + 1)
  
  ## jdhod
  statis <- jdcov.test(X_list, stat.type = "V", alpha = 0.05)
  power_hdcov_Ghod[k] <- statis$p.value
  
  
  print(power_rdhdcov_Ghod[k])
  print(power_hdcov_Ghod[k])
  print(k)
}

mix_2_power <- data.frame(rdhdcov = power_rdhdcov_Ghod, hdcov = power_hdcov_Ghod)

write.csv(mix_2_power, "C:/Users/ziangniu/Desktop/Simulation/result/mix_TT2_power.csv")


## test the third order dependence

n <- 500
set.seed(1)
emprical <- gensamdistrhodcov(500, dim_list = rep(3, 3), niter=1000)
power_rdhdcov_Ghod_3 <- rep(0, 500)
power_hdcov_Ghod_3 <- rep(0, 500)
for (k in 1:500) {
  X <- matrix(0, nrow=500, ncol = 9)
  a <- c(1, 4, 7)
  b <- c(2,3,5,6,8,9)
  X <- matrix(rt(500*9, 2), nrow = 500, ncol = 9)
  sign_list <- X[, 1]*X[, 4]*X[, 7]
  plus_index <- which(sign_list<=0)
  neg_index <- which(sign_list>0)
  X[, 7][plus_index] <- -X[, 7][plus_index]
  X[, 7][neg_index] <- X[, 7][neg_index]
  X_list <- list(X[,1:3], X[,4:6], X[,7:9])
  
  ## rdhod
  statis <- computestatisticrdcov(X, dim_list = rep(3, 3))
  power_rdhdcov_Ghod_3[k] <- length(which(emprical >= statis)) / (length(emprical) + 1)
  
  ## jdhod
  statis <- boot_hodcov(X_list, 500)
  power_hdcov_Ghod_3[k] <- statis
  
  
  print(power_rdhdcov_Ghod_3[k])
  print(power_hdcov_Ghod_3[k])
  print(k)
}

mix_3_power <- data.frame(rdhdcov = power_rdhdcov_Ghod_3, hdcov = power_hdcov_Ghod_3)

write.csv(mix_3_power, "C:/Users/ziangniu/Desktop/Simulation/result/mix_TT3_power.csv")


## test the joint independence
n <- 500
set.seed(1)
emprical <- gensamdistrjdcov(500, dim_list = rep(3, 3), niter=1000)
power_rdhdcov_P <- rep(0, 500)
power_hdcov_P <- rep(0, 500)
for (k in 1:500) {
  X <- matrix(0, nrow=500, ncol = 9)
  a <- c(1, 4, 7)
  b <- c(2,3,5,6,8,9)
  X <- matrix(rt(500*9, 2), nrow = 500, ncol = 9)
  sign_list <- X[, 1]*X[, 4]*X[, 7]
  plus_index <- which(sign_list<=0)
  neg_index <- which(sign_list>0)
  X[, 7][plus_index] <- -X[, 7][plus_index]
  X[, 7][neg_index] <- X[, 7][neg_index]
  X_list <- list(X[,1:3], X[, 4:6], X[, 7:9])
  
  ## rdhod
  statis <- computestatisticjdcov(X, dim_list = rep(3, 3))
  power_rdhdcov_P[k] <- length(which(emprical >= statis)) / (length(emprical) + 1)
  
  statis <- jdcov.test(X_list, stat.type = "V", alpha = 0.05)
  power_hdcov_P[k] <- statis$p.value
  
  
  print(power_rdhdcov_P[k])
  print(power_hdcov_P[k])
  print(k)
}

mix_P_power <- data.frame(rdhdcov = power_rdhdcov_P, hdcov = power_hdcov_P)

write.csv(mix_P_power, "C:/Users/ziangniu/Desktop/Simulation/result/mix_TTP_power.csv")