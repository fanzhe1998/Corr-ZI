# Title: Simulations to compare the performance of association analysis methods
# analysis of microbiome-seq data. Only one taxon are simulated for illustration purpose.
# Authors: Zhe Fan (fanzhe0308@163.com)
# Date: 2022/10/26

  require(psych)
  require(compositions)
  require(pscl)
  require(pscl)
  require(maigesPack)
  require(minerva)
  

  
#####Parameters####
X <- rnorm(100000, mean = 5, sd = 1)               # metabolite overall
para_beta <- c(seq(0.1, 0.9, by = 0.1))            # covariate effect on count component
para_gamma <- c(seq(-0.9, 0.9, by = 0.1))          # covariate effect on structural zero
para_zi <- seq(0.1, 0.9, by = 0.1)                 # zero inflation rate
para_dispersion <- c(0.1, 0.2, 0.5, 1, 2, 5, 10)   # dispersion
para_samplesize <- c(50,100,250,500,750,1000,2000) # sample size
n_var <- 200
n_noi <- 400

#####################################
# Generate linear correlated data
#####################################
sigmoid = function(x) {
  1 / (1 + exp(-x))
}

Generate_data <- function(n_var, n_noi, samplesize, zi, dispersion, beta, gamma){
  
  ##empty matrix
  data_micro <- matrix(0,samplesize,n_var)
  data_met <- matrix(0,samplesize,n_var)
  
  E <- rnorm(100000, mean = 0, sd = 0.4) # random error
  MU <- exp(beta*scale(X) + E + 5)
  Y <- rnbinom(100000, size = dispersion, mu = MU) #microbe overall
  PROB <- sigmoid(gamma*scale(X) + E + 5*zi-2)
  Y_0 <- rbinom(100000, 1, prob = PROB)
  Y[Y_0 == 0] <- 0
  
  Y2 <- rnbinom(100000, size = dispersion, mu = 5) #noisy microbe
  PROB2 <- sigmoid(E + 5*zi-2)
  Y2_0 <- rbinom(100000, 1, prob = PROB2)
  Y2[Y2_0 == 0] <- 0
  


  ##Sample
  for (i in 1:iter) {
    id <- sample(1:100000, samplesize, replace = FALSE)
    x1 <- X[id]
    y1 <- Y[id]
    data_met[,i] <- x1
    data_micro[,i] <- y1
  }
  
  ##Sample noisy variables
  for (j in (n_var+1):n_noi+n_var) {  
    id <- sample(1:100000, samplesize, replace = FALSE)
    x1 <- X[id]
    data_met[,j] <- x1
    y1 <- Y2[id]
    data_micro[,j] <- y1
  }
  
  ##Compositional data and clr transformation
  data_micro <- data.frame(data_micro)
  data_micro_comp <- t(apply(data_micro, 1, function(each_row){
    the_sum <- sum(each_row, na.rm = TRUE)
    each_row <- each_row/the_sum
    return(each_row)
  }))
  data_micro_clr <- t(apply(data_micro_comp, 1, function(each_row){compositions::clr(each_row)}))
  data_met <- data.frame(apply(data_met, 2, function(each_row){scale(each_row)}))
  data_raw <- cbind(data_micro, data_met)
  data_comp <- cbind(data_micro_comp, data_met)
  data_comp_clr <- cbind(data_micro_clr, data_met)
  data_list <- list(data_raw,data_comp,data_comp_clr)
  
  return(data_list)
}

#####################################
# Simulation 1. True positive rate of different zero inflation rates in linear correlation detection, for example
par <- para_zi
#####################################
Simulation_linear <- function(n_var, n_noi, samplesize,zi ,dispersion, beta, gamma, par, n_cores = 1){
  sim_res <- list(NULL)
  par_l <- length(par)
  length(sim_res) <- par_l
  for (k in 1:p1) {
    if(par == para_zi){
      zi <- par[k]
      data_list <- Generate_data(n1,n2,samplesize,zi,dispersion,beta,gamma)
    }else{
      if(par == para_dispersion){
        dispersion <- par[k]
        data_list <- Generate_data(n1,n2,samplesize,zi,dispersion,beta,gamma)
      }else{
        if(par == para_beta){
          beta <- par[k]
          data_list <- Generate_data(n1,n2,samplesize,zi,dispersion,beta,gamma)
        }else{
          samplesize <- par[k]
          data_list <- Generate_data(n1,n2,samplesize,zi,dispersion,beta,gamma)
        }
        
      }
    }
    # Initializes the parallel environment:
    myCluster = parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(myCluster)
    
    res_p <- data.frame(NULL)
    
    res <- foreach(i = c(1:n1),.combine='rbind') %dopar% {
      
      data_raw <- data_list[[1]]
      data_comp <- data_list[[2]]
      data_comp_clr <- data_list[[3]]
      zl <- pscl::zeroinfl(formula = data_raw[,i] ~ data_raw[,n2+i], data = data_raw, dist = "negbin")
      spe <- cor.test(data_raw[,i], data_raw[,n2+i], method = "spearman")
      lr <- lm(data_raw[,n2+i] ~ data_raw[,i])
      lr_comp <- lm(data_comp_clr[,n2+i] ~ data_comp_clr[,i])
      pear <- cor.test(data_raw[,i], data_raw[,n2+i], method = "pearson")
      pear_comp <- cor.test(data_comp_clr[,i], data_comp_clr[,n2+i], method = "pearson")
      sum_zl_v <- ifelse(rownames(summary(zl)[["coefficients"]][["count"]])[2] == "Log(theta)", NA,
                         summary(zl)[["coefficients"]][["count"]][2,1])
      sum_zl_p <- ifelse(rownames(summary(zl)[["coefficients"]][["count"]])[2] == "Log(theta)", NA,
                         summary(zl)[["coefficients"]][["count"]][2,4])
      sum_zl_zv <- ifelse(rownames(summary(zl)[["coefficients"]][["zero"]])[2] == "Log(theta)", NA,
                          summary(zl)[["coefficients"]][["zero"]][2,1])
      sum_zl_zp <- ifelse(rownames(summary(zl)[["coefficients"]][["zero"]])[2] == "Log(theta)", NA,
                          summary(zl)[["coefficients"]][["zero"]][2,4])
      sum_lr <- summary(lr)
      sum_lr_comp <- summary(lr_comp)
      
      res_uni_p <- data.frame(count = sum_zl_p, 
                              spearman = spe[["p.value"]],
                              pearson_p = pear[["p.value"]],
                              pearson_comp_p = pear_comp[["p.value"]],
                              lr_p = sum_lr[["coefficients"]][2,4],
                              lr_comp_p = sum_lr_comp[["coefficients"]][2,4],
                              zero = sum_zl_zp
                              )
      res_uni_est <- data.frame(zl_count = sum_zl_v, 
                                spearman_v = spe[["estimate"]],
                                pearson_v = pear[["estimate"]],
                                pearson_comp_v = pear_comp[["estimate"]],
                                lr_v = sum_lr[["coefficients"]][2,1],
                                lr_v2 = sum_lr_comp[["coefficients"]][2,1],
                                zl_zero = sum_zl_zv
                                )
      
      return(res_uni_p)  
    }
    res_p <- rbind(res,res_p)
    sim_res[[k]] <- res_p
  }
  
  parallel::stopCluster(myCluster)
  return(sim_res)
}  


#####################################
# Compute TPR of simulation 1
#####################################
TPR_com <- function(sim_res, par){
  
  TPR <- numeric(length(par)*ncol(sim_res[[1]]))
  for (j in 1:length(par)) {
    for (i in 1:ncol(sim_res[[1]])) {
      tmp <- ifelse(sim_res[[j]][,i] < 0.05, 1, 0)
      TPR[i + (j-1)*ncol(sim_res[[1]])] <- sum(tmp[which(!is.na(tmp))])/length(tmp[which(!is.na(tmp))])
    }
  }
  TPR <- data.frame(TPR = TPR, 
                    Method = rep(c("ZI_count", "Spearman","Pearson", "Pearson2", "lr1","lr_comp","ZI_zero"), times = length(par)),
                    para = rep(par, each = ncol(sim_res[[1]])))
  return(TPR)
}
 
#####################################
# Generate nonlinear correlated data
#####################################
sigmoid = function(x) {
  1 / (1 + exp(-x))
}

Generate_nonlinear_data <- function(n_var, n_noi, samplesize, zi, dispersion, beta, gamma, par, ncores = 1,
                             a = 30, b = 10, c = 4, d = 11){
  
  ##Empty matirx
  data_micro <- matrix(0,samplesize,4*n_var+n_noi)
  data_met <- matrix(0,samplesize,n_noi)
 
  
  E <- rnorm(100000, mean = 0, sd = 0.4)  #random error
  
  MU1 <- a*scale(X)^2 + E + 3                        # Parabolic 
  Y1 <- rnbinom(100000, size = dispersion, mu = MU1) 
  MU2 <- a*sigmoid((X-5) + E)                        # sigmoid
  Y2 <- rnbinom(100000, size = dispersion, mu = MU2)
  MU3 <- b*sin(c*X + E) +d                           # sine for different Amplitude, Period and Phase
  Y3 <- rnbinom(100000, size = dispersion, mu = MU3)
  MU4 <- b*sin(X + E) + d                            # sine
  Y4 <- rnbinom(100000, size = dispersion, mu = MU4)
  MU5 <- b*cos(X + E) + d                            # cosine
  Y5 <- rnbinom(100000, size = dispersion, mu = MU5)
  
  PROB <- sigmoid(gamma*scale(X) + E + 5*zi-2)# covariate-dependent zero inflation rate
  Y_0 <- rbinom(100000, 1, prob = PROB)#
  Y1[Y_0 == 0] <- 0
  Y2[Y_0 == 0] <- 0
  Y3[Y_0 == 0] <- 0
  Y4[Y_0 == 0] <- 0
  Y5[Y_0 == 0] <- 0
  
  Y02 <- rnbinom(100000, size = dispersion, mu = (exp(6) + E))
  PROB02 <- sigmoid(E + zi)
  Y02_0 <- rbinom(100000, 1, prob = PROB02)
  Y02[Y02_0 == 0] <- 0
  
  ##Sample
  for (i in 1:n_var) {
    id <- sample(1:100000, samplesize, replace = FALSE)
    x1 <- X[id]##代谢物数据
    y1 <- Y1[id]
    y2 <- Y2[id]
    y3 <- Y3[id]
    y4 <- Y4[id]
    y5 <- Y5[id]
    data_met[,i] <- x1
    data_micro[ ,5*(i-1)+1] <- y1
    data_micro[ ,5*(i-1)+2] <- y2
    data_micro[ ,5*(i-1)+3] <- y3
    data_micro[ ,5*(i-1)+4] <- y4
    data_micro[ ,5*(i-1)+5] <- y5
  }
  
  # Sample noisy varaibles
  for (j in (n_var+1):4*n_var + n_noi) {  
    id <- sample(1:100000, samplesize, replace = FALSE)
    x1 <- X[id]
    data_met[,j] <- x1
    y1 <- Y02[id]
    data_micro[,4*n1+j] <- y1
  }
  
  ##
  data_micro <- data.frame(data_micro)
  data_micro_c <- data.frame(t(apply(data_micro, 1, function(each_row){
    the_sum <- sum(each_row, na.rm = TRUE)
    each_row <- each_row/the_sum
    return(each_row)
  })))
  data_micro_clr <- data.frame(apply(data_micro_c, 2, function(each_row){compositions::clr(each_row)}))
  data_met <- data.frame(t(apply(data_met, 1, function(each_row){scale(each_row)})))
  data_raw <- cbind(data_micro, data_met)
  data_comp <- cbind(data_micro_c, data_met)
  data_comp_clr <- cbind(data_micro_clr, data_met)
  data_list <- list(data_raw,data_comp,data_comp_clr)
  return(data_list)
}

#####################################
# Simulation 2. True positive rate of different zero inflation rates in nonlinear correlation detection, for example
par <- para_zi
#####################################
Simulation_nonlinear <- function(n_var, n_noi, samplesize,zi ,dispersion, beta, gamma, par, n_cores = 1){
  sim_res <- list(NULL)
  par_l <- length(par)
  length(sim_res) <- par_l
  for (k in 1:p1) {
    if(par == para_zi){
      zi <- par[k]
      data_list <- Generate_nonlinear_data(n_var, n_noi, samplesize, zi, dispersion, beta, gamma, par, ncores = 1,
                                           a = 30, b = 10, c = 4, d = 11)
    }else{
      if(par == para_dispersion){
        dispersion <- par[k]
        data_list <- Generate_nonlinear_data(n_var, n_noi, samplesize, zi, dispersion, beta, gamma, par, ncores = 1,
                                             a = 30, b = 10, c = 4, d = 11)
      }else{
        if(par == para_beta){
          beta <- par[k]
          data_list <- Generate_nonlinear_data(n_var, n_noi, samplesize, zi, dispersion, beta, gamma, par, ncores = 1,
                                               a = 30, b = 10, c = 4, d = 11)
        }else{
          samplesize <- par[k]
          data_list <- Generate_nonlinear_data(n_var, n_noi, samplesize, zi, dispersion, beta, gamma, par, ncores = 1,
                                               a = 30, b = 10, c = 4, d = 11)
        }
        
      }
    }


    cl <- makeCluster(cores)
    registerDoParallel(cl)
    res_n1 <- data.frame(NULL)
    res <- foreach(i = c(1:n1),.combine='rbind') %dopar% {
      
      require(pscl)
      require(maigesPack)
      require(minerva)
      source("MIC.R", encoding = "utf-8")
      
      data_raw <- data_list[[1]]
      data_comp <- data_list[[2]]
      data_comp_clr <- data_list[[3]]
      
      mic1 <- MIC(data_raw[,4*n1+n2+i], data_raw[,1+5*(i-1)])
      MI1 <- bootstrapMI(data_raw[,4*n1+n2+i], data_raw[,1+5*(i-1)], bRep = brep, ret = "p-value")
      
      mic2 <- MIC(data_raw[,4*n1+n2+i], data_raw[,2+5*(i-1)])
      MI2 <- bootstrapMI(data_raw[,4*n1+n2+i], data_raw[,2+5*(i-1)], bRep = brep, ret = "p-value")
      
      mic3 <- MIC(data_raw[,4*n1+n2+i], data_raw[,3+5*(i-1)])
      MI3 <- bootstrapMI(data_raw[,4*n1+n2+i], data_raw[,3+5*(i-1)], bRep = brep, ret = "p-value")
      
      mic4 <- MIC(data_raw[,4*n1+n2+i], data_raw[,4 + 5*(i-1)])
      MI4 <- bootstrapMI(data_raw[,4*n1+n2+i], data_raw[,4 + 5*(i-1)], bRep = brep, ret = "p-value")
      
      mic5 <- MIC(data_raw[,4*n1+n2+i], data_raw[,5+5*(i-1)])
      MI5 <- bootstrapMI(data_raw[,4*n1+n2+i], data_raw[,5+5*(i-1)], bRep = brep, ret = "p-value")
      
      
      res_uni_p <- data.frame(mic_pvalue1 = mic1$'p-value',
                              MI_v1 = MI1,
                              
                              mic_pvalue2 = mic2$'p-value',
                              MI_v2 = MI2,
                              
                              mic_pvalue3 = mic3$'p-value',
                              MI_v3 = MI3,
                              
                              mic_pvalue4 = mic4$'p-value',
                              MI_v4 = MI4,
                              
                              mic_pvalue5 = mic5$'p-value',
                              MI_v5 = MI5
                              )
      
      return(res_uni_p)  
    }
    res_p <- rbind(res,res_p)
    sim_res[[k]] <- res_p
  }
  
  stopCluster(cl)
  return(sim_res)
}  

TPR_com_nonlinear <- function(sim_res, par){
  
  TPR <- numeric(length(par)*ncol(sim_res[[1]]))
  for (j in 1:length(par)) {
    for (i in 1:ncol(sim_res[[1]])) {
      tmp <- ifelse(sim_res[[j]][,i] < 0.05, 1, 0)
      TPR[i + (j-1)*ncol(sim_res[[1]])] <- sum(tmp[which(!is.na(tmp))])/length(tmp[which(!is.na(tmp))])
    }
  }
  TPR <- data.frame(TPR = TPR, 
                    Method = rep(c("MIC","MI"), times = length(par)),
                    para = rep(par, each = ncol(sim_res[[1]])))
  return(TPR)
}


