## data simulation
##########
library(doParallel)
detectCores()
cl <- makeCluster(40)
registerDoParallel(cl)
stopCluster(cl)
#########
MIC <- function(x,y,R=100,...) {
  #MIC
  mic<-mine(x=x,y=y,...)
  
  #For calculate p-value de MIC
  if (! is.null(R)) {
    R <- floor(R)
    if (R < 1) R <- 100
  } else {
    R <- 100
  }
  Rep<-as.data.frame(rep(as.data.frame(y),R))
  Rep2<-as.data.frame(apply(Rep,2,sample))
  permic<-matrix(NA,nrow=R,ncol=7)
  colnames(permic)<-c("MIC","MAS","MEV","MCN","MIC-R2", "GMIC","TIC")
  for (i in 1:R){
    p<-mine(x=x,y=Rep2[,i],...)
    permic[i,1:7]<-c(p$MIC,p$MAS,p$MEV,p$MCN,p$`MIC-R2`,p$GMIC,p$TIC)
  }
  
  permic<-as.data.frame(permic)
  pvalor<-nrow(permic[which(permic$MIC>=mic$MIC),])/nrow(permic)
  
  mic2<-as.data.frame(cbind(mic$MIC,pvalor))
  colnames(mic2)<-c("MIC","p-value")
  return(mic2)
}
sigmoid = function(x) {
  1 / (1 + exp(-x))
}
#####Parameters####
X <- rnorm(100000, mean = 5, sd = 1)##metabolite overall
para_beta <- c(seq(0.1, 0.9, by = 0.1))##covariate effect on count component
para_gamma <- c(seq(-0.9, 0.9, by = 0.1))##covariate effect on structural zero
para_zi <- seq(0.1, 0.9, by = 0.1)##zero inflation rate
para_dispersion <- c(0.1, 0.2, 0.5, 1, 2, 5, 10)##over-dispersion
para_samplesize <- c(50,100,250,500,750,1000,2000)##sample size
summary(X)
#########

## data generate
gedata <- function(num_vars1, n2, samplesize, zi, dispersion, beta, gamma){
  
  require(psych)
  require(compositions)
  
  ##empty matrix
  data_micro <- matrix(0,samplesize,n2)
  data_met <- matrix(0,samplesize,n2)
  
  E <- rnorm(100000, mean = 0, sd = 0.4)
  MU <- exp(beta*scale(X) + E + 5)
  Y <- rnbinom(100000, size = dispersion, mu = MU)
  PROB <- sigmoid(gamma*scale(X) + E + 5*zi-2)
  Y_0 <- rbinom(100000, 1, prob = PROB)
  Y[Y_0 == 0] <- 0
  
  Y2 <- rnbinom(100000, size = dispersion, mu = 5)
  PROB2 <- sigmoid(E + 5*zi-2)
  Y2_0 <- rbinom(100000, 1, prob = PROB2)
  Y2[Y2_0 == 0] <- 0
  
  
  ##sample
  for (i in 1:n1) {
    id <- sample(1:100000, samplesize, replace = FALSE)
    x1 <- X[id]
    y1 <- Y[id]
    data_met[,i] <- x1
    data_micro[,i] <- y1
  }
  
  ##no covariate effect
  for (j in (n1+1):n2) {  
    id <- sample(1:100000, samplesize, replace = FALSE)
    x1 <- X[id]
    data_met[,j] <- x1
    y1 <- Y2[id]
    data_micro[,j] <- y1
  }
  
  ##compositional data and clr transformation
  data_micro <- data.frame(data_micro)
  data_micro_c <- t(apply(data_micro, 1, function(each_row){
    the_sum <- sum(each_row, na.rm = TRUE)
    each_row <- each_row/the_sum
    return(each_row)
  }))
  data_micro_clr <- t(apply(data_micro_c, 1, function(each_row){compositions::clr(each_row)}))
  data_met <- data.frame(apply(data_met, 2, function(each_row){scale(each_row)}))
  data_raw <- cbind(data_micro, data_met)
  data_comp <- cbind(data_micro_c, data_met)
  data_comp_clr <- cbind(data_micro_clr, data_met)
  data_list <- list(data_raw,data_comp,data_comp_clr)
  return(data_list)
}

simu2 <- function(i,j,n,data_list){
  
  require(pscl)
  require(maigesPack)
  require(minerva)
  
  data_raw <- data_list[[1]]
  data_comp <- data_list[[2]]
  data_comp_clr <- data_list[[3]]
  zl <- zeroinfl(formula = data_raw[,i] ~ data_raw[,n+j], data = data_raw, dist = "negbin")
  spe <- cor.test(data_raw[,i], data_raw[,n+j], method = "spearman")
  lr <- lm(data_raw[,n+j] ~ data_raw[,i])
  lr_comp <- lm(data_comp_clr[,n+j] ~ data_comp_clr[,i])
  pear <- cor.test(data_raw[,i], data_raw[,n+j], method = "pearson")
  pear_comp <- cor.test(data_comp_clr[,i], data_comp_clr[,n+j], method = "pearson")
  sum_zl <- summary(zl)
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
  
  res_si <- data.frame(zl_count = sum_zl_v, 
                       spearman_v = spe[["estimate"]],
                       pearson_v = pear[["estimate"]],
                       pearson_comp_v = pear_comp[["estimate"]],
                       lr_v = sum_lr[["coefficients"]][2,1],
                       lr_v2 = sum_lr_comp[["coefficients"]][2,1],
                       zl_zero = sum_zl_zv,
                       
                       count = sum_zl_p, 
                       spearman = spe[["p.value"]],
                       pearson_p = pear[["p.value"]],
                       pearson_comp_p = pear_comp[["p.value"]],
                       lr_p = sum_lr[["coefficients"]][2,4],
                       lr_comp_p = sum_lr_comp[["coefficients"]][2,4],
                       zero = sum_zl_zp
                       
  )
  return(res_si)  
}

## 
Simulation_compare <- function(n1,n2,samplesize,zi,dispersion,beta,gamma,par,n_cores = 1){
  test <- list(NULL)
  p1 <- length(par)
  length(test) <- p1
  for (k in 1:p1) {
    if(par == para_zi){
      zi <- par[k]
      data_list <- gedata(n1,n2,samplesize,zi,dispersion,beta,gamma)
    }else{
      if(par == para_dispersion){
        dispersion <- par[k]
        data_list <- gedata(n1,n2,samplesize,zi,dispersion,beta,gamma)
      }else{
        if(par == para_beta){
          beta <- par[k]
          data_list <- gedata(n1,n2,samplesize,zi,dispersion,beta,gamma)
        }else{
          samplesize <- par[k]
          data_list <- gedata(n1,n2,samplesize,zi,dispersion,beta,gamma)
        }
        
      }
    }
    # Initializes the parallel environment:
    myCluster = parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(myCluster)

    res_n1 <- data.frame(NULL)
    
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
      
      res_uni <- data.frame(zl_count = sum_zl_v, 
                            spearman_v = spe[["estimate"]],
                            pearson_v = pear[["estimate"]],
                            pearson_comp_v = pear_comp[["estimate"]],
                            lr_v = sum_lr[["coefficients"]][2,1],
                            lr_v2 = sum_lr_comp[["coefficients"]][2,1],
                            zl_zero = sum_zl_zv,
                            
                            count = sum_zl_p, 
                            spearman = spe[["p.value"]],
                            pearson_p = pear[["p.value"]],
                            pearson_comp_p = pear_comp[["p.value"]],
                            lr_p = sum_lr[["coefficients"]][2,4],
                            lr_comp_p = sum_lr_comp[["coefficients"]][2,4],
                            zero = sum_zl_zp
                            
      )
      return(res_uni)  
    }
    res_n1 <- rbind(res,res_n1)
    test[[k]] <- res_n1
  }
  
  power <- numeric(p1*8)
  for (j in 1:p1) {
    for (i in 8:14) {
      a <- ifelse(test[[j]][,i] < 0.05, 1, 0)
      power[(i-7) + (j-1)*8] <- sum(a[which(!is.na(a))])/length(a[which(!is.na(a))])
    }
    a <- ifelse(test[[j]][,8] | test[[j]][,14] < 0.05, 1, 0)
    power[8 + (j-1)*8] <- sum(a[which(!is.na(a))])/length(a[which(!is.na(a))])
  }
  plot1 <- data.frame(Power = power, Method = rep(c("ZI_count", "Spearman","Pearson", "Pearson2",
                                                    "lr1","lr_comp","ZI_zero","ZI"), times = p1),
                      para = rep(par, each = 8))
  
  re <- list(plot1, test)
  parallel::stopCluster(myCluster)
  return(re)
}

alpha_size <- result(200,400,500,0,1,0,0,p_zero0)  


## 
gedata_nonlinear <- function(n1, n2, samplesize, zi, dispersion, beta, gamma,
                      a = 30, b = 10, c = 4, d = 11){
  
  
  require(pscl)
  require(psych)
  require(compositions)
  
  ##空矩阵
  data_micro <- matrix(0,samplesize,4*n1+n2)
  data_met <- matrix(0,samplesize,n2)
  
  
  E <- rnorm(100000, mean = 0, sd = 0.5)
  E2 <- rnorm(100000, mean = 0, sd = 0.2)
  
  MU1 <- a*scale(X)^2 + E + 3
  Y1 <- rnbinom(100000, size = dispersion, mu = MU1)
  MU2 <- a*sigmoid((X-5) + E)
  Y2 <- rnbinom(100000, size = dispersion, mu = MU2)
  MU3 <- b*sin(c*X + E) +d
  Y3 <- rnbinom(100000, size = dispersion, mu = MU3)
  MU4 <- b*sin(X + E) + d
  Y4 <- rnbinom(100000, size = dispersion, mu = MU4)
  MU5 <- b*cos(X + E) + d
  Y5 <- rnbinom(100000, size = dispersion, mu = MU5)
  
  PROB <- sigmoid(gamma*scale(X) + E + 5*zi-2)##零膨胀比例
  Y_0 <- rbinom(100000, 1, prob = PROB)##零膨胀比例
  Y1[Y_0 == 0] <- 0
  Y2[Y_0 == 0] <- 0
  Y3[Y_0 == 0] <- 0
  Y4[Y_0 == 0] <- 0
  Y5[Y_0 == 0] <- 0
  
  Y02 <- rnbinom(100000, size = dispersion, mu = (exp(6) + E))
  PROB02 <- sigmoid(E + zi)##零膨胀比例
  Y02_0 <- rbinom(100000, 1, prob = PROB02)##零膨胀比例
  Y02[Y02_0 == 0] <- 0
  
  ##抽出n1对相关变量
  for (i in 1:n1) {
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
  
  ##抽出n2-n1对不相关变量
  for (j in (n1+1):n2) {  
    id <- sample(1:100000, samplesize, replace = FALSE)
    x1 <- X[id]##代谢物数据
    data_met[,j] <- x1
    y1 <- Y02[id]
    data_micro[,4*n1+j] <- y1
  }
  
  ##微生物中心对数转化
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

simu_n2 <- function(i,j,n,data_list){
  
  require(pscl)
  require(maigesPack)
  require(minerva)
  
  data_raw <- data_list[[1]]
  data_comp <- data_list[[2]]
  data_comp_clr <- data_list[[3]]
  
  require(pscl)
  require(maigesPack)
  require(minerva)
  
  # data_raw <- data_list[[1]]
  # data_comp <- data_list[[2]]
  
  mic1 <- MIC(data_raw[,n+j], data_raw[,1+5*(i-1)])
  spe1 <- cor.test(data_raw[,n+j], data_raw[,1+5*(i-1)], method = "spearman")
  MI1 <- bootstrapMI(data_raw[,n+j], data_raw[,1+5*(i-1)], bRep = rep, ret = "p-value")
  
  mic2 <- MIC(data_raw[,n+j], data_raw[,2+5*(i-1)])
  spe2 <- cor.test(data_raw[,n+j], data_raw[,2+5*(i-1)], method = "spearman")
  MI2 <- bootstrapMI(data_raw[,n+j], data_raw[,2+5*(i-1)], bRep = rep, ret = "p-value")
  
  mic3 <- MIC(data_raw[,n+j], data_raw[,3+5*(i-1)])
  spe3 <- cor.test(data_raw[,n+j], data_raw[,3+5*(i-1)], method = "spearman")
  MI3 <- bootstrapMI(data_raw[,n+j], data_raw[,3+5*(i-1)], bRep = rep, ret = "p-value")
  
  mic4 <- MIC(data_raw[,n+j], data_raw[,4 + 5*(i-1)])
  spe4 <- cor.test(data_raw[,n+j], data_raw[,4 + 5*(i-1)], method = "spearman")
  MI4 <- bootstrapMI(data_raw[,n+j], data_raw[,4 + 5*(i-1)], bRep = rep, ret = "p-value")
  
  mic5 <- MIC(data_raw[,n+j], data_raw[,5+5*(i-1)])
  spe5 <- cor.test(data_raw[,n+j], data_raw[,5+5*(i-1)], method = "spearman")
  MI5 <- bootstrapMI(data_raw[,n+j], data_raw[,5+5*(i-1)], bRep = rep, ret = "p-value")
  
  res_si <- data.frame(mic_pvalue1 = mic1$'p-value',
                       spearman1 = spe1[["p.value"]],
                       MI_v1 = MI1,
                       
                       mic_pvalue2 = mic2$'p-value',
                       spearman2 = spe2[["p.value"]],
                       MI_v2 = MI2,
                       
                       mic_pvalue3 = mic3$'p-value',
                       spearman3 = spe3[["p.value"]],
                       MI_v3 = MI3,
                       
                       mic_pvalue4 = mic4$'p-value',
                       spearman4 = spe4[["p.value"]],
                       MI_v4 = MI4,
                       
                       mic_pvalue5 = mic5$'p-value',
                       spearman5 = spe5[["p.value"]],
                       MI_v5 = MI5,
                       
                       MIC1 = mic1$'MIC',
                       spearman_v1 = spe1[["estimate"]],
                       
                       
                       MIC2 = mic2$'MIC',
                       spearman_v2 = spe2[["estimate"]],
                       
                       MIC3 = mic3$'MIC',
                       spearman_v3 = spe3[["estimate"]],
                       
                       MIC4 = mic4$'MIC',
                       spearman_v4 = spe4[["estimate"]],
                       
                       MIC5 = mic5$'MIC',
                       spearman_v5 = spe5[["estimate"]]
  )
  return(res_si)  
}

##
result_n2 <- function(n1,n2,samplesize,zi,dispersion,beta,gamma,par,
                      cores,){
  test <- list(NULL)
  p1 <- length(par)
  length(test) <- p1
  for (k in 1:p1) {
    if(par == para_zi){
      zi <- par[k]
      data_list <- gedata_n2(n1, n2, samplesize, zi, dispersion, beta, gamma,
                           a = 30, b = 10, c = 4, d = 11)
    }else{
      if(par == para_dispersion){
        dispersion <- par[k]
        data_list <- gedata_n2(n1, n2, samplesize, zi, dispersion, beta, gamma,
                             a = 30, b = 10, c = 4, d = 11)
      }else{
        if(par == para_beta){
          beta <- par[k]
          data_list <- gedata_n2(n1, n2, samplesize, zi, dispersion, beta, gamma,
                               a = 30, b = 10, c = 4, d = 11)
        }else{
          samplesize <- par[k]
          data_list <- gedata_n2(n1, n2, samplesize, zi, dispersion, beta, gamma,
                               a = 30, b = 10, c = 4, d = 11)
        }
        
      }
    }
    
    require(doParallel)
    detectCores()
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
      spe1 <- cor.test(data_raw[,4*n1+n2+i], data_raw[,1+5*(i-1)], method = "spearman")
      MI1 <- bootstrapMI(data_raw[,4*n1+n2+i], data_raw[,1+5*(i-1)], bRep = brep, ret = "p-value")
      
      mic2 <- MIC(data_raw[,4*n1+n2+i], data_raw[,2+5*(i-1)])
      spe2 <- cor.test(data_raw[,4*n1+n2+i], data_raw[,2+5*(i-1)], method = "spearman")
      MI2 <- bootstrapMI(data_raw[,4*n1+n2+i], data_raw[,2+5*(i-1)], bRep = brep, ret = "p-value")
      
      mic3 <- MIC(data_raw[,4*n1+n2+i], data_raw[,3+5*(i-1)])
      spe3 <- cor.test(data_raw[,4*n1+n2+i], data_raw[,3+5*(i-1)], method = "spearman")
      MI3 <- bootstrapMI(data_raw[,4*n1+n2+i], data_raw[,3+5*(i-1)], bRep = brep, ret = "p-value")
      
      mic4 <- MIC(data_raw[,4*n1+n2+i], data_raw[,4 + 5*(i-1)])
      spe4 <- cor.test(data_raw[,4*n1+n2+i], data_raw[,4 + 5*(i-1)], method = "spearman")
      MI4 <- bootstrapMI(data_raw[,4*n1+n2+i], data_raw[,4 + 5*(i-1)], bRep = brep, ret = "p-value")
      
      mic5 <- MIC(data_raw[,4*n1+n2+i], data_raw[,5+5*(i-1)])
      spe5 <- cor.test(data_raw[,4*n1+n2+i], data_raw[,5+5*(i-1)], method = "spearman")
      MI5 <- bootstrapMI(data_raw[,4*n1+n2+i], data_raw[,5+5*(i-1)], bRep = brep, ret = "p-value")
      
      res_uni <- data.frame(mic_pvalue1 = mic1$'p-value',
                           spearman1 = spe1[["p.value"]],
                           MI_v1 = MI1,
                           
                           mic_pvalue2 = mic2$'p-value',
                           spearman2 = spe2[["p.value"]],
                           MI_v2 = MI2,
                           
                           mic_pvalue3 = mic3$'p-value',
                           spearman3 = spe3[["p.value"]],
                           MI_v3 = MI3,
                           
                           mic_pvalue4 = mic4$'p-value',
                           spearman4 = spe4[["p.value"]],
                           MI_v4 = MI4,
                           
                           mic_pvalue5 = mic5$'p-value',
                           spearman5 = spe5[["p.value"]],
                           MI_v5 = MI5,
                           
                           MIC1 = mic1$'MIC',
                           spearman_v1 = spe1[["estimate"]],
                           
                           
                           MIC2 = mic2$'MIC',
                           spearman_v2 = spe2[["estimate"]],
                           
                           MIC3 = mic3$'MIC',
                           spearman_v3 = spe3[["estimate"]],
                           
                           MIC4 = mic4$'MIC',
                           spearman_v4 = spe4[["estimate"]],
                           
                           MIC5 = mic5$'MIC',
                           spearman_v5 = spe5[["estimate"]]
      )
      return(res_uni)  
    }
    res_n1 <- rbind(res,res_n1)
    test[[k]] <- res_n1
  
  }
  power <- numeric(p1*15)
  for (j in 1:p1) {
    for (i in 1:15) {
      a <- ifelse(test[[j]][,i] < 0.05, 1, 0)
      power[i + (j-1)*15] <- sum(a[which(!is.na(a))])/length(a[which(!is.na(a))])
    }
  }
  plot1 <- data.frame(Power = power, Method = rep(c("MIC", "Spearman","KNN-MI"), times = p1*5),
                      para = rep(par, each = 15))
  plot1 <- aggregate(plot1$Power, by=list(plot1$Method, plot1$para), mean)
  colnames(plot1) <- c("Method", "para", "Power")
  #plot1 <- plot1[!(plot1$Method == "lr1" | plot1$Method == "lr2"),]
  re <- list(plot1, test)
  stopCluster(cl)
  return(re)
}


