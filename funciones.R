# if (!require("pacman")) install.packages("pacman")
# pacman::p_load(PropCIs, lme4, dplyr, msm, pracma, plotly, rjags, tibble, forcats)

library(PropCIs)
library(lme4)
library(dplyr)
library(msm)
library(pracma)
library(plotly)
library(rjags)
library(tibble)
library(forcats)



########
######## Format data table
######## 

formato_datos <- function(datos){
  Y = bind_rows(
    datos %>% mutate(
      sens = 0,
      n = tp + fn,
      true = tp,
      study = 1:nrow(datos)
    ),
    datos %>% mutate(
      sens = 1,
      n = fp + tn,
      true = tn,
      study = 1:nrow(datos)
    )
  )
  Y = Y[order(Y$id), ]
  Y$spec <- 1 - Y$sens
  return(Y)
}

#Y <- formato_datos(datos)



########
######## Format table for forest plots of sen and spec
########


format_forest <- function(X) {
  Y = X[order(X$id), ]
  Y = X %>% mutate(
    sens = tp/n1,
    spec = tn/n0
  )
  for (i in 1:nrow(Y)){
    tmp = exactci(Y$tp[i], Y$n1[i], 0.95)
    Y$CIinfsen[i] = round(tmp$conf.int[1],3)
    Y$CIsupsen[i] = round(tmp$conf.int[2],3)
    tmp2 = exactci(Y$tn[i], Y$n0[i], 0.95)
    Y$CIinfspe[i] = round(tmp2$conf.int[1],3)
    Y$CIsupspe[i] = round(tmp2$conf.int[2],3)
  }
    return(Y)
}

#Y_forest <-  format_forest(X)


########
######## Forest plot CI sen and spec
########


forest_sens <- function(Y_forest) {
  sens <- Y_forest %>% mutate(type = "Sensitivity", labels = "Sensitivity") %>% ggplot(aes(
    y = sens,
    fct_reorder(id, sens),
    ymin = CIinfsen,
    ymax = CIsupsen
  )) +  geom_pointrange(shape = 15, size = 0.5) + coord_flip() + theme_bw() +
    xlab("Study ID") + ylab("Sensitivity (95% confidence interval)")
  return(sens)
}

#CI_sens <- forest_sens(Y)

forest_spec <- function(Y_forest) {
  spec <- Y_forest %>% mutate(type = "Specificity", labels = "Specificity") %>% ggplot(aes(
    y = spec,
    fct_reorder(id, sens),
    ymin = CIinfspe,
    ymax = CIsupspe
  )) +  geom_pointrange(shape = 15, size = 0.5) +  coord_flip() + theme_bw() +
    xlab("Study ID") + ylab("Specificity (95% confidence interval)")
  return(spec)
}

#CI_spec <- forest_spec(Y)




########
######## GLM
######## 

modelo <- function(X){
  X_biv <-  reshape(X, direction = "long", varying = list(c("n1", "n0"), c("true1", "true0")), timevar = "sens", times = c(1,0),
                    v.names = c("n", "true"))
  X_biv %>% remove_rownames()
  
  X_biv <- X_biv[order(X_biv$id),]
  X_biv$spec <-  1-X_biv$sens
  
  # Meta-analysis:
  MA_biv <-  glmer(formula = cbind(  true , n - true ) ~ 0 + sens + spec + (0+sens + spec|study), data = X_biv, family = binomial, 
                   nAGQ = 1, verbose = 0)
  # summary of the model:
  ma_biv = summary(MA_biv)
  
  # logit sens and logit spec
  lsens = ma_biv$coeff[1,1]
  lspec = ma_biv$coeff[2,1]
  
  
  se.lsens = ma_biv$coeff[1,2]
  se.lspec = ma_biv$coeff[2,2]
  
  # 95% confidence intervals for logit sens and logit spec
  Sens = c(lsens, lsens-qnorm(0.975)*se.lsens, lsens+qnorm(0.975)*se.lsens)
  Spec = c(lspec, lspec-qnorm(0.975)*se.lspec, lspec+qnorm(0.975)*se.lspec)
  
  # sens and spec estimates in the raw scale
  s <-plogis( Sens ) 
  e <-plogis( Spec ) 
  
  # DOR and likelihood ratios
  DOR = exp(lsens+lspec )
  LRp = plogis(lsens)/(1-plogis(lspec))
  LRm = ((1-plogis(lsens))/plogis(lspec))
  
  # Confidence intervals with the delta method
  se.DOR = deltamethod (~ exp(x1+x2) , mean = c(lsens,lspec) , cov = ma_biv$vcov )
  
  se.LRp = deltamethod (~ (exp(x1)/(1+exp(x1)))/(1-(exp(x2)/(1+exp(x2)))) , 
                        mean = c(lsens,lspec) , cov = ma_biv$vcov )
  
  se.LRm = deltamethod (~ (1-(exp(x1)/(1+exp(x1))))/(exp(x2)/(1+exp(x2))) , 
                        mean = c(lsens,lspec) , cov = ma_biv$vcov ) 
  
  # random effects correlation
  corr <- attr(summary(MA_biv)$varcor$study,"correlation")[1,2]

  tablaOutPut = data.frame(
    Coefficient = c("Sensitivity", "Specificity", "DOR", "LR+" , "LR-", "RE Correlation"),
    Estimate = c(s[1], e[1], DOR , LRp , LRm, corr) ,
    LCI = c(
      round(s[2], 3),
      round(e[2], 3),
      round(DOR - qnorm(0.975) * se.DOR, 3) ,
      round(LRp - qnorm(0.975) * se.LRp, 3) ,
      round(LRm - qnorm(0.975) * se.LRm, 3),
      ""
    ),
    UCI = c(
      round(s[3], 3),
      round(e[3], 3),
      round(DOR + qnorm(0.975) * se.DOR,3) ,
      round(LRp + qnorm(0.975) * se.LRp,3) ,
      round(LRm + qnorm(0.975) * se.LRm,3),
      ""
    )
  )
  
  return(list(tabla=tablaOutPut, summary=ma_biv, MA_biv = MA_biv, corr = corr))
}

#ma_biv = modelo(X)$summary
#modelo(X)$tabla
#round(modelo(X)$corr, 3)
#MA_biv <- modelo(X)$MA_biv




########
######## SROC plane plot
######## 


plot_sroc <- function(Y, ma_biv){
  phi <- seq(0, 2 * pi, len = 10000)
  
  ## Extract the coefficients of the model ##
  
  u = ma_biv$coeff[1, 1]           ## E(logit(sen))
  v = ma_biv$coeff[2, 1]           ## E(logit(spe)
  SELogitSen = ma_biv$coeff[1, 2]  ## standard error E(logit(sen))
  SELogitSpe = ma_biv$coeff[2, 2]  ## standard error E(logit(spe))
  covAB = ma_biv$vcov[1, 2]        ##Covariance betwen E(logit(sen)) and E(logit(spe))
  
  VarLogitSen = ma_biv$varcor$study[1, 1]   ##Var(logit(sen))
  VarLogitSpe = ma_biv$varcor$study[2, 2]  ##Var(logit(spe))
  corAB = ma_biv$varcor$study[2]
  N = ma_biv$ngrps
  
  SEPA <- sqrt(VarLogitSen + SELogitSen ^ 2)
  SEPB <- sqrt(VarLogitSpe + SELogitSpe ^ 2)
  sAB <- corAB * sqrt(VarLogitSen) * sqrt(VarLogitSpe)
  rho <- (sAB + covAB) / (SEPA * SEPB)
  
  # boundary constant of the ellipse
  croot <- sqrt(2 * qf(0.95, 2, N - 2))
  
  a <- SEPA * croot
  b <- SEPB * croot
  sen <- exp(u) / (1 + exp(u))
  esp <- exp(v) / (1 + exp(v))
  
  y <-
    u + a * cos(phi + acos(pmin(pmax(rho,-1.0),1.0)))     ## Sensitivity (logit scale)
  x <-
    v + b * cos(phi)                 ## Specificity (logit scale)
  specificity <- 1 - (exp(x) / (1 + exp(x)))  ## ROC scale
  sensitivity <- exp(y) / (1 + exp(y))      ## ROC scale
  
  
  area <- polyarea(specificity, sensitivity)
  area <- round(area, 4)
  
  
  tabla = data.frame(sensitivity, specificity)
  

  
  sroc <- plot_ly(
    data = Y %>% mutate(x = 1 - (tn / (tn + fp)), y = tp / (tp + fn)),
    x =  ~ x,
    y =  ~ y,
    type = "scatter",
    hoverinfo = "text",
    text = ~ paste('ID:', id, '\n', '1-Specificity:', x, '\n', 'Sensitivity:' , y),
    marker = list(color = pal_npg("nrc")(1)),
    width = 450, 
    height = 450
    ) %>% 
    add_trace(
      data = tabla,
      x = ~ specificity,
      y = ~ sensitivity,
      text = NULL,
      mode = "lines",
      marker = list(size = 0.5, color = pal_npg("nrc")(10))
    ) %>% 
    layout(
      title = "SROC plane",
      xaxis = list(title = "1-Specificity", range = c(-0.05, 1.05)),
      yaxis = list(title = "Sensitivity", range = c(-0.05, 1.05)),
      showlegend = FALSE
    )
  
  tablaRevMan = data.frame(
    Coefficient = c("E(logitSe)", "E(logitSp)", "Var(logitSe)", "Var(logitSp)" ,  "SE(E(logitSe))",
                    "SE(E(logitSp))", "Cov(Es))", "Studies"),
    Estimate = c(u, v, VarLogitSen, VarLogitSpe, SELogitSen, SELogitSpe, covAB, nlevels(unique(Y$id))))
  
  
  return(list(sroc = sroc, area = area, tabla = tablaRevMan, tabla_ellipse = tabla))
}

#sroc <- plot_sroc(Y, ma_biv)$sroc
#area <- plot_sroc(Y, ma_biv)$area





########
######## Modelo Bayes I2-bivariante (heterogeneity)
########

heterogeneidad <- function(X){
  modelo.string <- "model  {
      for (i in 1:K){
  #Likelihood
  X1[i]~dbin(p[i,1],N1[i])
  X2[i]~dbin(p[i,2],N2[i])
  logit(p[i,1]) <- theta[i,1]
  logit(p[i,2]) <- theta[i,2]
  
  # reciprocal of sample size in ith study
  inv.N1[i] <- 1/N1[i]
  inv.N2[i] <- 1/N2[i]
  
  # Level-1 priors
  theta[i,1]~dnorm(thetaA, precision.thetaA)
  theta[i,2]~dnorm(theta2[i], precision.theta2)
  theta2[i] <- lambda0 + lambda1*(theta[i,1]-mean(theta[,1]))
  # See text for note on prior specification
      }

  # Level-2 priors
  thetaA~dunif(-5,5)
  precision.thetaA~dgamma(0.5,0.5)
  precision.theta2~dgamma(0.5,0.5)
  lambda1~dunif(-10,10)
  lambda0~dunif(-10,10)
  
  # Pooled sensitivity and specificity
  thetaB <- lambda0 + lambda1*(thetaA - mean(theta[,1]))
  S <- 1/(1+exp(-thetaA))
  C <- 1/(1+exp(-thetaB))
  
  # Between-study variance
  tausq[1,1] <- 1/precision.thetaA
  
  tausq[2,2] <- 1/precision.theta2+1/precision.thetaA*(lambda1*lambda1)
  
  # Between-study correlation
  rho<- lambda1*sqrt(tausq[1,1])/sqrt(1/precision.theta2+ tausq[1,1]*lambda1*lambda1)
  
  # Between-study covariance
  tausq[2,1] <- lambda1*tausq[1,1]
  tausq[1,2] <- tausq[2,1]
  
  # === I-squared E calculation ===
  # mean of inverse sample size in each study
  avgN1 <- mean(inv.N1[])
  avgN2 <- mean(inv.N2[])
  # calculating expected within-study variance
  sigmasqA <- (exp(tausq[1,1]/2-thetaA)+exp(tausq[1,1]/2+thetaA)+2)*avgN1
  sigmasqB <- (exp(tausq[2,2]/2-thetaB)+exp(tausq[2,2]/2+thetaB)+2)*avgN2
  
  # univariate I^2
  I2E[1] <- tausq[1,1]/(tausq[1,1]+sigmasqA)
  I2E[2] <- tausq[2,2]/(tausq[2,2]+sigmasqB)
  
  # bivariate I^2
  I2E.bivariate <-
  sqrt(exp(logdet(tausq[,])))/(sqrt(exp(logdet(tausq[,])))+sqrt(sigmasqA*sigmasqB))
  
  }"
  
  modelo <- textConnection(modelo.string)
  
  data2 <- list(K=nrow(X), X1=X$tp, X2=X$tn, N1=X$n1, N2=X$n0)
  
  ## Parametros iniciales
  
  set.seed(125)
  ##ni <- 100000; na <- 50000; nt <- 100
  ni <- 10000; na <- 500; nt <- 5
  
  inits <- function(){list(I2E.bivariate = runif(0,10), theta2 = rnorm(0,1), 
                           thetaA = rnorm(0,1), precision.theta2 = rnorm(0,1),
                           precision.thetaA = rnorm(0,1), lambda0 = rnorm(0,1), lambda1 = rnorm(0,1))}
  
  parameters <- c("I2E.bivariate", "I2E")
  
  modelo_bivariante <- jags.model(file=modelo, data = data2, 
                                  inits=inits, n.adapt=na, n.chains = 2, quiet = FALSE)
  
  modelo_bivariante_results <- coda.samples(modelo_bivariante, variable.names = parameters, n.iter = ni, thin=nt)
  
  
  
  # Intervalos modelo bivariante:
  inter_sen <- c(round(mean(modelo_bivariante_results[[1]][,c("I2E[1]")]) - sd(modelo_bivariante_results[[1]][,c("I2E[1]")]),2), 
                 round(mean(modelo_bivariante_results[[1]][,c("I2E[1]")]),2),
                 round(mean(modelo_bivariante_results[[1]][,c("I2E[1]")]) + sd(modelo_bivariante_results[[1]][,c("I2E[1]")]),2))
  
  inter_spe <- c(round(mean(modelo_bivariante_results[[1]][,c("I2E[2]")]) - sd(modelo_bivariante_results[[1]][,c("I2E[2]")]),2), 
                 round(mean(modelo_bivariante_results[[1]][,c("I2E[2]")]),2),
                 round(mean(modelo_bivariante_results[[1]][,c("I2E[2]")]) + sd(modelo_bivariante_results[[1]][,c("I2E[2]")]),2))
  
  inter_biv <- c(round(mean(modelo_bivariante_results[[1]][,c("I2E.bivariate")]) - sd(modelo_bivariante_results[[1]][,c("I2E.bivariate")]),2), 
                 round(mean(modelo_bivariante_results[[1]][,c("I2E.bivariate")]),2),
                 round(mean(modelo_bivariante_results[[1]][,c("I2E.bivariate")]) + sd(modelo_bivariante_results[[1]][,c("I2E.bivariate")]),2))
  
  
  tablaBivariante <- data.frame(
    Coefficient = c("I2E Sensitivity", "I2E Specificity", "I2E Bivariate"),
    Estimate = c(inter_sen[2], inter_spe[2], inter_biv[2]),
    LCI = c(inter_sen[1], inter_spe[1], inter_biv[1]),
    UPI = c(inter_sen[3], inter_spe[3], inter_biv[3])
  )
  colnames(tablaBivariante) <- c("Coefficient", "Estimate", "2.5 Percentile", "97.5 Percentile")
  
  return(tablaBivariante)
}

#heterogen <- heterogeneidad(X)






########
######## ANALYSIS OF SUBGROUPS
########



modelo_sg <- function(X3) {

  MA_Y = list()

  for (i in 1:length(X3)){

    X3[[i]] <- X3[[i]] %>% mutate(
      true1 = tp,
      true0 = tn,
      study = 1:nrow(X3[[i]])

    )

    X3[[i]] <- reshape(X3[[i]], direction = "long", varying = list(c("n1", "n0"), c("true1", "true0")), timevar = "sens", times = c(1,0),
                       v.names = c("n", "true"))

    X3[[i]]$spec <-  1-X3[[i]]$sens

    #names(MA_Y) = names(X3)

    MA_Y[[i]] = list(resultado = glmer(formula = cbind(true, n - true) ~ 0 + sens + spec + (0 + sens + spec|study), data = X3[[i]], family = binomial,
                                       nAGQ = 1, verbose = 0)
    )
    MA_Y[[i]]$summary = summary(MA_Y[[i]]$resultado)


    MA_Y[[i]]$tabla <- list(
      lsens = MA_Y[[i]]$summary$coeff[1,1],
      lspec = MA_Y[[i]]$summary$coeff[2,1],
      se.lsens = MA_Y[[i]]$summary$coeff[1,2],
      se.lspec = MA_Y[[i]]$summary$coeff[2,2],
      Sens = c(MA_Y[[i]]$summary$coeff[1,1],
               MA_Y[[i]]$summary$coeff[1,1]-qnorm(0.975)*MA_Y[[i]]$summary$coeff[1,2],
               MA_Y[[i]]$summary$coeff[1,1]+qnorm(0.975)*MA_Y[[i]]$summary$coeff[1,2]),

      Spec = c(MA_Y[[i]]$summary$coeff[2,1],
               MA_Y[[i]]$summary$coeff[2,1]-qnorm(0.975)*MA_Y[[i]]$summary$coeff[2,2],
               MA_Y[[i]]$summary$coeff[2,1]+qnorm(0.975)*MA_Y[[i]]$summary$coeff[2,2]),

      s = plogis(c(MA_Y[[i]]$summary$coeff[1,1],
                   MA_Y[[i]]$summary$coeff[1,1]-qnorm(0.975)*MA_Y[[i]]$summary$coeff[1,2],
                   MA_Y[[i]]$summary$coeff[1,1]+qnorm(0.975)*MA_Y[[i]]$summary$coeff[1,2])),
      e = plogis(c(MA_Y[[i]]$summary$coeff[2,1],
                   MA_Y[[i]]$summary$coeff[2,1]-qnorm(0.975)*MA_Y[[i]]$summary$coeff[2,2],
                   MA_Y[[i]]$summary$coeff[2,1]+qnorm(0.975)*MA_Y[[i]]$summary$coeff[2,2])),
      DOR = exp(MA_Y[[i]]$summary$coeff[1,1]+MA_Y[[i]]$summary$coeff[2,1] ),
      LRp = plogis(MA_Y[[i]]$summary$coeff[1,1])/(1-plogis(MA_Y[[i]]$summary$coeff[2,1])),
      LRm = ((1-plogis(MA_Y[[i]]$summary$coeff[1,1]))/plogis(MA_Y[[i]]$summary$coeff[2,1])),

      se.DOR = deltamethod (~ exp(x1+x2) , mean = c(MA_Y[[i]]$summary$coeff[1,1],MA_Y[[i]]$summary$coeff[2,1]) , cov = MA_Y[[i]]$summary$vcov ),

      se.LRp = deltamethod (~ (exp(x1)/(1+exp(x1)))/(1-(exp(x2)/(1+exp(x2)))) ,
                            mean = c(MA_Y[[i]]$summary$coeff[1,1],MA_Y[[i]]$summary$coeff[2,1]) , cov = MA_Y[[i]]$summary$vcov ),

      se.LRm = deltamethod (~ (1-(exp(x1)/(1+exp(x1))))/(exp(x2)/(1+exp(x2))) ,
                            mean = c(MA_Y[[i]]$summary$coeff[1,1],MA_Y[[i]]$summary$coeff[2,1]) , cov = MA_Y[[i]]$summary$vcov),


      phi = (seq(0, 2*pi, len=10000)),
      u = MA_Y[[i]]$summary$coeff[1,1],      ## E(logit(sen))
      v = (MA_Y[[i]]$summary$coeff[2,1]),          ## E(logit(spe)
      SELogitSen = MA_Y[[i]]$summary$coeff[1,2], ## standard error E(logit(sen))
      SELogitSpe = MA_Y[[i]]$summary$coeff[2,2],  ## standard error E(logit(spe))
      covAB = (MA_Y[[i]]$summary$vcov[1,2]),     ##Covariance betwen E(logit(sen)) and E(logit(spe))

      VarLogitSen =  MA_Y[[i]]$summary$varcor$study[1,1],   ##Var(logit(sen))
      VarLogitSpe =  MA_Y[[i]]$summary$varcor$study[2,2], ##Var(logit(spe))
      corAB =  MA_Y[[i]]$summary$varcor$study[2],
      N = MA_Y[[i]]$summary$ngrps,

      SEPA = (sqrt(MA_Y[[i]]$summary$varcor$study[1,1]+MA_Y[[i]]$summary$coeff[1,2]^2)),
      SEPB = (sqrt(MA_Y[[i]]$summary$varcor$study[2,2]+MA_Y[[i]]$summary$coeff[2,2]^2)),
      sAB = ((MA_Y[[i]]$summary$varcor$study[2])*(sqrt(MA_Y[[i]]$summary$varcor$study[1,1]))*(sqrt(MA_Y[[i]]$summary$varcor$study[2,2]))),
      rho = ((((MA_Y[[i]]$summary$varcor$study[2])*(sqrt(MA_Y[[i]]$summary$varcor$study[1,1]))*(sqrt(MA_Y[[i]]$summary$varcor$study[2,2])))+(MA_Y[[i]]$summary$vcov[1,2]))/((sqrt(MA_Y[[i]]$summary$varcor$study[1,1]+MA_Y[[i]]$summary$coeff[1,2]^2))*(sqrt(MA_Y[[i]]$summary$varcor$study[2,2]+MA_Y[[i]]$summary$coeff[2,2]^2)))),


      croot = (sqrt(2*qf(0.95,2,MA_Y[[i]]$summary$ngrps-2))),

      a=(sqrt(MA_Y[[i]]$summary$varcor$study[1,1]+MA_Y[[i]]$summary$coeff[1,2]^2))*(sqrt(2*qf(0.95,2,MA_Y[[i]]$summary$ngrps-2))),
      b=(sqrt(MA_Y[[i]]$summary$varcor$study[2,2]+MA_Y[[i]]$summary$coeff[2,2]^2))*(sqrt(2*qf(0.95,2,MA_Y[[i]]$summary$ngrps-2))),
      sen = exp(MA_Y[[i]]$summary$coeff[1,1])/(1+exp(MA_Y[[i]]$summary$coeff[1,1])),
      esp = exp(MA_Y[[i]]$summary$coeff[2,1])/(1+exp(MA_Y[[i]]$summary$coeff[2,1])),

      y = MA_Y[[i]]$summary$coeff[1,1] + (sqrt(MA_Y[[i]]$summary$varcor$study[1,1]+MA_Y[[i]]$summary$coeff[1,2]^2))*(sqrt(2*qf(0.95,2,MA_Y[[i]]$summary$ngrps-2))) * cos((seq(0, 2*pi, len=10000)) + acos(((((MA_Y[[i]]$summary$varcor$study[2])*(sqrt(MA_Y[[i]]$summary$varcor$study[1,1]))*(sqrt(MA_Y[[i]]$summary$varcor$study[2,2])))+(MA_Y[[i]]$summary$vcov[1,2]))/((sqrt(MA_Y[[i]]$summary$varcor$study[1,1]+MA_Y[[i]]$summary$coeff[1,2]^2))*(sqrt(MA_Y[[i]]$summary$varcor$study[2,2]+MA_Y[[i]]$summary$coeff[2,2]^2)))))),    ## Sensitivity (logit scale)
      x = (MA_Y[[i]]$summary$coeff[2,1]) + (sqrt(MA_Y[[i]]$summary$varcor$study[2,2]+MA_Y[[i]]$summary$coeff[2,2]^2))*(sqrt(2*qf(0.95,2,MA_Y[[i]]$summary$ngrps-2))) * cos((seq(0, 2*pi, len=10000)))          ,       ## Specificity (logit scale)
      specificity = (1-(exp((MA_Y[[i]]$summary$coeff[2,1]) + (sqrt(MA_Y[[i]]$summary$varcor$study[2,2]+MA_Y[[i]]$summary$coeff[2,2]^2))*(sqrt(2*qf(0.95,2,MA_Y[[i]]$summary$ngrps-2))) * cos((seq(0, 2*pi, len=10000))))/(1+exp((MA_Y[[i]]$summary$coeff[2,1]) + (sqrt(MA_Y[[i]]$summary$varcor$study[2,2]+MA_Y[[i]]$summary$coeff[2,2]^2))*(sqrt(2*qf(0.95,2,MA_Y[[i]]$summary$ngrps-2))) * cos((seq(0, 2*pi, len=10000))))))) , ## ROC scale
      sensitivity = (exp( MA_Y[[i]]$summary$coeff[1,1] + (sqrt(MA_Y[[i]]$summary$varcor$study[1,1]+MA_Y[[i]]$summary$coeff[1,2]^2))*(sqrt(2*qf(0.95,2,MA_Y[[i]]$summary$ngrps-2))) * cos((seq(0, 2*pi, len=10000)) + acos(((((MA_Y[[i]]$summary$varcor$study[2])*(sqrt(MA_Y[[i]]$summary$varcor$study[1,1]))*(sqrt(MA_Y[[i]]$summary$varcor$study[2,2])))+(MA_Y[[i]]$summary$vcov[1,2]))/((sqrt(MA_Y[[i]]$summary$varcor$study[1,1]+MA_Y[[i]]$summary$coeff[1,2]^2))*(sqrt(MA_Y[[i]]$summary$varcor$study[2,2]+MA_Y[[i]]$summary$coeff[2,2]^2)))))))/(1+exp( MA_Y[[i]]$summary$coeff[1,1] + (sqrt(MA_Y[[i]]$summary$varcor$study[1,1]+MA_Y[[i]]$summary$coeff[1,2]^2))*(sqrt(2*qf(0.95,2,MA_Y[[i]]$summary$ngrps-2))) * cos((seq(0, 2*pi, len=10000)) + acos(((((MA_Y[[i]]$summary$varcor$study[2])*(sqrt(MA_Y[[i]]$summary$varcor$study[1,1]))*(sqrt(MA_Y[[i]]$summary$varcor$study[2,2])))+(MA_Y[[i]]$summary$vcov[1,2]))/((sqrt(MA_Y[[i]]$summary$varcor$study[1,1]+MA_Y[[i]]$summary$coeff[1,2]^2))*(sqrt(MA_Y[[i]]$summary$varcor$study[2,2]+MA_Y[[i]]$summary$coeff[2,2]^2)))))))))  ,    ## ROC scale


      area = polyarea((1-(exp((MA_Y[[i]]$summary$coeff[2,1]) + (sqrt(MA_Y[[i]]$summary$varcor$study[2,2]+MA_Y[[i]]$summary$coeff[2,2]^2))*(sqrt(2*qf(0.95,2,MA_Y[[i]]$summary$ngrps-2))) * cos((seq(0, 2*pi, len=10000))))/(1+exp((MA_Y[[i]]$summary$coeff[2,1]) + (sqrt(MA_Y[[i]]$summary$varcor$study[2,2]+MA_Y[[i]]$summary$coeff[2,2]^2))*(sqrt(2*qf(0.95,2,MA_Y[[i]]$summary$ngrps-2))) * cos((seq(0, 2*pi, len=10000))))))),
                      (exp( MA_Y[[i]]$summary$coeff[1,1] + (sqrt(MA_Y[[i]]$summary$varcor$study[1,1]+MA_Y[[i]]$summary$coeff[1,2]^2))*(sqrt(2*qf(0.95,2,MA_Y[[i]]$summary$ngrps-2))) * cos((seq(0, 2*pi, len=10000)) + acos(((((MA_Y[[i]]$summary$varcor$study[2])*(sqrt(MA_Y[[i]]$summary$varcor$study[1,1]))*(sqrt(MA_Y[[i]]$summary$varcor$study[2,2])))+(MA_Y[[i]]$summary$vcov[1,2]))/((sqrt(MA_Y[[i]]$summary$varcor$study[1,1]+MA_Y[[i]]$summary$coeff[1,2]^2))*(sqrt(MA_Y[[i]]$summary$varcor$study[2,2]+MA_Y[[i]]$summary$coeff[2,2]^2)))))))/(1+exp( MA_Y[[i]]$summary$coeff[1,1] + (sqrt(MA_Y[[i]]$summary$varcor$study[1,1]+MA_Y[[i]]$summary$coeff[1,2]^2))*(sqrt(2*qf(0.95,2,MA_Y[[i]]$summary$ngrps-2))) * cos((seq(0, 2*pi, len=10000)) + acos(((((MA_Y[[i]]$summary$varcor$study[2])*(sqrt(MA_Y[[i]]$summary$varcor$study[1,1]))*(sqrt(MA_Y[[i]]$summary$varcor$study[2,2])))+(MA_Y[[i]]$summary$vcov[1,2]))/((sqrt(MA_Y[[i]]$summary$varcor$study[1,1]+MA_Y[[i]]$summary$coeff[1,2]^2))*(sqrt(MA_Y[[i]]$summary$varcor$study[2,2]+MA_Y[[i]]$summary$coeff[2,2]^2))))))))))

    )


    MA_Y[[i]]$tabla2 <- list(
      tabla2 = data.frame(
        sensitivity = MA_Y[[i]]$tabla$sensitivity, specificity = MA_Y[[i]]$tabla$specificity)
    )


    MA_Y[[i]]$resultados_final <- list(
      df = data.frame(
        Estimate = c(round(MA_Y[[i]]$tabla$s[1],3), round(MA_Y[[i]]$tabla$e[1],3), round(MA_Y[[i]]$tabla$DOR,3) , round(MA_Y[[i]]$tabla$LRp,3) , round(MA_Y[[i]]$tabla$LRm,3)) ,
        LCI = c(round(MA_Y[[i]]$tabla$s[2],3), round(MA_Y[[i]]$tabla$e[2],3), round(MA_Y[[i]]$tabla$DOR-qnorm(0.975)*MA_Y[[i]]$tabla$se.DOR,3) , round(MA_Y[[i]]$tabla$LRp-qnorm(0.975)*MA_Y[[i]]$tabla$se.LRp,3) , round(MA_Y[[i]]$tabla$LRm-qnorm(0.975)*MA_Y[[i]]$tabla$se.LRm,3)) ,
        UCI = c(round(MA_Y[[i]]$tabla$s[3],3), round(MA_Y[[i]]$tabla$e[3],3), round(MA_Y[[i]]$tabla$DOR+qnorm(0.975)*MA_Y[[i]]$tabla$se.DOR,3) , round(MA_Y[[i]]$tabla$LRp+qnorm(0.975)*MA_Y[[i]]$tabla$se.LRp,3) , round(MA_Y[[i]]$tabla$LRm+qnorm(0.975)*MA_Y[[i]]$tabla$se.LRm,3)),
        row.names = c("Sensitivity","Specificity", "DOR", "LR+" , "LR-" )
      )
    )
  }

  return(MA_Y = MA_Y)
}

#modelo_covariables <-  modelo_sg(X3)




########
########  META-REGRESSION TABLE
######## 

mr_table <- function(model_A, model_final){
  
  # summary of model A:
  ma_biv_A = summary(model_A)
  
  # logit sens and logit spec
  lsens_A = ma_biv_A$coeff[1,1]
  lspec_A = ma_biv_A$coeff[2,1]
  
  
  se.lsens_A = ma_biv_A$coeff[1,2]
  se.lspec_A = ma_biv_A$coeff[2,2]
  
  # 95% confidence intervals for logit sens and logit spec
  Sens_A = c(lsens_A, lsens_A-qnorm(0.975)*se.lsens_A, lsens_A+qnorm(0.975)*se.lsens_A)
  Spec_A = c(lspec_A, lspec_A-qnorm(0.975)*se.lspec_A, lspec_A+qnorm(0.975)*se.lspec_A)
  
  # sens and spec estimates in the raw scale
  s_A <-plogis( Sens_A ) 
  e_A <-plogis( Spec_A ) 
  
  
  
  
  # summary of the final model:
  ma_biv = summary(model_final)
  
  
  if (nrow(ma_biv$coeff) == 2) {
    # logit sens and logit spec
    lsens = ma_biv$coeff[1,1]
    lspec = ma_biv$coeff[2,1]
    
    
    se.lsens = ma_biv$coeff[1,2]
    se.lspec = ma_biv$coeff[2,2]
    
    # 95% confidence intervals for logit sens and logit spec
    Sens = c(lsens, lsens-qnorm(0.975)*se.lsens, lsens+qnorm(0.975)*se.lsens)
    Spec = c(lspec, lspec-qnorm(0.975)*se.lspec, lspec+qnorm(0.975)*se.lspec)
    
    # sens and spec estimates in the raw scale
    s <-plogis( Sens ) 
    e <-plogis( Spec ) 
    
    
    
    tablaOutPut = data.frame(
      Coefficient = c("Sensitivity", "Specificity"),
      Estimate_A = c(s_A[1], e_A[1]),
      LCI_A = c(
        round(s_A[2], 3),
        round(e_A[2], 3)
        
      ),
      UCI_A = c(
        round(s_A[3], 3),
        round(e_A[3], 3)
      ),
      Estimate_final = c(s[1], e[1]) ,
      LCI_final = c(
        round(s[2], 3),
        round(e[2], 3)
        
      ),
      UCI_final = c(
        round(s[3], 3),
        round(e[3], 3)
      )
    )
  } else if (nrow(ma_biv$coeff) == 4){
    # coeficientes logit
    lsens_0 = ma_biv$coefficients[1,1]
    lsens_1 = ma_biv$coefficients[2,1]
    lspec_0 = ma_biv$coefficients[3,1]
    lspec_1 = ma_biv$coefficients[4,1]
    
    # errores estandar logit
    se_lsens_0 = ma_biv$coefficients[1,2]
    se_lsens_1 = ma_biv$coefficients[2,2]
    se_lspec_0 = ma_biv$coefficients[3,2]
    se_lspec_1 = ma_biv$coefficients[4,2]
    
    # 95% confidence intervals for logit sens and logit spec
    Sens_0 = c(lsens_0, lsens_0-qnorm(0.975)*se_lsens_0, lsens_0+qnorm(0.975)*se_lsens_0)
    Sens_1 = c(lsens_1, lsens_1-qnorm(0.975)*se_lsens_1, lsens_1+qnorm(0.975)*se_lsens_1)
    Spec_0 = c(lspec_0, lspec_0-qnorm(0.975)*se_lspec_0, lspec_0+qnorm(0.975)*se_lspec_0)
    Spec_1 = c(lspec_1, lspec_1-qnorm(0.975)*se_lspec_1, lspec_1+qnorm(0.975)*se_lspec_1)
    
    # sens and spec estimates in the raw scale
    s_0 <-plogis( Sens_0 ) 
    s_1  <-plogis( Sens_1 ) 
    e_0 <-plogis( Spec_0 ) 
    e_1 <-plogis( Spec_1 ) 
    
    
    tablaOutPut = data.frame(
      Coefficient = c("Sensitivity", "Specificity"),
      Estimate_A = c(s_A[1], e_A[1]),
      LCI_A = c(
        round(s_A[2], 3),
        round(e_A[2], 3)
        
      ),
      UCI_A = c(
        round(s_A[3], 3),
        round(e_A[3], 3)
      ),
      Estimate_final_0 = c(s_0[1], e_0[1]) ,
      LCI_final_0 = c(
        round(s_0[2], 3),
        round(e_0[2], 3)
        
      ),
      UCI_final_0 = c(
        round(s_0[3], 3),
        round(e_0[3], 3)
      ),
      Estimate_final_1 = c(s_1[1], e_1[1]) ,
      LCI_final_1 = c(
        round(s_1[2], 3),
        round(e_1[2], 3)
        
      ),
      UCI_final_1 = c(
        round(s_1[3], 3),
        round(e_1[3], 3)
      )
    )
    
  } else {
    print("NOT YET CONFIGURED")
  }
  
  
  
  return(list(tabla=tablaOutPut))
}

#mr_table(selected)$tabla



