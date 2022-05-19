library(deSolve)
library(ggplot2)
library(dplyr)
library(ggpmisc)
library(patchwork)

# IFFL model
model_IFFL = function (t, y, parms){
  with(as.list(c(y,parms)),{
  dCRP = beta_CRP*(AMP^n)/(k_AMP^n+AMP^n)-gamma_CRP*CRP
  dgalS = beta_galS*(CRP^n)/(k_CRP^n+CRP^n)-gamma_galS*galS
  dgalE = beta_galE*(k_galS^n)*(CRP^n)/(k_galS^n+galS^n)*(k_CRP2^n+CRP^n)-gamma_galE*galE
  list(c(dCRP, dgalS, dgalE))})
  }
# FFL model
model_FFL = function (t, y, parms){
  with(as.list(c(y,parms)),{
    dCRP = beta_CRP*(AMP^n)/(k_AMP^n+AMP^n)-gamma_CRP*CRP
    daraC = beta_araC*(CRP^n)/(k_CRP^n+CRP^n)-gamma_araC*araC
    dBAD = beta_BAD*(araC^n)*(CRP^n)/(k_araC^n+araC^n)*(k_CRP2^n+CRP^n)-gamma_BAD*BAD
    list(c(dCRP, daraC, dBAD))})
}
# function
solveode_FFL<-function(yini_FFL,times_FFL,model_FFL,parms){
  for (index in 1:length(parms)) {
    if(length(parms[[index]])!=1){
      target_index=index
    }
  }
  bad_vals<-c()
  for(parm in parms[[target_index]]){
    parmstmp<-parms
    parmstmp[[target_index]]<-parm
    parms_new<-unlist(parmstmp)
    out_FFL<-ode(y = yini_FFL,times = times_FFL, func = model_FFL, parms = parms_new)
    outFFL_TIME20<-data.frame(out_FFL)%>%filter(time==20)
    bad_val<-outFFL_TIME20$BAD
    bad_vals<-c(bad_vals,bad_val)
  }
  result<-cbind(bad_vals, parms[[target_index]])
  return(result)
}
solveode_IFFL<-function(yini_IFFL,times_IFFL,model_IFFL,parms){
  for (index in 1:length(parms)) {
    if(length(parms[[index]])!=1){
      target_index=index
    }
  }
  galE_vals<-c()
  for(parm in parms[[target_index]]){
    parmstmp<-parms
    parmstmp[[target_index]]<-parm
    parms_new<-unlist(parmstmp)
    out_IFFL<-ode(y = yini_IFFL,times = times_IFFL, func = model_IFFL, parms = parms_new)
    outIFFL_TIME20<-data.frame(out_IFFL)%>%filter(time==20)
    galE_val<-outIFFL_TIME20$galE
    galE_vals<-c(galE_vals,galE_val)
  }
  result<-cbind(galE_vals, parms[[target_index]])
  return(result)
}

# set up parameters for IFFL 
yini_IFFL = c(CRP=0, galS=0, galE=0)
parms_IFFL <- c(AMP=10, beta_CRP=0.5, beta_galS=0.5, beta_galE=1, 
                n=1, k_AMP=0.1, k_CRP=0.2, k_galS=0.3, k_CRP2=0, 
                gamma_CRP=1, gamma_galS=1, gamma_galE=0.5)
times_IFFL <- seq(from = 0, to = 25, by = 0.01)
out_IFFL <- ode(y = yini_IFFL,times = times_IFFL, func = model_IFFL, parms = parms_IFFL)
matplot.deSolve(out_IFFL)


# set up parameters for FFL
yini_FFL = c(CRP=0, araC=0, BAD=0)
parms_FFL <- c(AMP=10, beta_CRP=0.5, beta_araC=0.5, beta_BAD=1, 
               n=1, k_AMP=0.1, k_CRP=0.2, k_araC=0.3, k_CRP2=0, 
               gamma_CRP=1, gamma_araC=1, gamma_BAD=0.5)
times_FFL <- seq(from = 0, to = 25, by = 0.01)
out_FFL <- ode(y = yini_FFL,times = times_FFL, func = model_FFL, parms = parms_FFL)
matplot.deSolve(out_FFL)

# make the matrix
parms_FFL <- list(AMP=10, beta_CRP=0.5, beta_araC=0.5, beta_BAD=1, 
                  n=1, k_AMP=0.1, k_CRP=0.2, k_araC=0.3, k_CRP2=seq(0,1,0.1), 
                  gamma_CRP=1, gamma_araC=1, gamma_BAD=0.5)
result_k_CRP2<-solveode_FFL(yini_FFL,times_FFL,model_FFL,parms_FFL)
colnames(result_k_CRP2)<- c("BAD_vals","k_CRP2")

parms_IFFL <- list(AMP=10, beta_CRP=0.5, beta_galS=0.5, beta_galE=1, 
                   n=1, k_AMP=0.1, k_CRP=0.2, k_galS=0.3, k_CRP2=seq(0,1,0.1), 
                   gamma_CRP=1, gamma_galS=1, gamma_galE=0.5)
result_k_CRP<-solveode_IFFL(yini_IFFL,times_IFFL,model_IFFL,parms_IFFL)
colnames(result_k_CRP)<- c("galE_vals","k_CRP2")

parms_IFFL <- list(AMP=10, beta_CRP=0.5, beta_galS=0.2, beta_galE=1, 
                   n=1, k_AMP=0.1, k_CRP=0.2, k_galS=0.3, k_CRP2=0.1, 
                   gamma_CRP=1, gamma_galS=1, gamma_galE=0.5)
result_gamma_galS<-solveode_IFFL(yini_IFFL,times_IFFL,model_IFFL,parms_IFFL)
colnames(result_gamma_galS)<- c("galE_vals","gamma_galS")

# plot
p1 <- ggplot(as.data.frame(result_k_CRP2),aes(k_CRP2,BAD_vals))+
  geom_point()+
  stat_smooth(method = lm,level=0.99)+
  stat_poly_eq(aes(label =paste(..eq.label.., ..adj.rr.label.., sep = '~~')),
               formula = y ~ x,  parse = TRUE,
               family="serif",
               size = 3.6,
               color="black",
               label.x = 0.1,
               label.y = 1)


p2 <- ggplot(as.data.frame(result_k_CRP),aes(k_CRP2,galE_vals))+
  geom_point()+
  stat_smooth(method = lm,level=0.99)+
  stat_poly_eq(aes(label =paste(..eq.label.., ..adj.rr.label.., sep = '~~')),
               formula = y ~ x,  parse = TRUE,
               family="serif",
               size = 3.6,
               color="black",
               label.x = 0.1,
               label.y = 1)

ggplot(as.data.frame(result_gamma_galS),aes(gamma_galS,galE_vals))+
  geom_point()+
  stat_smooth(method = lm,level=0.99)+
  stat_poly_eq(aes(label =paste(..eq.label.., ..adj.rr.label.., sep = '~~')),
               formula = y ~ x,  parse = TRUE,
               family="serif",
               size = 3.6,
               color="black",
               label.x = 0.1,
               label.y = 1)

data2<-as.data.frame(data2)
data2[,1] <- as.numeric(data2[,1])
data2[,2] <- as.numeric(data2[,2])
colnames(data2)<- c("relative_concentration","k_CRP2","group")
ggplot(data2,aes(k_CRP2,relative_concentration))+
  geom_point(aes(col = group))

