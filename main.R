## This is the main wrapper function. It performs all analyses specified
main = function(age,depth,proxy, events=NULL,nsims=10000, eventmeasure = "depth",reference.label=NULL,proxy.type="d18O",
                transform="identity",reg.model = list(
  const=FALSE,depth1=FALSE,depth2=TRUE,proxy=TRUE,psi0=TRUE,psi1=TRUE), noise="ar1", method="inla", CI.type="quantiles",
  event.estimation = NULL, bias = NULL,store.everything=FALSE,print.progress=FALSE){

  time.start = Sys.time()
  main.call = sys.call(which=1)
  if(print.progress) cat("Initiating data formatting...",sep="")

  #prepare object by formatting dataset, writing formulastrings etc
  object = prepare(age,depth,proxy, events=events,nsims=nsims, eventmeasure = eventmeasure,reg.model = reg.model,
                         noise=noise, method=method,reference.label=reference.label,transform=transform,proxy.type=proxy.type)
  if(print.progress) cat(" completed!\n",sep="")
  #fit the data, first by least squares, then by INLA (if specified)
  object = modelfitter(object, method=method,print.progress=print.progress)

  #produce samples from the chronologies
  object = chronology_simulation(object, nsims=nsims, method=method,store.means=store.everything,print.progress=print.progress)
  
  #compute posterior marginal mean, quantiles and other summary statistics
  object = simulationsummarizer(object,CI.type=CI.type,print.progress=print.progress)
  
  #if event.estimation list object (containing specifications) is included, perform dating estimation
  if(!is.null(event.estimation)){
    #find onset depth posterior by fitting linear ramp model with INLA
    object = linrampfitter(object,interval=event.estimation$interval,h=event.estimation$h,t1.sims=event.estimation$t1.sims,
               rampsims=event.estimation$rampsims,label=event.estimation$label,
               depth.reference = event.estimation$depth.reference,
               print.progress=print.progress,log.ramp=event.estimation$log.ramp)

    #perform Monte Carlo simulations to produce samples for onset age of warming transition
    object = event_depth_to_age(object, nsims = nsims, print.progress=print.progress,label=event.estimation$label,age.reference = event.estimation$age.reference)
  }
  #if bias list object is included, perform this analysis
  if(!is.null(bias)){
    object = biased_chronologies(object,bias.model=bias$bias.model,biasparams = bias$biasparams,nsims=nsims,store.samples=store.everything)
  }
  time.total = difftime(Sys.time(), time.start,units="secs")[[1]]
  object$.args$call = main.call
  object$time$total = time.total
  
  if(print.progress){
    cat("Completed everything in ",time.total, "seconds\n",sep="")
  }
  return(object)
}

#remember to set proper directory
setwd("/Submission files/Code/")

#load functions, data and packages
source("biased_chronologies.R"); source("chronology_simulation.R"); source("event_depth_to_age.R"); source("helpfulfunctions.R"); source("linrampfitter.R")
source("modelfitter.R"); source("plot_results.R"); source("prepare.R"); source("rgeneric_model.R"); source("summary_results.R"); source("rgeneric.heteroskedastic.R")
library("INLA"); library("matrixStats"); library("numDeriv")

library(readODS)
library(readxl)
#
maindata = read_excel("/Submission files/Code/datasets_used/NGRIP_d18O_and_dust_5cm.xls",
                      sheet="NGRIP-2 d18O and Dust",col_types="numeric" )



startindex = 2921 #11703yb2k (holocene onset)

depth = maindata$`NGRIP-2 depth (m)`[startindex:nrow(maindata)]
age = maindata$`GICC05 age (yr b2k)`[startindex:nrow(maindata)]
MCE = maindata$`GICC05 MCE (yr)`[startindex:nrow(maindata)]
water = maindata$`Delta O18 (permil)`[startindex:nrow(maindata)]#dust, NA until 1167
dust0 = maindata$`Dust count (ml^-1)`[startindex:nrow(maindata)]#dust, NA until 1167

library(zoo)
dust = na.approx(dust0) #Impute missing data using linear interpolation

dust = log(dust)

do.dust = FALSE #set FALSE for d18O proxy, and true for Ca2+ proxy

if(do.dust){
  proxy.type = "ca"
  proxy=dust
  
}else{
  proxy.type = "d18O"
  proxy=water
}


#load abrupt warming transitions
eventdata = read_ods("datasets_used/GISevents.ods")
GISevents = eventdata[(eventdata$`NGRIP depth (m)`>min(depth))&(eventdata$`NGRIP depth (m)`<max(depth)),]
event_intervals = read.table("datasets_used/event_intervals.txt")


eventdepths = GISevents$`NGRIP depth (m)`
eventindexes = c(1,which.index(eventdepths, depth[2:length(depth)]) )


#plot d18O proxies as a function of depth and age (GICC05), respectively
par(mfrow=c(1,1),mar=c(5,4.5,2,2)+0.1)
plot(depth,water,type="l",xlab="Depth (m)",ylab=expression(paste(delta^18,"O (permil)")),xlim=rev(range(depth))); abline(v=eventdepths,lwd=0.7,col="gray")
plot(age,water,type="l",xlab="Age (yb2k)",ylab=expression(paste(delta^18,"O (permil)")),xlim=rev(range(age))); abline(v=age[eventindexes],lwd=0.7,col="gray")


nsims=10000

#events used in regression model
events=eventdepths #locations of transitions used in regression model. Pairs with 'eventmeasure' for finding the corresponding indices
eventmeasure="depth"
#reg.model specifies the structure for which the formula string should be created
reg.model = list(
  const=FALSE,depth1=FALSE,depth2=TRUE,proxy=TRUE,psi0=TRUE,psi1=TRUE); method="inla";
plots = list(posteriors=TRUE) 


eventnumber=13 #number between 1 and 29. specifies which transition to consider

#load data window and specifics to transition
lowerints = which.index(event_intervals$depth_int_lower.m, depth[2:length(depth)])
upperints = which.index(event_intervals$depth_int_upper.m, depth[2:length(depth)])
depth.reference = event_intervals$NGRIP_depth_m[eventnumber]
age.reference = event_intervals$GICC_age.yb2k[eventnumber]
interval = lowerints[eventnumber]:upperints[eventnumber]



#list object specifying parameters for abrupt warming transition dating
event.estimation = list(interval=interval,t1.sims=50000,rampsims=50000,label="GI-11",
                     depth.reference=event_intervals$NGRIP_depth_m[eventnumber],
                     age.reference=event_intervals$NGRIP_agee_m[eventnumber])

reference.label="GICC05" #this string is used in plot_results

method="inla" #should INLA be run in fitting procedure? alternatively "LS" for no (not quite as tested)
CI.type="quantiles" #method used for credible intervals. Alternatives are 'quantiles' and 'hpd'. For Gaussian processes (or any unimodal and symmetric) these are equal


print.progress=TRUE

transform = "identity" #set equal to "log" for logarithmic transformation
noise = "ar1" #iid, ar1 and ar2 supported


bias = list(bias.model="uniform",biasparams=cbind( c(1,1),c(0.98,1.02),c(0.96,1.04) ),
            store.samples=FALSE)
store.everything=FALSE #should fixed model components and stochastic both be stored (TRUE), or just the sum (FALSE)

#run main function wrapper
object = main(age,depth,proxy, events=eventdepths,nsims=nsims, eventmeasure = eventmeasure,proxy.type=proxy.type,
              reference.label=reference.label,transform=transform,reg.model = reg.model, 
              noise=noise, method=method, CI.type=CI.type,
  event.estimation = event.estimation,store.everything=store.everything,
  print.progress=print.progress,bias = bias
  )

plot_results(object) #plot results


summary_results(object) #print out summary statistics

  
##############################################
### Comparing models: iid, AR(1) and AR(2) ###
##############################################

transform = "identity" #set equal to "log" for logarithmic transformation, "identity" for no transformation
do.dust = TRUE #TRUE for Ca2+ proxy, FALSE for d18O proxy
if(do.dust){
  proxy.type = "ca"
  proxy=dust
  
}else{
  proxy.type = "d18O"
  proxy=water
}


#perform analysis on iid, ar1 and ar2 models
objectiid= main(age,depth,proxy,eventdepths,transform=transform,noise="iid")
objectar1= main(age,depth,proxy,eventdepths,transform=transform,noise="ar1")
objectar2= main(age,depth,proxy,eventdepths,transform=transform,noise="ar2")


#summary results
summary_results(objectiid); summary_results(objectar1); summary_results(objectar2)


## plot uncertainties and compare
xlim = rev(range(objectiid$data$z))
gicc05=objectiid$data$y
ylim=range(objectar2$simulation$summary$lower-gicc05,objectar2$simulation$summary$upper-gicc05)

plot(x=objectiid$data$z,y=objectiid$simulation$summary$mean-gicc05,xlim=xlim,ylim=ylim,xlab="Depth (m)",ylab="Simulated time scale - GICC05 (years)",
     type="l",col="gray",lwd=1.5)
lines(x=objectiid$data$z,y=objectiid$simulation$summary$lower-gicc05,col="black",lwd=1.5)
lines(x=objectiid$data$z,y=objectiid$simulation$summary$upper-gicc05,col="black",lwd=1.5)
lines(x=objectiid$data$z,y=objectar1$simulation$summary$lower-gicc05,col="blue",lwd=1.5)
lines(x=objectiid$data$z,y=objectar1$simulation$summary$upper-gicc05,col="blue",lwd=1.5)
lines(x=objectiid$data$z,y=objectar2$simulation$summary$lower-gicc05,col="red",lwd=1.5)
lines(x=objectiid$data$z,y=objectar2$simulation$summary$upper-gicc05,col="red",lwd=1.5)
abline(h=0,lty=3,col="blue",lwd=0.8)


if(transform == "log"){
  legend(1600,310,legend=c("Mean","iid CI", "AR(1) CI", "AR(2) CI"),
         col=c("gray","black","blue","red"),
         lty=c(1,1), cex=0.5)
}else{
  legend(1570,250,legend=c("Mean","iid CI", "AR(1) CI", "AR(2) CI"),
         col=c("gray","black","blue","red"),
         lty=c(1,1), cex=0.4)
}



#plot posterior marginal distributions of model parameters for the different models
layout(mat=matrix(c(1,2,4,7,3,5,8,9,6),nrow=3))
par(mar=c(4,4,2,1))
xrange_sigma = c(0.420,0.442)
xrange_sigma = range(objectiid$fitting$hyperparameters$posteriors$sigma_epsilon[,1],
                     objectar1$fitting$hyperparameters$posteriors$sigma_epsilon[,1],
                     objectar2$fitting$hyperparameters$posteriors$sigma_epsilon[,1])
xrange_phi = range(objectar1$fitting$hyperparameters$posteriors$phi[,1],
                   objectar2$fitting$hyperparameters$posteriors$phi1[,1],
                   objectar2$fitting$hyperparameters$posteriors$phi2[,1])

plot(objectiid$fitting$hyperparameters$posteriors$sigma_epsilon,xlab=expression(paste(sigma[epsilon])),ylab="Density",lwd=2,type="l",main="(a) Independent residuals",xlim=xrange_sigma)
abline(v=objectiid$fitting$hyperparameters$results$sigma_epsilon$mean);abline(v=c(objectiid$fitting$hyperparameters$results$sigma_epsilon$quant0.025,objectiid$fitting$hyperparameters$results$sigma_epsilon$quant0.975),col="gray")

plot(objectar1$fitting$hyperparameters$posteriors$sigma_epsilon,xlab=expression(paste(sigma[epsilon])),ylab="Density",lwd=2,type="l",main="(b) AR(1) residuals",xlim=xrange_sigma)
abline(v=objectar1$fitting$hyperparameters$results$sigma_epsilon$mean);abline(v=c(objectar1$fitting$hyperparameters$results$sigma_epsilon$quant0.025,objectar1$fitting$hyperparameters$results$sigma_epsilon$quant0.975),col="gray")
plot(objectar1$fitting$hyperparameters$posteriors$phi,xlab=expression(phi),ylab="Density",lwd=2,type="l",main="(c) AR(1) residuals",xlim=xrange_phi)
abline(v=objectar1$fitting$hyperparameters$results$phi$mean);abline(v=c(objectar1$fitting$hyperparameters$results$phi$quant0.025,objectar1$fitting$hyperparameters$results$phi$quant0.975),col="gray")

plot(objectar2$fitting$hyperparameters$posteriors$sigma_epsilon,xlab=expression(paste(sigma[epsilon])),ylab="Density",lwd=2,type="l",main="(d) AR(2) residuals",xlim=xrange_sigma)
abline(v=objectar2$fitting$hyperparameters$results$sigma_epsilon$mean);abline(v=c(objectar2$fitting$hyperparameters$results$sigma_epsilon$quant0.025,objectar2$fitting$hyperparameters$results$sigma_epsilon$quant0.975),col="gray")
plot(objectar2$fitting$hyperparameters$posteriors$phi1,xlab=expression(phi[1]),ylab="Density",lwd=2,type="l",main="(e) AR(2) residuals",xlim=xrange_phi)
abline(v=objectar2$fitting$hyperparameters$results$phi1$mean);abline(v=c(objectar2$fitting$hyperparameters$results$phi1$quant0.025,objectar2$fitting$hyperparameters$results$phi1$quant0.975),col="gray")
plot(objectar2$fitting$hyperparameters$posteriors$phi2,xlab=expression(phi[2]),ylab="Density",lwd=2,type="l",main="(f) AR(2) residuals",xlim=xrange_phi)
abline(v=objectar2$fitting$hyperparameters$results$phi2$mean);abline(v=c(objectar2$fitting$hyperparameters$results$phi2$quant0.025,objectar2$fitting$hyperparameters$results$phi2$quant0.975),col="gray")

par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1)



#############################
### Abrupt warming events ###
#############################

transform = "identity" ##"log" for logaritmic transformation, "identity" for ordinary  
nsims = 10000

eventdata = read_ods("datasets_used/GISevents.ods")
event_intervals = read.table("datasets_used/event_intervals.txt")


#dust

  depth = maindata$`NGRIP-2 depth (m)`[1167:nrow(maindata)]
  proxy0 = maindata$`Dust count (ml^-1)`[1167:nrow(maindata)]#dust, NA until 1167
  age = maindata$`GICC05 age (yr b2k)`[1167:nrow(maindata)]
  MCE = maindata$`GICC05 MCE (yr)`[1167:nrow(maindata)]
  library(zoo)
  proxy = log(na.approx(proxy0))
  
  GISevents = eventdata[(eventdata$`NGRIP depth (m)`>min(depth))&(eventdata$`NGRIP depth (m)`<max(depth)),]
  
  eventdepths = GISevents$`NGRIP depth (m)`
  eventindexes = c(1,which.index(eventdepths, depth[2:length(depth)]) )
  
  plot(depth,proxy,type="l",xlab="Depth",ylab="log(Ca2+)"); abline(v=eventdepths,lwd=0.7,col="gray")
  
  reg.model = list(const=FALSE,depth1=FALSE,depth2=TRUE,proxy=TRUE,psi0=TRUE,psi1=TRUE)
  
  
  object_dust= main(age,depth,proxy,eventdepths,nsims=nsims,transform=transform,
                    noise="ar1",reg.model=reg.model,proxy.type="calcium")
  

  #d18O
  depth = maindata$`NGRIP-2 depth (m)`[2979:nrow(maindata)]
  proxy = maindata$`Delta O18 (permil)`[2979:nrow(maindata)]#dust, NA until 1167
  age = maindata$`GICC05 age (yr b2k)`[2979:nrow(maindata)]
  MCE = maindata$`GICC05 MCE (yr)`[2979:nrow(maindata)]
  
  GISevents = eventdata[(eventdata$`NGRIP depth (m)`>min(depth))&(eventdata$`NGRIP depth (m)`<max(depth)),]
  
  eventdepths = GISevents$`NGRIP depth (m)`
  eventindexes = c(1,which.index(eventdepths, depth[2:length(depth)]) )
  
  plot(depth,proxy,type="l",xlab="Depth",ylab="d18O"); abline(v=eventdepths,lwd=0.7,col="gray")
  
  reg.model = list( const=FALSE,depth1=FALSE,depth2=TRUE,proxy=TRUE,psi0=TRUE,psi1=TRUE)
  
  object_d18O= main(age,depth,proxy,eventdepths,nsims=nsims,transform=transform,
                    noise="ar1",reg.model=reg.model,proxy.type="d18O")



#allocate data
agesimmatrix_dust = matrix(NA,29,nsims)
agesimmatrix_d18O = matrix(NA,29,nsims)
depthlist_dust = c()
depthlist_d18O = c()
do.plot.rampfit = TRUE
do.plot.agehist = FALSE
depthstats = as.data.frame(matrix(NA,7,29))
colnames(depthstats) = c("true", "mean", "CIL","CIU","dustmean","dustCIL","dustCIU")
agestats = as.data.frame(matrix(NA,7,29))
colnames(agestats) = c("true", "mean", "CIL","CIU","dustmean","dustCIL","dustCIU")

steplengths = rep(0.01,29)
steplengths[21]=0.001 #sometimes changing the steplength slightly might improve convergence (default is 0.01)

#iterate over all 29 transitions
for(eventnumber in 1:29){
  #find data windows
  lowerints = which.index(event_intervals$depth_int_lower.m, object_d18O$data$z)
  upperints = which.index(event_intervals$depth_int_upper.m, object_d18O$data$z)
  interval_d18O = lowerints[eventnumber]:upperints[eventnumber]
  
  lowerints = which.index(event_intervals$depth_int_lower.m, object_dust$data$z)
  upperints = which.index(event_intervals$depth_int_upper.m, object_dust$data$z)
  interval_dust = lowerints[eventnumber]:upperints[eventnumber]
  
  depth.reference = event_intervals$NGRIP_depth_m[eventnumber]
  age.reference = event_intervals$GICC_age.yb2k[eventnumber]
  
  #specify parameters for transition onset age estimation
  event.estimation = list(interval_dust=interval_dust,interval_d18O=interval_d18O,
                          t1.sims=50000,rampsims=50000,h=steplengths[eventnumber],
                          label=paste0(eventnumber," (",event_intervals[eventnumber,1],"): ",event_intervals[eventnumber,2]),
                          depth.reference=event_intervals$NGRIP_depth_m[eventnumber],
                          age.reference=event_intervals$NGRIP_agee_m[eventnumber])
  
  #fit linear ramp
  eventobject_dust = linrampfitter(object_dust,interval=event.estimation$interval_dust,
                           h=event.estimation$h,t1.sims=event.estimation$t1.sims,
                           rampsims=event.estimation$rampsims,
                           label=paste0("Dust: ",event.estimation$label),
                           depth.reference = event.estimation$depth.reference,
                           print.progress=print.progress)
  
  if(do.plot.rampfit){ #plot linear ramp fit
    plot_results(eventobject_dust,plot.proxydata=NULL,plot.ls = NULL,
                 plot.inla.posterior = NULL,plot.inlasims = NULL,plot.bias = NULL,
                 plot.linramp = list(depth.reference=NULL,show.t0=TRUE,show.t1=TRUE,xrev=TRUE,
                                     paste0("Dust: ",event.estimation$label)),
                 plot.event_depth = NULL,plot.event_age = NULL)
  }
  eventobject_d18O = linrampfitter(object_d18O,interval=event.estimation$interval_d18O,
                                h=event.estimation$h,t1.sims=event.estimation$t1.sims,
                                rampsims=event.estimation$rampsims,
                                label=paste0("d18O: ",event.estimation$label),
                                log.ramp=FALSE,
                                depth.reference = event.estimation$depth.reference,
                                print.progress=print.progress)
  if(do.plot.rampfit){
    plot_results(eventobject_d18O,plot.proxydata=NULL,plot.ls = NULL,
                 plot.inla.posterior = NULL,plot.inlasims = NULL,plot.bias = NULL,
                 plot.linramp = list(depth.reference=NULL,show.t0=TRUE,show.t1=TRUE,xrev=TRUE,
                                     paste0("d18O: ",event.estimation$label)),
                 plot.event_depth = NULL,plot.event_age = NULL)
  }
  
  #sample onset age
  eventobject_dust = event_depth_to_age(eventobject_dust, nsims = nsims, print.progress=print.progress,
                             label=event.estimation$label,
                             age.reference = event.estimation$age.reference)
  agesimmatrix_dust[eventnumber,] = eventobject_dust$event_dating$samples
  if(do.plot.agehist){ #plot histogram of simulated onset ages
    hist(agesimmatrix_dust[eventnumber,],freq=0,breaks=50,col="orange",main=paste0("Dust: ",event.estimation$label))
  }
  eventobject_d18O = event_depth_to_age(eventobject_d18O, nsims = nsims, print.progress=print.progress,
                                  label=event.estimation$label,
                                  age.reference = event.estimation$age.reference)
  agesimmatrix_d18O[eventnumber,] = eventobject_d18O$event_dating$samples
  if(do.plot.agehist){
    hist(agesimmatrix_d18O[eventnumber,],freq=0,breaks=50,col="orange",main=paste0("d18O: ",event.estimation$label))
  }
  
  #compute summary statistics
  depthstats[1,eventnumber] = event_intervals$NGRIP_depth_m[eventnumber]
  depthstats[2,eventnumber] = eventobject_d18O$linramp$param$t0$mean
  depthstats[3,eventnumber] = eventobject_d18O$linramp$param$t0$q0.025
  depthstats[4,eventnumber] = eventobject_d18O$linramp$param$t0$q0.975
  depthstats[5,eventnumber] = eventobject_d18O$linramp$param$t0$mean
  depthstats[6,eventnumber] = eventobject_d18O$linramp$param$t0$q0.025
  depthstats[7,eventnumber] = eventobject_d18O$linramp$param$t0$q0.975
  
  densage1 = density(agesimmatrix_d18O[eventnumber,]); densage1 = data.frame(x=densage1$x,y=densage1$y)
  densage2 = density(agesimmatrix_dust[eventnumber,]); densage2 = data.frame(x=densage2$x,y=densage2$y)
  zage1 = inla.zmarginal(densage1,silent=TRUE)
  zage2 = inla.zmarginal(densage2,silent=TRUE)
  
  agestats[1,eventnumber] = event_intervals$GICC_age.yb2k[eventnumber]
  agestats[2,eventnumber] = mean(agesimmatrix_d18O[eventnumber])
  agestats[3,eventnumber] = zage1$quant0.025
  agestats[4,eventnumber] = zage1$quant0.975
  agestats[5,eventnumber] = mean(agesimmatrix_dust[eventnumber])
  agestats[6,eventnumber] = zage1$quant0.025
  agestats[7,eventnumber] = zage1$quant0.975
  
  
}

#import onsets from other sources to compare with
rasmussen_ages_full = read_excel("datasets_used/Rasmussen_et_al_2014_QSR_Table_2.xlsx",
                                 skip=21)

rasmussen_ages = event_intervals$GICC_age.yb2k #GISevents$`Age (a b2k)`

buizert_onsets = read_excel("datasets_used/Buizert_onsets.xlsx",
                            col_types=c("text","numeric","numeric","numeric","numeric"))
buizert_depths = buizert_onsets$Depth; buizert_ages = buizert_onsets$Age
#buizert_depths[c(2:6,8:9,12:13,15:17,19:21)]
buizert_ages = buizert_ages[c(2:6,8:9,12:13,15:17,19:21)]

capron_onsets_NGRIP_d18O = read_excel("datasets_used/Capron_onsets.xls",
                           sheet="d18O",n_max=25)
capron_ages_NGRIP_d18O = capron_onsets_NGRIP_d18O$`t1 (50%)`

capron_onsets_NGRIP_dust = read_excel("datasets_used/Capron_onsets.xls",
                                sheet="Ca2+",n_max=24)
capron_ages_NGRIP_dust = capron_onsets_NGRIP_dust$`t1 (50%)`


capron_onsets_NEEM_d18O = read_excel("datasets_used/Capron_onsets.xls",
                                      sheet="d18O",skip=26)
capron_ages_NEEM_d18O = capron_onsets_NEEM_d18O$`t1 (50%)`

capron_onsets_NEEM_dust = read_excel("datasets_used/Capron_onsets.xls",
                                      sheet="Ca2+",skip=25)
capron_ages_NEEM_dust = capron_onsets_NEEM_dust$`t1 (50%)`


#visualize dating uncertainty with comparisons to Rasmussen, Buizert and Capron
library(stringr)
par(mfrow=c(5,6),mar=c(2,2,1.5,1))
capind = numeric(29); buiind = numeric(29)
event_intervals[,2]
t(capron_onsets_NGRIP_d18O[,1])
capind = c(NA,2,3, 4,5,6,NA,NA,7,8,NA,9,10,11,NA,NA,NA,NA,NA,12,13,14,NA,NA,15,NA,NA,16,17)
buiind = c(NA,2,NA,3,4,6,NA,NA,8,9,NA,11,12,13,NA,NA,NA,NA,NA,15,16,17,NA,NA,19,NA,NA,20,21)
for(i in 1:29){
  dens1 = density(agesimmatrix_dust[i,]); dens1 = data.frame(x=dens1$x,y=dens1$y)#/max(dens1$y))
  dens2 = density(agesimmatrix_d18O[i,]); dens2 = data.frame(x=dens2$x,y=dens2$y)#/max(dens2$y))
  xlim = range(dens1$x,dens2$x);ylim = range(dens1$y,dens2$y)
  #par(mar=c(2,2,1.5,1))
  
  
  plot(dens1,type="l",col=1,xlab="Onset year (b2k)",ylab="Density",
       main=str_sub(event_intervals[i,2],str_locate(event_intervals[i,2],"GI-")[1]),
       xlim=xlim,ylim=ylim)
  #Axis(side=1,labels=FALSE,at=dens1$x,tck=0.1)
  
  lines(dens2,type="l",col="gray")
  if(FALSE){#95% credible intervals
    abline(v=mean(agesimmatrix_dust[i,]),col=1)
    zdens1 = inla.zmarginal(dens1,silent = TRUE)
    abline(v=c(zdens1$quant0.025,zdens1$quant0.975),col=1,lty=3,lwd=0.8)
    abline(v=mean(agesimmatrix_d18O[i,]),col=2)
    zdens2 = inla.zmarginal(dens2,silent = TRUE)
    abline(v=c(zdens2$quant0.025,zdens2$quant0.975),col=2,lty=3,lwd=0.8)
  }
  
  abline(v=event_intervals$GICC_age.yb2k[i],col="blue")  
  if(!is.na(buiind[i])){
    abline(v=buizert_onsets$Age[buiind[i]],col="green")
    #abline(v=buizert_ages[buiind[i]],col="green")
  }
  
  
  # abline(v=capron_ages_NGRIP_d18O,col="orange")
  # abline(v=capron_ages_NGRIP_dust,col="orange",lty=3)
  if(!is.na(capind[i])){
    abline(v=capron_ages_NGRIP_d18O[capind[i]],col="red")
    abline(v=capron_ages_NGRIP_dust[capind[i]],col="pink",lty=1)
  }
  
  
}
par(mar=c(0,0,0,0));plot(-1,axes=FALSE,xlab="",ylab="",xlim=c(0,1),ylim=c(0,1))
legend(0,1,legend=c("d18O onset posterior","Ca2+ onset posterior",
                    "Rasmussen onset","Buizert onset",
                    "Capron d18O onset","Capron Ca2+ onset"),
       col=c("black","gray", "blue", "green", "red", "pink"),
       lty=c(1,1,1,1,1,1), cex=0.65,bty="n")

par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1)



