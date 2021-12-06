
main = function(age,depth,proxy, events=NULL,nsims=10000, eventmeasure = "depth",reference.label=NULL,reg.model = list(
  const=FALSE,depth1=FALSE,depth2=TRUE,proxy=TRUE,psi0=TRUE,psi1=TRUE), noise="ar1", method="inla", CI.type="quantiles",
  DO.estimation = NULL, bias = NULL,store.everything=FALSE,print.progress=FALSE){
  #DO.estimation = list(interval=...,h=0.1,t1.sims=50000,rampsims=50000,label="GI-11,depth.reference=NULL,age.reference=NULL)
  #bias = list(bias.model="uniform",biasparams=cbind(c(1,1),c(0.98,1.02),c(0.96,1.04)),store.samples=FALSE)

  time.start = Sys.time()
  main.call = sys.call(which=1)
  if(print.progress) cat("Initiating data formatting...",sep="")

  object = prepare(age,depth,proxy, events=events,nsims=nsims, eventmeasure = eventmeasure,reg.model = reg.model,
                         noise=noise, method=method,reference.label=reference.label)
  if(print.progress) cat(" completed!\n",sep="")
  object = modelfitter(object, method=method,print.progress=print.progress)


  object = chronology_simulation(object, nsims=nsims, method=method,store.means=store.everything,print.progress=print.progress)
  object = simulationsummarizer(object,CI.type=CI.type,print.progress=print.progress)
  if(!is.null(DO.estimation)){
    object = linrampfitter(object,interval=DO.estimation$interval,h=DO.estimation$h,t1.sims=DO.estimation$t1.sims,
               rampsims=DO.estimation$rampsims,label=DO.estimation$label,depth.reference = DO.estimation$depth.reference,print.progress=print.progress)

    object = DO_depth_to_age(object, nsims = nsims, print.progress=print.progress,label=DO.estimation$label,age.reference = DO.estimation$age.reference)
  }
  if(!is.null(bias)){
    object = biased_chronologies(object,bias.model=bias$bias.model,biasparams = bias$biasparams,nsims=nsims,store.samples=store.everything)
  }
  time.total = difftime(Sys.time(), time.start,units="secs")[[1]]
  object$.args$call = main.call
  object$time$total = time.total

  return(object)
}

#remember to set proper directory
setwd("/home/myrvoll/Dropbox/Postdoc/Submission files/Uncertaintypaper/Code/")

source("biased_chronologies.R"); source("chronology_simulation.R"); source("DO_depth_to_age.R"); source("helpfulfunctions.R"); source("linrampfitter.R")
source("modelfitter.R"); source("plot_results.R"); source("prepare.R"); source("rgeneric_model.R"); source("summary_results.R")
library("INLA"); library("matrixStats"); library("numDeriv")


GISevents = read.table("datasets_used/GISevents.txt")
DO_intervals = read.table("datasets_used/DO_intervals.txt")
NGRIP5cm_d18O_GICC05 = read.table("datasets_used/NGRIP5cm_d18O_GICC05.txt")


depth = NGRIP5cm_d18O_GICC05$NGRIP_depth.m
proxy = NGRIP5cm_d18O_GICC05$delta_O18.permil
age = NGRIP5cm_d18O_GICC05$GICC05_age.yb2k
MCE= NGRIP5cm_d18O_GICC05$GICC05_MCE.yr



eventdepths = GISevents$NGRIP.depth..m.
eventindexes = which.index(eventdepths, depth[2:length(depth)])

nsims=10000; eventmeasure = "depth";reg.model = list(
  const=FALSE,depth1=FALSE,depth2=TRUE,proxy=TRUE,psi0=TRUE,psi1=TRUE); method="inla";
plots = list(posteriors=TRUE)

DOnumber=13

lowerints = which.index(DO_intervals$depth_int_lower.m, depth[2:length(depth)])
upperints = which.index(DO_intervals$depth_int_upper.m, depth[2:length(depth)])
depth.reference = DO_intervals$NGRIP_depth_m[DOnumber]
age.reference = DO_intervals$GICC_age.yb2k[DOnumber]
interval = lowerints[DOnumber]:upperints[DOnumber]

DO.estimation=list(interval=interval,make.t1=TRUE,t1.sims=50000,h=0.1,make.rampfit=TRUE,rampsims=50000)
oldwd=getwd()

nsims=10000

DO.estimation = list(interval=interval,t1.sims=50000,rampsims=50000,label="GI-11",depth.reference=DO_intervals$NGRIP_depth_m[DOnumber],age.reference=DO_intervals$NGRIP_agee_m[DOnumber])

print.progress=TRUE

noise = "ar1"

object = main(age,depth,proxy, events=eventdepths,nsims=10000, eventmeasure = "depth",reference.label="gicc05",reg.model = list(
  const=FALSE,depth1=FALSE,depth2=TRUE,proxy=TRUE,psi0=TRUE,psi1=TRUE), noise=noise, method="inla", CI.type="quantiles",
  DO.estimation = list(interval=interval,t1.sims=50000,rampsims=50000,label="GI-11",depth.reference=depth.reference,age.reference=age.reference),
  bias = list(bias.model="uniform",biasparams=cbind( c(1,1),c(0.98,1.02),c(0.96,1.04) ),store.samples=FALSE),store.everything=FALSE,print.progress=TRUE)

plot_results(object)
summary_results(object)

objectiid = modelfitter(object, method=method,print.progress=TRUE,noise="iid")
objectar1 = modelfitter(object, method=method,print.progress=TRUE,noise="ar1")
objectar2 = modelfitter(object, method=method,print.progress=TRUE,noise="ar2")

summary_results(objectiid); summary_results(objectar1); summary_results(objectar2)

layout(mat=matrix(c(1,2,4,7,3,5,8,9,6),nrow=3))

xrange_sigma = c(0.420,0.442)
xrange_phi = c(0.07,0.24)

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
