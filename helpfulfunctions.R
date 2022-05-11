# Contains functions used by other functions


meanmaker = function(coefs,reg.model,nevents=69,data){
  ## Computes mean vector from given fixed effects 'coefs'.
  ## Requires specification of which effects to include ('reg.model'), the number of climate transitions ('nevents') and a data.frame with covariates ('data')
  coefcounter=1
  fitted=numeric(dim(data)[1])
  if(reg.model$const){
    fitted=coefs[1]
    coefcounter=coefcounter+1
  }
  if(reg.model$depth1){
    fitted = fitted + coefs[coefcounter]*data$z
    coefcounter=coefcounter+1
  }
  if(reg.model$depth2){
    fitted = fitted + coefs[coefcounter]*data$z2
    coefcounter=coefcounter+1
  }
  if(reg.model$proxy){
    fitted=fitted + coefs[coefcounter]*data$x
    coefcounter=coefcounter+1
  }
  if(nevents>0){
    for(i in 2:nevents){
      if(reg.model$psi0){
        
        fitted = fitted + coefs[coefcounter]*data[[paste0("a",i-1)]]
        coefcounter=coefcounter+1
      }
      if(reg.model$psi1){
        
        fitted = fitted + coefs[coefcounter]*data[[paste0("c",i-1)]]
        coefcounter=coefcounter+1
      }
    }
  }
  return(fitted)
}

# Computes the (noiseless) linear ramp function.
linramp = function(t,t0=0,dt=1,y0=0,dy=1){
  y = numeric(length(t))
  y = y0 + dy*(t-t0)/dt
  y[t<t0]=y0
  y[t>t0+dt]=y0+dy
  return(y)
}


# Computes the (noiseless) linear ramp function, but in reverse.
linramprev = function(t,t0=0,dt=1,y0=0,dy=1){
  y = numeric(length(t))
  y = y0 + dy*(t-t0)/dt

  y[t>t0]=y0
  y[t<t0+dt]=y0+dy
  return(y)
}


# Finds the indices where 'events' best match values in a given 'record'.
which.index = function(events, record){
  eventindexes = numeric(length(events))
  for(i in 1:length(events)){
    if(events[i] < min(record) || events[i] > max(record)){ #Gives NA if located outside range of 'record'
      warning(paste0("Event ",i,", located at ",events[i]," is outside the interval covered by 'record' (",min(record),", ",max(record),"). The event will be omitted!"))
      eventindexes[i] = NA
    }else{
      eventindexes[i] = which(abs(events[i]-record) == min(abs(events[i]-record)))
    }

  }
  #eventindexes = unique(c(1,eventindexes[!is.na(eventindexes)])) #Placing transition at the start of record. Removing NA and duplicates
  return(eventindexes)
}

# Computes posterior marginal mean and uncertainty intervals from simulations.
simulationsummarizer = function(object,CI.type="hpd",print.progress=FALSE){
  if(print.progress) cat("Computing posterior marginal mean and 95% credible intervals from chronology samples...\n",sep="")
  time.start = Sys.time()
  n = dim(object$simulation$age)[1]
  nsims = dim(object$simulation$age)[2]

  meanvek = rowMeans2(object$simulation$age)
  sdvek = sqrt(rowVars(object$simulation$age))

  lower = numeric(n); upper = numeric(n)
  if(CI.type=="hpd"){
    modevek = numeric(n)
    for(i in 1:n){
      dens = density(object$simulation$age[i,])
      modevek[i]=dens$x[which(dens$y == max(dens$y))]

      lower[i] = inla.hpdmarginal(0.95,dens)[1]
      upper[i] = inla.hpdmarginal(0.95,dens)[2]

    }
  }else{
    lower = meanvek-1.96*sdvek
    upper = meanvek+1.96*sdvek
  }

  time.summary = Sys.time()
  object$simulation$summary = list(mean=meanvek,sd=sdvek,lower=lower,upper=upper,
                                   .args=list(interval=cbind(lower,upper),print.progress=print.progress,CI.type=CI.type))
  if(CI.type=="hpd") object$simulation$summary$mode = modevek
  if(print.progress) cat(" completed in ",difftime(time.summary,time.start,units="secs")[[1]],"\n",sep="")

  object$time$samplesummary = list(total=difftime(time.summary,time.start,units="secs")[[1]])
  object$simulation$summary$sim.sum.time = difftime(time.summary,time.start,units="secs")[[1]]

  return(object)
}



## sets initial values for fixed parameters equal to least squares 'fit'
control.fixed.priors = function(reg.model, fit, nevents){

  my.control.fixed = list(mean=list(  ))

  if(reg.model$depth1) my.control.fixed$mean[["z1"]] = fit$coefficients[["z1"]]
  if(reg.model$depth2) my.control.fixed$mean[["z2"]] = fit$coefficients[["z2"]]
  if(reg.model$proxy) my.control.fixed$mean[["x"]] = fit$coefficients[["x"]]

  if(reg.model$psi0 || reg.model$psi1){
    for(i in 1:(nevents-1)){
      if(reg.model$psi0){
        my.control.fixed$mean[[paste0("a",i)]] = fit$coefficients[[paste0("a",i)]]
      }
      if(reg.model$psi0){
        my.control.fixed$mean[[paste0("c",i)]] = fit$coefficients[[paste0("c",i)]]
      }
    }
  }

  return(my.control.fixed)
}








