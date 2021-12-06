# Prepares input and formats data for 'main.R' and functions therein.

prepare = function(age,depth,proxy, events=NULL,nsims=10000, eventmeasure = "depth",reg.model = list(
  const=FALSE,depth1=FALSE,depth2=TRUE,proxy=TRUE,psi0=TRUE,psi1=TRUE), noise="ar1", method="inla",reference.label=NULL){
  n = length(age)
  y = age[2:n]
  dy = diff(age)
  z = depth[2:n]

  data = data.frame(y=y,dy=dy,z=z)

  formulastring = "dy ~ "
  if(!reg.model$const) formulastring=paste0(formulastring,"-1") else formulastring=paste0(formulastring,"1")
  if(reg.model$depth1) formulastring=paste0(formulastring," + z")
  if(reg.model$depth2) {
    formulastring=paste0(formulastring," + z2")
    data[["z2"]]=z^2
  }
  if(reg.model$proxy) {
    if(missing(proxy)){
      stop("'proxy' is missing. Shutting down...")
    }
    formulastring=paste0(formulastring," + x")
    x = proxy[2:n]
    data[["x"]]=x
  }
  if(!is.null(events)){
    eventindexes = numeric(length(events))

    if(tolower(eventmeasure) %in% c("depth","z")){
      eventindexes = which.index (events, z)
    }else if(tolower(eventmeasure) %in% c("age","time","y")){
      eventindexes = which.index (events, y)
    }

    eventindexes = unique(c(1,eventindexes[!is.na(eventindexes)]))
    nevents = length(eventindexes)
  }
  if(reg.model$psi0) {
    if(is.null(events)){
      stop("'events' are missing. Shutting down...")
    }


    for(i in 2:nevents){
      ev = numeric(n-1)
      ev[ eventindexes[i-1]:(eventindexes[i]-1) ] = data$z[eventindexes[i-1]:(eventindexes[i]-1)]
      data[[paste0("a",i-1)]] = ev
      konst = numeric(n-1)
      konst[eventindexes[i-1]:(eventindexes[i]-1)] = 1
      data[[paste0("c",i-1)]] = konst
      formulastring = paste0(formulastring, " + a",i-1," + c",i-1)
    }
  }
  if(tolower(method)=="inla") data$idy=data$z

  object = list(data=data)
  object$.args = list(reg.model = reg.model, noise=noise,method=method, eventmeasure=eventmeasure)
  object$.args$ls.formulastring = formulastring
  object$.args$eventindexes = eventindexes
  object$ls.formula = as.formula(formulastring)
  object$.args$nevents=nevents
  object$preceeding = list(y0 = age[1],z0=depth[1],x0=proxy[1])

  if(tolower(method)=="inla"){
    if(tolower(noise) %in% c("iid","independent","ar(0)")){
      formulastring = paste0(formulastring, " + f(idy,model=\"iid\")")
    }else if(tolower(noise) %in% c("ar1","ar(1)")){
      formulastring = paste0(formulastring, " + f(idy,model=\"ar1\")")
    }else if(tolower(noise) %in% c("ar2","ar(2)")){
      formulastring = paste0(formulastring, " + f(idy,model=\"ar\",order=2)")
    }
  }
  formula = as.formula(formulastring)
  object$formula = formula
  object$.args$formulastring=formulastring
  object$.args$reference.label=reference.label

  ## move to separate function

  return(object)

}
