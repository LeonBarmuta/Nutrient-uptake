ResidualTS <- function(parm){
  simu <- BTCTS(parm)
  resid <- c(data$Cl, data$N) - simu
  return(resid)
}
