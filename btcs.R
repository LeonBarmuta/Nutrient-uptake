BTCTS <- function(parm){
  # Fit Cl and N together, assuming same residual variance
  require(deSolve)
  parm.Cl <- c(parm[1],parm[3],parm[4],0,parm[5],parm[6])
  parm.N <- c(parm[2],parm[3],parm[4],parm[7],parm[5],parm[6])
  yini.Cl <- c(parm.Cl[1],rep(Cl.ini,initial.length),
               rep(parm.Cl[1],(199-initial.length)),
               rep(parm.Cl[1],200))
  yini.N <- c(parm.N[1],rep(N.ini,initial.length),
              rep(parm.N[1],(199-initial.length)),
              rep(parm.N[1],200))
  Cl.conc <- ode.1D(func = NutUpTS,y = yini.Cl,
                    parms = parm.Cl, times = time, nspec = 2)[,2*distance]
  N.conc <- ode.1D(func = NutUpTS,y = yini.N,
                   parms = parm.N, times = time, nspec = 2)[,2*distance]
  return(c(Cl.conc,N.conc))
}
