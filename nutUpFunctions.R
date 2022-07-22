# function code from 'Fite experimenta data.r'

if (!require(ReacTran)) install.packages('ReacTran')
library(ReacTran)
if (!require(nlmrt)) install.packages('nlmrt')
library(nlmrt)
if (!require(minpack.lm)) install.packages('minpack.lm')
library(minpack.lm)

NutUpTS <- function (times, y, parm){
  # Advection-dispersion with transient storage and first order uptake in main channel
  require(ReacTran)
  Cs <- y[1:200]
  Cts <- y[201:400]
  AFDW <- fiadeiro(v = parm[3], D = parm[2], grid = grid)
  trans <- tran.1D(C = Cs,
                   C.up = parm[1],
                   C.down = parm[1],
                   D = parm[2],
                   v = parm[3],
                   dx = grid,
                   AFDW = AFDW)
  ts.exchange <- parm[5] * (Cs-Cts)
  uptake <- parm[4] * Cs
  input <- parm[4] * parm[1]
  dCs <- trans$dC - ts.exchange - uptake + input
  dCts <- parm[6] * ts.exchange
  list(c(dCs, dCts))
}

BTCTS <- function(parm){
  # Fit Cl and N together, assuming same residual variance
  # require(deSolve)
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

ResidualTS <- function(parm){
  simu <- BTCTS(parm)
  resid <- c(data$Cl, data$N) - simu
  return(resid)
}


NutUp <- function (times, y, parm){
  #  Advection-dispersion with first order uptake in main channel
  Cs <- y
  AFDW <- fiadeiro(v = parm[3],D = parm[2],grid = grid)
  trans <- tran.1D(C = Cs, C.up = parm[1],
                   C.down = parm[1],
                   D = parm[2],
                   v = parm[3],
                   dx = grid,
                   AFDW = AFDW)
  uptake <- parm[4]*Cs
  input <- parm[4]*parm[1]
  dCs <- trans$dC - uptake + input
  list(c(dCs))
}


BTC <- function(parm){
  parm.Cl <- c(parm[1],parm[3],parm[4],0)
  parm.N <- c(parm[2],parm[3],parm[4],parm[5])
  yini.Cl <- c(parm.Cl[1],rep(Cl.ini,initial.length),rep(parm.Cl[1],(199-initial.length)))
  yini.N <- c(parm.N[1],rep(N.ini,initial.length),rep(parm.N[1],(199-initial.length)))
  Cl.conc <- ode.1D(func = NutUp,y = yini.Cl, parms = parm.Cl, times = time, nspec = 1)[,2*distance]
  N.conc <- ode.1D(func = NutUp,y = yini.N, parms = parm.N, times = time, nspec = 1)[,2*distance]
  return(c(Cl.conc,N.conc))
}

Residual <- function(parm){
  simu <- BTC(parm)
  resid <- c(data$Cl, data$N) - simu
  return(resid)
}


NutUpMM <- function (times, y, parm){
  # Advection-dispersion with transient storage and M-M uptake in main channel
  Cs <- y
  AFDW <- fiadeiro(v = parm[3],D = parm[2],grid = grid)
  trans <- tran.1D(C = Cs, C.up = parm[1], C.down = parm[1], D = parm[2], v = parm[3], dx = grid, AFDW = AFDW)
  uptake <- parm[4]*Cs/(Cs+parm[5])
  input <- parm[4]*parm[1]/(parm[1]+parm[5])
  dCs <- trans$dC - uptake + input
  list(c(dCs))
}


BTCMM <- function(parm){
  parm.Cl <- c(parm[1],parm[3],parm[4],0,0)
  parm.N <- c(parm[2],parm[3],parm[4],parm[5],parm[6])
  yini.Cl <- c(parm.Cl[1],rep(Cl.ini,initial.length),rep(parm.Cl[1],(199-initial.length)))
  yini.N <- c(parm.N[1],rep(N.ini,initial.length),rep(parm.N[1],(199-initial.length)))
  Cl.conc <- ode.1D(func = NutUpMM,y = yini.Cl, parms = parm.Cl, times = time, nspec = 1)[,2*distance]
  N.conc <- ode.1D(func = NutUpMM,y = yini.N, parms = parm.N, times = time, nspec = 1)[,2*distance]
  return(c(Cl.conc,N.conc))
}

ResidualMM <- function(parm){
  simu <- BTCMM(parm)
  resid <- c(data$Cl, data$N) - simu
  return(resid)
}


# Advection-dispersion with transient storage and M-M uptake in main channel #
NutUpTSMM <- function (times, y, parm){
  Cs <- y[1:200]
  Cts <- y[201:400]
  AFDW <- fiadeiro(v = parm[3],D = parm[2],grid = grid)
  trans <- tran.1D(C = Cs, C.up = parm[1], C.down = parm[1], D = parm[2], v = parm[3], dx = grid, AFDW = AFDW)
  ts.exchange <- parm[5]*(Cs-Cts)
  uptake <- parm[4]*Cs/(Cs+parm[7])
  input <- parm[4]*parm[1]/(parm[1]+parm[7])
  dCs <- trans$dC - ts.exchange - uptake + input
  dCts <- parm[6]*ts.exchange
  list(c(dCs, dCts))
}

BTCTSMM <- function(parm){
  parm.Cl <- c(parm[1],parm[3],parm[4],0,parm[5],parm[6],0)
  parm.N <- c(parm[2],parm[3],parm[4],parm[7],parm[5],parm[6],parm[8])
  yini.Cl <- c(parm.Cl[1],rep(Cl.ini,initial.length),rep(parm.Cl[1],(199-initial.length)), rep(parm.Cl[1],200))
  yini.N <- c(parm.N[1],rep(N.ini,initial.length),rep(parm.N[1],(199-initial.length)), rep(parm.N[1],200))
  Cl.conc <- ode.1D(func = NutUpTSMM,y = yini.Cl, parms = parm.Cl, times = time, nspec = 2)[,2*distance]
  N.conc <- ode.1D(func = NutUpTSMM,y = yini.N, parms = parm.N, times = time, nspec = 2)[,2*distance]
  return(c(Cl.conc,N.conc))
}

ResidualTSMM <- function(parm){
  simu <- BTCTSMM(parm)
  resid <- c(data$Cl, data$N) - simu
  return(resid)
}
