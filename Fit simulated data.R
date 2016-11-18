#Simulate data set. Demonstrate issues with Covino method#
library(ReacTran)
library(nlmrt)
library(minpack.lm)

#Advection-dispersion with transient storage and first order uptake in main channel#
NutUpTS = function (times, y, parm){
	Cs = y[1:200]
	Cts = y[201:400]
	AFDW = fiadeiro(v=parm[3],D=parm[2],grid=grid)
	trans = tran.1D(C=Cs, C.up=parm[1], C.down=parm[1], D=parm[2], v=parm[3], dx=grid, AFDW=AFDW)
	ts.exchange = parm[5]*(Cs-Cts)
	uptake = parm[4]*Cs
	input = parm[4]*parm[1]
	dCs = trans$dC - ts.exchange - uptake + input
	dCts = parm[6]*ts.exchange
	list(c(dCs, dCts))
}

#Fit Cl and N together, assuming same residual variance#
BTCTS = function(parm){
	parm.Cl = c(parm[1],parm[3],parm[4],0,parm[5],parm[6])
	parm.N = c(parm[2],parm[3],parm[4],parm[7],parm[5],parm[6])
	yini.Cl = c(parm.Cl[1],rep(Cl.ini,initial.length),rep(parm.Cl[1],(199-initial.length)), rep(parm.Cl[1],200))
	yini.N = c(parm.N[1],rep(N.ini,initial.length),rep(parm.N[1],(199-initial.length)), rep(parm.N[1],200))
	Cl.conc = ode.1D(func=NutUpTS,y=yini.Cl, parms=parm.Cl, times=time, nspec=2)[,2*distance]
	N.conc = ode.1D(func=NutUpTS,y=yini.N, parms=parm.N, times=time, nspec=2)[,2*distance]
	return(c(Cl.conc,N.conc))
}

ResidualTS = function(parm){
	simu = BTCTS(parm)
	resid = c(data$Cl, data$N) - simu
	return(resid)
}


#Advection-dispersion with first order uptake in main channel#
NutUp = function (times, y, parm){
	Cs = y
	AFDW = fiadeiro(v=parm[3],D=parm[2],grid=grid)
	trans = tran.1D(C=Cs, C.up=parm[1], C.down=parm[1], D=parm[2], v=parm[3], dx=grid, AFDW=AFDW)
	uptake = parm[4]*Cs
	input = parm[4]*parm[1]
	dCs = trans$dC - uptake + input
	list(c(dCs))
}


BTC = function(parm){
	parm.Cl = c(parm[1],parm[3],parm[4],0)
	parm.N = c(parm[2],parm[3],parm[4],parm[5])
	yini.Cl = c(parm.Cl[1],rep(Cl.ini,initial.length),rep(parm.Cl[1],(199-initial.length)))
	yini.N = c(parm.N[1],rep(N.ini,initial.length),rep(parm.N[1],(199-initial.length)))
	Cl.conc = ode.1D(func=NutUp,y=yini.Cl, parms=parm.Cl, times=time, nspec=1)[,2*distance]
	N.conc = ode.1D(func=NutUp,y=yini.N, parms=parm.N, times=time, nspec=1)[,2*distance]
	return(c(Cl.conc,N.conc))
}

Residual = function(parm){
	simu = BTC(parm)
	resid = c(data$Cl, data$N) - simu
	return(resid)
}

#Advection-dispersion with transient storage and M-M uptake in main channel#
NutUpMM = function (times, y, parm){
	Cs = y
	AFDW = fiadeiro(v=parm[3],D=parm[2],grid=grid)
	trans = tran.1D(C=Cs, C.up=parm[1], C.down=parm[1], D=parm[2], v=parm[3], dx=grid, AFDW=AFDW)
	uptake = parm[4]*Cs/(Cs+parm[5])
	input = parm[4]*parm[1]/(parm[1]+parm[5])
	dCs = trans$dC - uptake + input
	list(c(dCs))
}


BTCMM = function(parm){
	parm.Cl = c(parm[1],parm[3],parm[4],0,0)
	parm.N = c(parm[2],parm[3],parm[4],parm[5],parm[6])
	yini.Cl = c(parm.Cl[1],rep(Cl.ini,initial.length),rep(parm.Cl[1],(199-initial.length)))
	yini.N = c(parm.N[1],rep(N.ini,initial.length),rep(parm.N[1],(199-initial.length)))
	Cl.conc = ode.1D(func=NutUpMM,y=yini.Cl, parms=parm.Cl, times=time, nspec=1)[,2*distance]
	N.conc = ode.1D(func=NutUpMM,y=yini.N, parms=parm.N, times=time, nspec=1)[,2*distance]
	return(c(Cl.conc,N.conc))
}

ResidualMM = function(parm){
	simu = BTCMM(parm)
	resid = c(data$Cl, data$N) - simu
	return(resid)
}


#Advection-dispersion with transient storage and M-M uptake in main channel #
NutUpTSMM = function (times, y, parm){
	Cs = y[1:200]
	Cts = y[201:400]
	AFDW = fiadeiro(v=parm[3],D=parm[2],grid=grid)
	trans = tran.1D(C=Cs, C.up=parm[1], C.down=parm[1], D=parm[2], v=parm[3], dx=grid, AFDW=AFDW)
	ts.exchange = parm[5]*(Cs-Cts)
	uptake = parm[4]*Cs/(Cs+parm[7])
	input = parm[4]*parm[1]/(parm[1]+parm[7])
	dCs = trans$dC - ts.exchange - uptake + input
	dCts = parm[6]*ts.exchange
	list(c(dCs, dCts))
}

BTCTSMM = function(parm){
	parm.Cl = c(parm[1],parm[3],parm[4],0,parm[5],parm[6],0)
	parm.N = c(parm[2],parm[3],parm[4],parm[7],parm[5],parm[6],parm[8])
	yini.Cl = c(parm.Cl[1],rep(Cl.ini,initial.length),rep(parm.Cl[1],(199-initial.length)), rep(parm.Cl[1],200))
	yini.N = c(parm.N[1],rep(N.ini,initial.length),rep(parm.N[1],(199-initial.length)), rep(parm.N[1],200))
	Cl.conc = ode.1D(func=NutUpTSMM,y=yini.Cl, parms=parm.Cl, times=time, nspec=2)[,2*distance]
	N.conc = ode.1D(func=NutUpTSMM,y=yini.N, parms=parm.N, times=time, nspec=2)[,2*distance]
	return(c(Cl.conc,N.conc))
}

ResidualTSMM = function(parm){
	simu = BTCTSMM(parm)
	resid = c(data$Cl, data$N) - simu
	return(resid)
}


initial.length=1
Cl.ini = 1000
N.ini = 2000
grid = setup.grid.1D(L=100, N=200)
time = c(3,14, seq(from=18, to=38, by=0.5), seq(from=40,to=160,by=5))
parameter = c(Cl=6,N=3,D=0.3,U=1.5,alpha=0.03,AsA=2,K=0.05)
distance = 40
simu.data = BTCTS(parameter) + rnorm(n=2*length(time),mean=0, sd=0.5)
data = data.frame(Cl=simu.data[1:length(time)], N=simu.data[-(1:length(time))])


mod1 = nlfb(start=c(Cl=6,N=3,D=0.3,U=1.5,K=0.05), resfn=Residual, lower=c(0,0,0,0,0), upper=c(12,5,Inf,Inf,Inf))
mod2 = nlfb(start=c(Cl=6,N=3,D=0.3,U=1.5,alpha=0.03,AsA=2,K=0.05), resfn=ResidualTS, lower=c(0,0,0,0,0,0,0), upper=c(12,5,Inf,Inf,Inf,Inf,Inf))
mod3 = nlfb(start=c(Cl=6,N=3,D=0.3,U=1.5,Vmax=24000, Km=480000), resfn=ResidualMM, lower=c(0,0,0,0,0,0), upper=c(12,5,Inf,Inf,Inf,Inf))
mod4 = nlfb(start=c(Cl=6,N=3,D=0.3,U=1.5,alpha=0.03,AsA=2,Vmax=24000,Km=480000), resfn=ResidualTSMM, upper=c(10,5,Inf,Inf,Inf,Inf,Inf,Inf), lower=c(0,0,0,0,0,0,0,0))
