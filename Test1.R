source("nutUpFunctions.R")


# import to field data#
data.all <- read.csv("LUQ13E01TPost.csv", header = T)
width <- data.all$AvgWettedWidth_m[1]
depth <- data.all$AvgWettedDepth_cm[1] / 100
distance <- data.all$Reach.Length_meters[1]
discharge <- data.all$Discharge_LitersPerSec[1]

# Calculate injected Cl in g and N in mg#
inject.Cl <- data.all$Injected_NaCl_g[1] * 35.5 / (35.5 + 23) + data.all$Injected_NH4Cl_g[1] * 35.5 / (14 + 4 + 35.5)
inject.N <- data.all$Injected_NH4Cl_g[1] * 14 / (14 + 4 + 35.5) * 1000


# Assumping stock solution is released over one segment 0.5m#
initial.length <- 1

# Injectate concentration immediately after releasing. Cl in g/m3 and N in mg/m3#
Cl.ini <- inject.Cl / (width * depth * initial.length * 0.5)
N.ini <- inject.N / (width * depth * initial.length * 0.5)

# Setting up grid, calculate time since released for each sample#
grid <- setup.grid.1D(L = 100, N = 200)
collect.time <- strptime(data.all$CollectionTime, format = "%H:%M:%S")
start.time <- strptime(data.all$InjectionTime[1], format = "%H:%M:%S")
time <- as.numeric(collect.time - start.time)
data <- data.frame(Cl = data.all$ObservedCl_mgL, N = data.all$ObservedNH4N_ugL)

# Fit different models to data, calculate AIC for each model

# resfn argument gives model formula


mod1 <- nlfb(start = c(Cl = 11.7221, N = 4.81,
                       D = 1.1167, U = 1.107, K = 0.04515),
             resfn = Residual,
             lower = c(0, 0, 0, 0, 0),
             upper = c(12, 5, Inf, Inf, Inf))

mod2 <- nlfb(start = c(Cl = 9.18, N = 0.875,
                       D = 0.07537, U = 1.62454,
                       alpha = 0.1423,
                       AsA = 1.69, K = 0.0564),
             resfn = ResidualTS,
             lower = c(0, 0, 0, 0, 0, 0, 0),
             upper = c(12, 5, Inf, Inf, Inf, Inf, Inf))

# this takes a while
mod3 <- nlfb(start = c(Cl = 11.7221, N = 4.81,
                       D = 1.1167, U = 1.107,
                       Vmax = 45176, Km = 1000551),
             resfn = ResidualMM,
             lower = c(0, 0, 0, 0, 0, 0),
             upper = c(12, 5, Inf, Inf, Inf, Inf))

mod4 <- nlfb(start = c(Cl = 9.178, N = 0.87,
                       D = 0.075156, U = 1.6238,
                       alpha = 0.1419,
                       AsA = 1.69043,
                       Vmax = 56420,
                       Km = 1000000),
             resfn = ResidualTSMM,
             upper = c(12, 5, Inf, Inf, Inf, Inf, Inf, Inf),
             lower = c(0, 0, 0, 0, 0, 0, 0, 0))


str(mod1)

summary(mod1) # doesn't do anything?
class(mod1)

class(mod4)

summary.nls(mod4)
