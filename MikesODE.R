# A translation from python of Mike's simple calculation of genetic diversity
# given different deme structure and migration regularities.
plot.new()
frame()

## Parameters
# Number of demes
L = 5
# Individuals per deme
N = 50
# Original population size
Ntot0 = 500*100
# Mutation rate
mu = 1e-8
# Original heterozygosity 
H0 = 4 * Ntot0 * mu
# Array of times to plot
times <- seq(0, 2000, by=1)
# yini <- c(H0,H0) # initial conditions - within-deme and general H are the steady state
# before the population crash.
yini <- c(H0,H0,H0,H0)


## Calculation of heterozygosity trajectories
  # Function that returns (dHs/dt, dHt/dt)
  # Assumptions: infinite alleles per site.
  # asserions: two lineages can coalesce only if they are in the same deme.
  # Derivation: backwards-coalescence frame of thought. 
  # Within-demes:
  # H(t+1) = Ht(1-1/2N-2m-2Miu) + 2mHgen(t) + 2miu(*), whose reordering gives dH/timestep,
  # As apears in the function below. [note that the case of coalescence doesn't have a
  # term of its own, because its contribution to H is zero]
  # * : If we assume an infinite allele model, any mutation always leads to a contribution 
  # to H. If we assume a 2-allele model with back mutation, we need to add the terms:
  # 2miu(1-Ht) - 2miu*Ht, where the first is a mutation in a pair that were identical, 
  # and the second a back-mutation in an initially non-identical pair.
  # Hgeneral (HG): we assume no two events occur in a single timestep, so ignore the option
  # of migration and coalescence or mig+mutation, and we even ignore the option that both
  # individuals happen to be from the same deme and one of them had a mutation (although
  # this last thing isn't completely justified, I think).
  # HG(t+1) = HG(t)(1-1/L-2mu) + 1/L * H(t) + 2mu
  # HG(t+1)-HG(t) = 1/L(H(t)-HG(t)) + 2mu(1-HG(t))
  # Which looks, admittedly, quite different from Mike's derivation which is based on a
  # forward-in-time-in-his-head intuition/derivation-that-I-didn't-follow.
  library(deSolve)
  heterozygosities <- function(times, y, parms){
  L = parms[1]; N = parms[2]; mu = parms[3]; m = parms[4];
  dH_within      = -y[1]/(2*N) + 2*m*(y[2]-y[1]) + 2*mu*(1-y[1])
  dH_generalMike = -y[1]/(2*N*L) + 2*mu*(1-y[2])
  dH_withinOren  = -y[3]/(2*N) + 2*m*(y[4]-y[3]) + 2*mu*(1-y[3])
  dH_generalOren = (y[3]-y[4])/L + 2*mu*(1-y[4])
  # list(c(dH_within,dH_generalMike))
  list(c(dH_within,dH_generalMike,dH_withinOren,dH_generalOren))
  }

m = 0
out <- ode (times = times, y = yini, func = heterozygosities, parms = c(L,N,mu,m))
timeAxis = out[,1]; Hwithin = out[,2]; Hgeneral = out[,3]
plot(timeAxis,Hwithin, type="n", main='Heterozygosity in time') 
ylab="Heterozygosity"
par(col="blue", lwd = 1)
lines(timeAxis,Hwithin, lwd=1, lty = 2)
lines(timeAxis,Hgeneral, lwd = 1, lty = 1)

# timeAxis = out[,1]; HwithinOren = out[,3]; HgeneralOren = out[,4]
# par(col="yellow", lwd = 1)
# lines(timeAxis,HwithinOren, pch=22, lwd=1, lty = 3)
# lines(timeAxis,HgeneralOren, pch = 22, lwd = 1, lty = 3)
# # Seems like there's very little actual difference between the results for my 
# derivation of Hgeneral and Mike's.

m = 1/(10*N)
out <- ode (times = times, y = yini, func = heterozygosities, parms = c(L,N,mu,m))
timeAxis = out[,1]; Hwithin = out[,2]; Hgeneral = out[,3]
par(col="green", lwd = 1)
lines(timeAxis,Hwithin, lwd=1, lty = 2)
lines(timeAxis,Hgeneral, lwd = 1, lty = 1)


m = 1/N
out <- ode (times = times, y = yini, func = heterozygosities, parms = c(L,N,mu,m))
timeAxis = out[,1]; Hwithin = out[,2]; Hgeneral = out[,3]
par(col="red", lwd = 1)
lines(timeAxis,Hwithin, lwd=1, lty = 2)
lines(timeAxis,Hgeneral, lwd = 1, lty = 1)

m = 1
out <- ode (times = times, y = yini, func = heterozygosities, parms = c(L,N,mu,m))
timeAxis = out[,1]; Hwithin = out[,2]; Hgeneral = out[,3]
par(col="black", lwd = 1)
lines(timeAxis,Hwithin, lwd=1, lty = 2)
lines(timeAxis,Hgeneral, lwd = 1, lty = 1)

