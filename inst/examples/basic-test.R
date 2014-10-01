## test code for a simple model
ncity <- 50
eventnames <- c("birth", "latent", "infect", "recover", "deltaR"); 
statenames <- c("S", "E", "I", "R", "N")
## initial conditions, vary I
initstates = data.frame(S=7.5e4, E=0, I=2*(1:ncity), R=1e6)
initstates$N = rowSums(initstates)
## states by row
initstates = t(as.matrix(initstates))
## markov state transition matrix 
transmat <- matrix(0, nrow=length(statenames), ncol=length(eventnames)); 
rownames(transmat) <- rownames(initstates) <- statenames
colnames(transmat) <- eventnames
## markov state transition matrix definitions
transmat["S", "birth"] = 1
transmat["S", "latent"] = -1
transmat["E", "latent"] = 1
transmat["E", "infect"] = -1
transmat["I", "infect"] = 1
transmat["I", "recover"] = -1
transmat["R", "recover"] = 1
transmat["R", "deltaR"] = 1
transmat["N", "birth"] = 1
transmat["N", "deltaR"] = 1

accumvars = c("latent", "infect")
obsall = T
obs_nstep = 7  ## per week 
nsteps = 30*365  ## 30 years
deltat = 1/365  ## timestep, not used

## initialize model
## see function newSEIRModel for a pre-made version
## here everything is done by hand
mymod <- newModel(initstates, transmat, accumvars, obsall, nsteps, obs_nstep, deltat)

## each city gets its own list of model parameters
modlist <- lapply( 1:ncity, function(x) list( R0=16,
                    probs=0, distmethod=0,
                    betaforce=0.1,
                    imports=10^-5.5, schoollag=0,
                    ## schooltype 0 = sin, type 1 = term
                    schooltype=c(0), importmethod=c(1),
                    birth=0.05/365, deltaR=0.01/365, 
                    deltat=1, 
                    sigma=1/8, gamma=1/5
))

## initialize metapop and population parameters
mymod$set_metapop(list(dummy=1))
mymod$set_pop(modlist)

## run the simulation
mymod$steps(nsteps)
## pull out the first 7 model states into a list
basic.results <- lapply(1:7, function(x) mymod$get_metapop_state( x));

## 4 years worth of weekly observations
.plot.obs <- 1:(52*4)
.plot.state <- 1
## for each city, 
## plot a timeseries of state .plot.state over time .plot.obs
require(lattice)
plot(
    levelplot(basic.results[[.plot.state]][.plot.obs, ], aspect='fill')
)
