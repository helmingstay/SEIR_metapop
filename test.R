## test code
## 10 pops, 4 states

eventnames <- c("birth", "latent", "infect", "recover", "deltaR"); 
statenames <- c("S", "E", "I", "R", "N")
initstate = c(S=5e5, E=0, I=0, R=5e5, N=1e6)
initstates <- matrix(initstate, nrow=length(statenames), ncol=10)
transmat <- matrix(0, nrow=length(statenames), ncol=length(eventnames)); 
rownames(transmat) <- rownames(initstates) <- statenames
colnames(transmat) <- eventnames
transmat["S", "birth"] = 1
transmat["S", "latent"] = -1
transmat["E", "latent"] = 1
transmat["E", "infect"] = -1
transmat["I", "infect"] = 1
transmat["I", "recover"] = -1
transmat["R", "recover"] = 1
transmat["R", "deltaR"] = 1

accumvars = c("latent", "infect")
obsall = TRUE
nobs = 365*50
obs_nstep = 7
deltat = 1/365

mymod <- newModel(initstates, transmat, accumvars, obsall, nobs, obs_nstep, deltat)

modlist <- lapply( 1:ncol(initstates), function(x) list( R0=16,
                    probs=0, distmethod=c('null'),
                    betaforce=0.25,
                    imports=10^-5.5, schoollag=0,
                    ## schooltype 0 = sin, type 1 = term
                    schooltype=c(0), importmethod=c(1),
                    birth=0.02/365, deltaR=0.01/365, 
                    deltat=1, 
                    sigma=1/8, gamma=1/5
))

mymod$setpars(list(dummy=1), modlist)
 mymod$steps(7)
