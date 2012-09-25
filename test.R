## test code

#ncity <- 2
ncity <- 50
eventnames <- c("birth", "latent", "infect", "recover", "deltaR"); 
statenames <- c("S", "E", "I", "R", "N")
initstate = c(S=5e4, E=0, I=10, R=5e5)
initstate['N'] = sum(initstate)
initstates <- matrix(initstate, nrow=length(statenames), ncol=ncity)
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
transmat["N", "birth"] = 1
transmat["N", "deltaR"] = 1

accumvars = c("latent", "infect")
obsall = T
obs_nstep = 7
#nsteps = 30
nsteps = 30*365
nobs = 1+(nsteps/obs_nstep)
deltat = 1/365

#!! add debug
mymod <- newModel(initstates, transmat, accumvars, obsall, nobs, obs_nstep, deltat)

modlist <- lapply( 1:ncity, function(x) list( R0=16,
                    probs=0, distmethod=0,
                    betaforce=0.25,
                    imports=10^-5.5, schoollag=0,
                    ## schooltype 0 = sin, type 1 = term
                    schooltype=c(0), importmethod=c(1),
                    birth=0.02/365, deltaR=0.01/365, 
                    deltat=1, 
                    sigma=1/8, gamma=1/5
))

mymod$set_metapop(list(dummy=1))
mymod$set_pop(modlist)

ps(mymod$steps(nsteps))
aa <- lapply(1:7, function(x) mymod$get_metapop_state( x));
# mymod$steps(nobs-365)
plot(
    levelplot(aa[[1]][1:100, ], aspect='fill')
)
