## test code
if (T) {
## 10 pops, 4 states
initstates <- matrix(10, nrow=4, ncol=10)
rownames(initstates) <- c("S", "E", "I", "R") 
transmat <- matrix(0, nrow=4, ncol=5); 
colnames(transmat) <- c("birth", "latent", "infect", "recover", "deltaR"); 
rownames(transmat) <- c("S", "E", "I", "R")
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
}
