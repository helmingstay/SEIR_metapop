rcpp_hello_world <- function(){
	.Call( "rcpp_hello_world", PACKAGE = "SEIR" )
}


newModel <- function(metapars, initstates, transmat, accumvars, obsall, nobs, obs_nstep, deltat) {
    ## model construction function
    ## metapars:    parameter list for the metapopulation (prestep, poststep)
    ## initstates:  states x population (integer) matrix of initial values of states
    ## transmat:    states x events (integer) matrix of update coefficients mapping events to states
    ## accumvar:    character vector naming events to accumulate
    ## obsall:      bool, observe states as well as accumvars
    ## nobs:        int, total number of observations

    ## errorchecking
    if ( !all( (initstates%1) == 0 )) stop("Only integer values allowed for initstates")
    if ( !all( (transmat%1) == 0 )) stop("Only integer values allowed for initstates")
    if ( !is.character(accumvars))  stop("Accumvars is a character vector naming events to accumulate")
    if ( !all( accumvars %in% colnames(transmat))) stop("accumvars must be valid event names (see colnames(transmat))") 
    if ( !all.equal( rownames(initstates), rownames(transmat))) {
        stop("Rownames (states) of initstates and transmat must match")
    }

    ## these names are important
    poplist <- list(accum=accumvars, obsall = obsall, nobs=nobs, obs_nstep=obs_nstep, deltat=deltat)

    npop <- ncol(initstates)
    model <- new(Metapop, npop, metapars, initstates, transmat, poplist)
    return(model)
}


## test code


## 10 pops, 4 states
initstates = =matrix(10, nrow=4, ncol=10)
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

accumvars = c("latent")
obsall = TRUE
nobs = 365*50
obs_nstep = 7
deltat = 1/365

mymod <- newModel(metapars, initstates, transmat, accumvars, obsall, nobs, obs_nstep, deltat)
