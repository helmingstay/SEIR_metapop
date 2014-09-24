

newModel <- function(initstates, transmat, accumvars, obsall, nobs, obs_nstep, deltat) {
    ## model construction function
    ## initstates:  states x population (integer) matrix of initial values of states
    ## transmat:    states x events (integer) matrix of update coefficients mapping events to states
    ## accumvar:    character vector naming events to accumulate
    ## obsall:      bool, observe states as well as accumvars
    ## nobs:        int, total number of observations
    ##
    ## errorchecking
    if ( !all( (initstates %% 1) == 0 )) stop("Only integer values allowed for initstates")
    if ( !all( (transmat %% 1) == 0 )) stop("Only integer values allowed for initstates")
    if ( !is.character(accumvars))  stop("Accumvars is a character vector naming events to accumulate")
    if ( !all( accumvars %in% colnames(transmat))) stop("accumvars must be valid event names (see colnames(transmat))") 
    if ( !identical( rownames(initstates), rownames(transmat))) {
        stop("Rownames (states) of initstates and transmat must match")
    }
    ##
    ## these names are important
    poplist <- list(accum=accumvars, obsall = obsall, nobs=nobs, obs_nstep=obs_nstep, deltat=deltat)
    npop <- ncol(initstates)
    model <- new(Metapop, npop, initstates, transmat, poplist)
    return(model)
}

newSEIRModel <- function(initstates, nsteps=30*365, accumvars=c("latent","imports", "infect"), obsall=T, obs_nstep=7, deltat=1/365) {
    ## model construction function
    ## initstates:  states x population (integer) matrix of initial values of states
    ## transmat:    states x events (integer) matrix of update coefficients mapping events to states
    ## accumvar:    character vector naming events to accumulate
    ## obsall:      bool, observe states as well as accumvars
    ## nobs:        int, total number of observations
    ##
    nobs <- 1+(nsteps/obs_nstep)
    eventnames <- c("dS", "latent", "imports", "infect", "recover", "dR");
    statenames <- c("S", "E", "I", "R", "N")
    transmat <- matrix(0, nrow=length(statenames), ncol=length(eventnames)); 
    rownames(transmat) <- rownames(initstates) <- statenames
    colnames(transmat) <- eventnames
    ## state transition definitions
    transmat["S", "dS"] = 1
    transmat["S", "latent"] = -1
    transmat["E", "latent"] = 1
    transmat["E", "imports"] = 1
    transmat["E", "infect"] = -1
    transmat["I", "infect"] = 1
    transmat["I", "recover"] = -1
    transmat["R", "recover"] = 1
    transmat["R", "dR"] = 1
    transmat["N", "dS"] = 1
    transmat["N", "dR"] = 1
    model <- newModel(initstates, transmat, accumvars, obsall, nobs, obs_nstep, deltat)
    return(model)
}


