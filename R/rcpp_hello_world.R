rcpp_hello_world <- function(){
	.Call( "rcpp_hello_world", PACKAGE = "SEIR" )
}


newModel <- function(initstates, transmat, accumvars, obsall, nobs, obs_nstep, deltat) {
    ## model construction function
    ## initstates:  states x population (integer) matrix of initial values of states
    ## transmat:    states x events (integer) matrix of update coefficients mapping events to states
    ## accumvar:    character vector naming events to accumulate
    ## obsall:      bool, observe states as well as accumvars
    ## nobs:        int, total number of observations

    ## errorchecking
    if ( deltat < 0 )  stop("deltat must be >  0")
    if ( !all( (initstates %% 1) == 0 )) stop("Only integer values allowed for initstates")
    if ( !all( (transmat %% 1) == 0 )) stop("Only integer values allowed for initstates")
    if ( !is.character(accumvars))  stop("Accumvars is a character vector naming events to accumulate")
    if ( !all( accumvars %in% colnames(transmat))) stop("accumvars must be valid event names (see colnames(transmat))") 
    if ( !identical( rownames(initstates), rownames(transmat))) {
        stop("Rownames (states) of initstates and transmat must match")
    }

    ## these names are important
    poplist <- list(accum=accumvars, obsall = obsall, nobs=nobs, obs_nstep=obs_nstep, deltat=deltat)

    npop <- ncol(initstates)
    model <- new(Metapop, deltat, npop, initstates, transmat, poplist)
    return(model)
}


