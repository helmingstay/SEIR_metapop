#' convenience function, find total zeros
sum0 <- function(x, ...) sum(x==0, ...)

#' for each timeseries in case.list
#' apply function to each col 
#' bind together and aggregate results with aggfun
mk.bycol.bind <- function( cases.list, .fun=sum0, .na.rm=T){
    ret <- lapply(cases.list, function(x) { 
                ## column / city summaries
                apply( x, 2, .fun, na.rm=.na.rm)
    })
    ret <- do.call( rbind, ret)
    return(ret)
}

#' For each timeseries in list, sum each column over period weeks
#' Return modified list of timeseries
mk.resample <- function( cases.list, period=2) {
    ## set up index using first timeseries in list 
    index = c(0, seq(from=2, to=nrow(cases.list[[1]]), by=period))
    ret <- llply(cases.list, function(.cases) { 
        ## for each column, do period.apply, return list
        new.cases = alply(.cases, 2, function(x) period.apply(x, index, sum)); 
        ## combine list into timeseries
        new.cases = do.call(cbind, new.cases)
        ## clean up and return
        colnames(new.cases) <- colnames(.cases)
        new.cases
    }, .progress='none', .parallel=T)
    ret
}

mk.city.stats <- function(cases.list, .period=1) {
    if (.period !=1) {
        cases.list  <- mk.resample(cases.list, .period)
    }
    ## for each model run
    ## get aggregate prop zeros and max for each city
    .zeros <- mk.bycol.bind(cases.list, .fun=sum0)
    ## get prop zeros
    .zeros <- .zeros/nrow(cases.list[[1]])
    .zeros <- melt( .zeros, varnames=c('nrep','placename'))
    .zeros$variable = 'zeros'
    ## need to divide by pop and reprate later to get percap incidence
    .max <- mk.bycol.bind(cases.list, .fun=max)
    .max <- melt( .max, varnames=c('nrep','placename'))
    .max$variable = 'max'
    ## combine results, cleanup & return
    ret <- rbind(.zeros, .max)
    ret$variable <- factor(ret$variable)
    return(ret)
}
