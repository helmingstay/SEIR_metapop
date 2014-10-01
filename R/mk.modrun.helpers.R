#' convenience function
#' turn each matrix into a timeseries and then subset by the limits
mk.mat.to.xts <- function(mymat, nstep, index.limits, my.names) { 
    ## model outputs weekly, as shown by nstep (number of days per step)
    my.index <- as.Date( (0:(nrow(mymat)-1)*nstep+1), origin=index.limits[1]);
    ## trim to match data time window
    ret <- xts(mymat, my.index)[index.limits,]
    colnames(ret) <- my.names
    ret
}

#' Error checking function
#' Compare final population size between model and demog
mk.modrun.check.popsize <- function(modpop, city.names, 
    city.rates, run.limits, thresh=0.1
) {
    mod.popsize <-  data.frame(placename=city.names, modpop=as.vector(modpop))
    data.popsize = subset(city.rates, year==run.limits[2], 
            select=c(placename, pop))
    fin.popsize = merge(data.popsize, mod.popsize) 
    fin.popsize$diff = with(fin.popsize, pop-modpop)
    fin.popsize$percent_diff = with(fin.popsize, diff/pop)
    .over <- subset(fin.popsize, abs( percent_diff)  >thresh) 
    if( nrow(.over) >0 ) { 
        warning(
            sprintf('Final model populations differ by more than %d%% from truth', thresh*100)
        )
        print(.over)
        print(pardf)
    }
}

#' if debugging mode, plot actual and modeled cases
#' Use list debug.info to get relevant vars
#' this is sloooooow
mk.modrun.diag.plots <- function(debug.info, cases.obs, modrun.xts, .plotdir = 'modfigs') {
    with(debug.info, { 
        ## hardcoded cities??
        which.cities <- seq(from=1, to=81, by=4)
        which.layout <- c(3,7)
        this.modname  <-  sprintf('imports_%s,betaforce_%s,R0_%s', imports, betaforce,  R0)
        pdf(file=sprintf('%s/xts-%s.pdf', .plotdir, this.modname), width=40, heigh=50 )
        ## on one page, plot timeseries of cases and possible states
        ## e.g. apparent extinctions, rescues
        ## 
        ## store plots in list
        .mk.plot <- function(.dat, .xlab, .layout, .main=NULL) {
            ret <- xyplot(
                .dat, type=c('l','h','g'), 
                scales=list(x=list(alternating=F)), 
                layout=.layout, main=.main, xlab=.xlab
            )
        }
        .p <- list()
        .p[[1]] <- .mk.plot(cases.obs[, which.cities], .layout=which.layout, main=this.modname, .xlab='Actual cases')
        .p[[2]] <- .mk.plot(modrun.xts[, which.cities], .layout=which.layout, .xlab='Model observed cases')
        .p[[3]] <- .mk.plot(xts(is.infinite(obsratio.xts), index(obsratio.xts))[, which.cities], .layout=which.layout, .xlab='Apparent extinctions')
        .p[[4]] <- .mk.plot(xts(is.infinite(rescue.xts), index(rescue.xts))[, which.cities], .layout=which.layout, .xlab='Rescues')
        .p[[5]] <- .mk.plot(xts(is.infinite(rescue.xts) & is.infinite(obsratio.xts), index(rescue.xts))[, which.cities], .layout=which.layout, .xlab='Apparent extinctions and rescues co-occur') 
        .np <- length(.p)
        for (ii in 1:.np) {
            print(.p[[ii]], split = c(1,ii,1,.np), more=(ii !=.np))
        }
        dev.off()
    })
}

#' for each simulation, 
#' infer distribution of log case reports
mk.modrun.distrib <- function(simlist) {
    ret <- ldply(simlist, .id='nrep', function(.rep) {
            ## make the fpc from the model timeseries
            cases.obs <- .rep[['cases.obs']]
            ## fit distribution of log case reports
            ret <- mk.all.distrib(cases.obs)
            ret
    })
    return(ret)
}

#' for each simulation, compute spectrum
mk.modrun.spec <- function(simlist, spec.names) {
    cases.list <- llply(simlist, function(.rep) {
        ## pull out all the lists of timeseries
        cases <- .rep[['cases.obs']]
    })
    ret <- mk.spec.model(cases.list, spec.names)
    ret
}

#' deprecated??
#' add all the columns in newdf that are missing to mydf with values from newdf, unless value is specified
#' NOTE only the first row of newdf is used
mk.addcols.df = function(mydf, newdf, value=F) {
    oldnames = colnames(mydf)
    newnames = colnames(newdf)
    ## all the ones in new that are *not* in old
    addnames = newnames[ !( newnames %in% oldnames) ]
    ## if no value is given, use those in newdf
    if (identical(value, FALSE)) {
        for (ii in addnames) {
                mydf[,ii] = newdf[1,ii]
        }
    } else {
        for (ii in addnames) {
                mydf[,ii] = value
        }
    }
    mydf
}



#' deprecated??
#' keeping for reference
#' different import methods, 
#' label numeric factors of fpc (see SEIR C++ model for final definition of numeric codes)
mk.label.fpc <- function(myfpc) { 
   within(myfpc, {
        if( is.numeric( importmethod)){
            importmethod <- factor(importmethod, levels=0:4, labels=c('No imports', 'Constant imports', 'Pulsed imports', 
                                                            'Constant pop-proport imports', 'Pulsed pop-proport imports'))
        }
        if( is.numeric( schooltype)){
            schooltype <- factor(schooltype, levels=0:1, labels=c('Sin forcing', 'Term forcing'))
        }
        if(is.numeric( betaforce)  ) {
            betaforce <- factor(betaforce, levels= unique(betaforce), labels = sprintf("betaforce: %s", unique(betaforce)))
        }
        if(is.numeric( R0)  ) {
            R0 <- factor(R0, levels= unique(R0), labels = sprintf("R0: %s", unique(R0)))
        }
        ## check to see if myhigh is in colnames first...
        #if(is.numeric( myhigh)  & length(unique(myhigh)) ) {
            #myhigh <- factor(myhigh, levels= unique(myhigh), labels = sprintf("high: %2.2e",   unique(myhigh)))
        #}
    }) -> myfpc
    myfpc
}
