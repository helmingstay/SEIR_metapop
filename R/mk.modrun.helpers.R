## TODO!! see R/mk.modrun.R for use example
newSEIRModel <- function(initmat, nsteps, obs_nstep){
    ret = new(Metapop, ncity, geogmat, connectmat, nobs, nstate, nstep)
    return(ret)
}


#' convenience function
#' turn each matrix into a timeseries and then subset by the limits
mk.modrun.mat.to.xts <- function(mymat, .nstep=obs_nstep, .index.limits=index.limits, my.names=city.names) { 
    ## model outputs weekly, as shown by nstep (number of days per step)
    my.index <- as.Date( (0:(nrow(mymat)-1)*.nstep+1), origin=index.limits[1]);
    ## trim to match data time window
    ret <- xts(mymat, my.index)[index.limits,]
    colnames(ret) <- my.names
    ret
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
