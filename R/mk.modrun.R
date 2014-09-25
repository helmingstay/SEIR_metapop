#' Main wrapper function to run stochastic metapop sim
mk.modrun = function(parlist, index.limits, demog, 
                    city.rates, cases.obs,
                    spec.names, obs_nstep=7,
                    .parallel=F, debug=F, nreps=1,
                    do.dist=T, do.ratios=F, do.obs=F, 
                    do.spec=F, 
                    check.popsize=T, quiet=F)
{
    ## index constants
    i.Estep <- 1
    i.Eimport <- 2
    i.Istep <- 3
    i.N <- 8
    ## set-up
    if(debug) {debug.info <- moddf} else { debug.info <- FALSE }
    ##
    ## range to run the model == time we have demographics for
    run.limits = range(city.rates$year)
    ##
    ## mk.modrun assumes a particular order
    if (!is.ordered(demog$placename)) stop( 'Placename must be ordered by population size')
    demog <- demog[order(demog$placename),]
    city.names = levels(demog$placename)
    reprates = subset(demog, select=c(placename,reprate))
    ## get population-ordered reprates
    reprates = with(reprates, reprate[ order(placename)])
    ##
    ncity <- length(city.names)
    ## still required by model, but not used
    connectmat <- matrix( 1, nrow=ncity, ncol=ncity)
      
    schoolterms <- lapply(1:ncity, function(x){
        ## Termtime Ref: Seasonally forced disease dynamics explored as switching between attractors
        ## Matt J. Keeling∗ , Pejman Rohani, Bryan T. Grenfell
        ## Physica D 2001
        schooldays = rep(1, 365)  ## 1 == term
        schooldays[c(xmas=1:6, easter=100:115, summer=200:251, fall=300:307, xmas1=356:365)] = -1  ## vacation
        ## US? Ref??
        #schooldays = rep(-1, 365)  ## -1 == vacation
        #schooldays[c(fall1=240:330, fall2=334:350, spring1=5:90, spring2=100:150)] = 1  ## school term
        ## Rohani, opposite patterns of synchrony
        ##
        schooldays
    })
    ## schooltype: 
    ## 0 == sin, 1 == step (calculate beta low/hi from R0 & schoolyear)
    #schooltype = 1
    #N0 = 1e6
    ## connection matrix?
    #phi = 1e-4
    #diag(connectmat) = 1
    ####
    ## invariant over all models, needed for initialization
    ## pull out demographics only for these cities
    city.rates = subset(city.rates, placename %in% city.names)
    N0 <- subset(city.rates, year==run.limits[1])$pop
    ncity = length(N0)
    ##dummy geographic matrix
    ## remove geogmat from model code??
    geogmat = matrix(rnorm(ncity*2), ncol=2) 
    with(parlist, if(R0<2) stop('Equilib calculated from R0, R0<2 not stable'))
    ## parameter list -- ok if there's extra stuff in it    
    ## params can change, one list per city
    par.template <- parlist
    ## return list
    finret <- list()
    ## copy model specs into return list
    finret$params <- par.template
    finret$connectmat <- connectmat
    ## dummy list
    ## pull out reporting rates for each city and add to 
    ## paramater list (one sublist for each city),
    ## starts the same for each city, 
    ## with the exception of pobs, which stays the same for the full sim
    ##
    city.parlist = llply(1:ncity, function(x) {
        ret <- par.template; 
        ret$pobs <- demog$reprate[x]
        ## adjust import rate by connection
        #ret$imports = with(demog[x,], connect/(ret$imports*pop)^ret$pop.pow)
        schooldays <- sum(schoolterms[[x]]==1)
        ret$beta0 <- with( ret, 
            (R0*gamma)/(1+betaforce)^((2*schooldays-365)/365)
        )
        #beta0 = (rR0*rgamma) / pow(1.0+rbetaforce, (2.0*rschooldays - 365.0)/365.0);
        ##// See keeling et al 2001, "seasonally forced disease dynamics",
        ##// physica D, eqn 3.
        ## things we tried
        #ret$imports = ret$imports*demog$connect[x]
        #ret$imports = with(demog[x,], (connect*ret$imports)/sqrt(pop))
        ## model is way too high for the below
        #ret$imports = with(demog[x,], (1)/sqrt(pop))  
        if( ret$imports <0 ) {warning("Negative imports!"); browser()}
        return(ret)
    })
    ##
    #cat(sprintf('starting %s sims\n', nreps))
    #cat(paste(rep('.', nreps), collapse=''))
    #cat('\n')
    llply(1:nreps, function(nrep) {
        ##
        ## initial states as data.frame
        ## start at (approx) equilibrium -- same birth for everyone
        initdf = with(par.template, data.frame(S=N0/R0, E=(N0*dS)/sigma, I=(N0*dS)/gamma))
        initdf=round(initdf)
        ## compute R, N
        initdf = cbind(initdf, with(initdf, 
                        data.frame(R = N0-rowSums(initdf), N=N0)
        ))
        ## turn to matrix for SEIR
        initmat = t(as.matrix(initdf))
        nyears = diff(run.limits)+1
        mymod = newSEIRModel(initmat, obs_nstep=obs_nstep, nsteps=nyears*365)
        mymod$set_school(schoolterms)
        ## add metapop pars here
        mymod$set_metapop(list(dummy=1))

        ##
        if (!quiet) cat(sprintf('\n Sim #%d', nrep))
        ## for each year in the run.limits
        ## run the actual model
        for (yy in seq(from=run.limits[1], to=run.limits[2], by=1)) {
        ## !!fixme??
        ## conflict on ocimum from head
        ## for (yy in seq(from=run.limits[1], to=run.limits[2]-1, by=1)) {
            ## update paramlist for this year
            ## get actual dS/dR values and set them
            ## get birthrate to set equilibrium?? 
            mydemog <- subset(city.rates, year==yy)
            city.parlist <- llply(1:ncity, function(ncity) {
                    ret <- city.parlist[[ncity]]
                    ## should move this outside the loop...
                    ret$dS <- mydemog[ncity,'dS']/365
                    ret$dR <- mydemog[ncity,'dR']/365
                    return(ret)
                    }
            )
            mymod$set_pop(city.parlist)
            mymod$steps(365);
            #if(yy %%10 ==0 )browser()
            if (!quiet) cat('.')
        }
        ## get the (almost) final modeled population size 
        ## e.g. a year before finish
        ## and compare it with actual
        ## we have pop data for beginning of model year, so trim a year off end
        ## and take last row
        if (check.popsize) {
            modpop <- tail(head(mymod$get_metapop_state(i.N),-52),1)
            mk.modrun.check.popsize(modpop, city.names, city.rates, run.limits)
        }
        ## preallocate and assign by column
        ## pull out weekly true and observed state variables as matrices
        Eimport.mat <-  mymod$get_metapop_state(i.Eimport)
        Estep.mat <-   mymod$get_metapop_state(i.Estep)
        Istep.mat <-   mymod$get_metapop_state(i.Istep)
        ## Total 
        EIstep.mat <- Estep.mat + Istep.mat
        .nrows = nrow(Istep.mat)
        Iobs.mat <-   do.call( 
            cbind, lapply(1:ncity, function(ii) {
                ret = rbinom(.nrows, Istep.mat[,ii], reprates[ii])
            })
        )
        #colnames(Eimport.mat ) <- colnames(Estep.mat ) <- colnames(EIstep.mat )
            #colnames(Istep.mat ) <- colnames(Iobs.mat ) <- city.names
        events <- within(list(),{
            `Rescue event` <- (EIstep.mat == Eimport.mat) & (Eimport.mat >0)
            `True extinction` <- (EIstep.mat == 0)
            `Apparent extinction` <- (EIstep.mat > 0) & (Iobs.mat == 0)
        })
        ratio.df <- ldply( events, function(x) {
            colnames(x) <- city.names
            ## get proportion of time by col(city)
            adply(x, 2, function(y) { prop=sum(y)/length(y)})
        })
        colnames(ratio.df) <- c('event', 'placename', 'value')
        ratio.df$event <- factor(ratio.df$event, levels = names(events))
        ## full output matrix as xts
        modrun.xts = mk.mat.to.xts(Iobs.mat, nstep=obs_nstep, index.limits=index.limits, my.names=city.names)
        ## number of weeks after subsetting
        n.cut.weeks <- nrow(modrun.xts)
        ## plot timeseries of cases & state transitions to pdfs
        if (!identical(debug.info, F)){
            mk.modrun.diag.plots(debug.info, cases.obs, modrun.xts)
        }
        return(list(ratio.df=ratio.df, cases.obs=modrun.xts))
    }, .parallel=.parallel) -> simlist 
    ## used for nrep .id column
    names(simlist) <- 1:length(simlist)
#, .progress='text')
    ## various reports, store results in return list
    if(do.dist) {
        finret$dist <- mk.modrun.dist(simlist)
    }
    if(do.spec) {
        finret$spec  <- mk.modrun.spec(simlist, spec.names)
    }
    if(do.ratios) {
        finret$ratio.df <- ldply(simlist, .id='nrep', function(.rep) {
            ## pull out all the ratiodfs, let plyr join them
            ret <- .rep[['ratio.df']]
        }) 
    } 
    if(do.obs) {
        finret$cases.list <- llply(simlist, function(.rep) {
            ## pull out all the lists of timeseries
            cases <- .rep[['cases.obs']]
        })
    }
    return(finret)
}
 
