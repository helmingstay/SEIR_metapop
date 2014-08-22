require(foreach)
require(multicore)
require(doMC)
## required but implicit:
## load('full.RData')
#load('dat.demog.cityrates.RData')
#load('dat.city.pop.RData')
#load('dat.reprates.RData')
## force reload of SEIR
require(SEIR)
#detach(package:SEIR, unload=TRUE)
#require(SEIR)
source('mk.paperfuns.R')
source('mk.scalefuns.R')
source('mk.newfpc.R')
source('mk.ensemble.fit.R')


mk.modrun = function(pardf, index.limits, newdemog, 
                    city.rates, spec.names, obs_nstep=7,
                    .parallel=F, debug=F, nreps=1,
                    do.fpc=T, do.ratios=F, do.obs=F, do.spec=F){
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
    if (!is.ordered(newdemog$placename)) stop( 'Placename must be ordered by population size')
    newdemog <- newdemog[order(newdemog$placename),]
    city.names = levels(newdemog$placename)
    reprates = subset(newdemog, select=c(placename,reprate))
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
    #phi = 1e-4
    #diag(connectmat) = 1
    ## these 2 lines should be the same...
    #N0 = mk.pop(lcasedata$Measles$city.order, ldemog$city.dec, .year=1930, do.cum=F)
    ##
    ##
    ####
    ## invariant over all models, needed for initialization
    ##
    ## pull out demographics only for these cities
    city.rates = subset(city.rates, placename %in% city.names)
    N0 <- subset(city.rates, year==run.limits[1])$pop
    ncity = length(N0)
    ##dummy geographic matrix
    ## remove this from model code!!
    geogmat = matrix(rnorm(ncity*2), ncol=2) 
    with(pardf, if(R0<2) stop('Equilib calculated from R0, R0<2 not stable'))
    ## parameter list -- ok if there's extra stuff in it    
    ## params can change, one list per city
    par.template = as.list(pardf)
    within(par.template, {
        ## hardcoded initial births to calc equilib S/E/I/R to start with
        dS=0.02/365
        ## for model initialization
        dR = 0.0
        ##
        ## incubation
        sigma = 1/8
        ## infectious period
        gamma = 1/5

    }) -> par.template
    ##
    finret = list()
    finret$params=par.template
    finret$connectmat=connectmat
    ## dummy list
    ## pull out reporting rates for each city and add to 
    ## paramater list (one sublist for each city),
    ## starts the same for each city, 
    ## with the exception of pobs, which stays the same for the full sim
    ##
    parlist = llply(1:ncity, function(x) {
        ret = par.template; 
        ret$pobs = newdemog$reprate[x]
        ## adjust import rate by connection
        #ret$imports = with(newdemog[x,], connect/(ret$imports*pop)^ret$pop.pow)
        schooldays = sum(schoolterms[[x]]==1)
        ret$beta0 <- with( ret, 
            (R0*gamma)/(1+betaforce)^((2*schooldays-365)/365)
        )
        #beta0 = (rR0*rgamma) / pow(1.0+rbetaforce, (2.0*rschooldays - 365.0)/365.0);
        ##// See keeling et al 2001, "seasonally forced disease dynamics",
        ##// physica D, eqn 3.
        ## things we tried
        #ret$imports = ret$imports*newdemog$connect[x]
        #ret$imports = with(newdemog[x,], (connect*ret$imports)/sqrt(pop))
        ## model is way too high for the below
        #ret$imports = with(newdemog[x,], (1)/sqrt(pop))  
        if( ret$imports <0 ) {print("Negative imports!"); browser()}
        #if( ret$imports >1e-3 ) stop('imports very high!!')
        return(ret)
    })
    ##
    #cat(sprintf('starting %s sims\n', nreps))
    #cat(paste(rep('.', nreps), collapse=''))
    #cat('\n')
    llply(1:nreps, function(nrep) {
        #if( progress != 'none') cat('.')
        ##
        # not needed -- do in loop??
        # mymod$setpars(parlist)
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
        mymod = newSEIRModel(initmat, nsteps=nyears*365, obs_nstep=obs_nstep )
        mymod$set_school(schoolterms)
        ## add metapop pars here
        mymod$set_metapop(list(dummy=1))

        ##
        ## this is the order in helen's demog file...
        #ii = order(subset(demog.cityrates.us, year==1910)$pop)
        #cat('stepping through years')
        ## for each year in the run.limits
        ## run the actual model
        for (yy in seq(from=run.limits[1], to=run.limits[2], by=1)) {
        ## !!fixme??
        ## conflict on ocimum from head
        ## for (yy in seq(from=run.limits[1], to=run.limits[2]-1, by=1)) {
            ## update paramlist for this year
            ## get actual dS/dR values and set them
            ## get birthrate to set equilibrium?? 
            mydemog = subset(city.rates, year==yy)
            parlist = llply(1:ncity, function(ncity) {
                    ret = parlist[[ncity]]
                    ## should move this outside the loop...
                    ret$dS = mydemog[ncity,'dS']/365
                    ret$dR = mydemog[ncity,'dR']/365
                    return(ret)
                    }
            )
            mymod$set_pop(parlist)
            mymod$steps(365);
            #if(yy %%10 ==0 )browser()
        }
        ## get the (almost) final modeled population size 
        ## e.g. a year before finish
        ## and compare it with actual
        ## we have pop data for beginning of model year, so trim a year off end
        ## and take last row
        modpop <- tail(head(mymod$get_metapop_state(i.N),-52),1)
        mod.popsize <-  data.frame(placename=city.names, modpop=as.vector(modpop))
        data.popsize = subset(city.rates, year==run.limits[2], 
                select=c(placename, pop))
        fin.popsize = merge(data.popsize, mod.popsize) 
        fin.popsize$diff = with(fin.popsize, pop-modpop)
        fin.popsize$percent_diff = with(fin.popsize, diff/pop)
        if( with(fin.popsize, max(abs( percent_diff)) >0.10)) { 
        #if( with(fin.popsize, abs(median( percent_diff))  >0.01)) { 
            print(subset(fin.popsize, abs( percent_diff)  >0.10)) 
            print(pardf)
            print('Final model populations differ by more than 5% from truth')
            browser()
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
        
        ## check to make sure nothing's negative and we've dont the math right
        #if ( sum(rescue.mat< 0, na.rm=T) > 0 ) stop("Negative values in rescue mat (Eimport/(Estep-Eimport)), something's wrong")
        
        mk.xts <- function(mymat, .nstep=obs_nstep, .index.limits=index.limits, my.names=city.names) { 
            ## convenience function
            ## turn each matrix into a timeseries and then subset by the limits
            ## model outputs weekly, as shown by nstep (number of days per step)
            my.index <- as.Date( (0:(nrow(mymat)-1)*.nstep+1), origin=index.limits[1]);
            ## trim to match data time window
            ret <- xts(mymat, my.index)[index.limits,]
            colnames(ret) <- my.names
            ret
        }
        ## full output matrix as xts
        modrun.xts = mk.xts(Iobs.mat)
        ## number of weeks after subsetting
        n.cut.weeks <- nrow(modrun.xts)
        
        ##
        ## if debugging mode, plot shit
        ## this is sloooooow
        if (!identical(debug.info, F)){
            with(debug.info, {
                which.cities=seq(from=1, to=81, by=4)
                which.layout=c(3,7)
                this.modname = sprintf('imports_%s,betaforce_%s,R0_%s', imports, betaforce,  R0)
                pdf(file=sprintf('modfigs/xts-%s.pdf', this.modname), width=40, heigh=50 )
                ## be careful, cases.obs is external
                cat('\n\n\n#############################\nNaively using cases.obs to plot cases... are you sure this is the right object??\n\n')
                print( xyplot(cases.obs[, which.cities], type=c('l','h','g'), scales=list(x=list(alternating=F)), layout=which.layout, main=this.modname, 
                            xlab='Actual cases'), 
                        split=c(1,1,1,5), more=T)
                print( xyplot(modrun.xts[, which.cities], type=c('l','h','g'), scales=list(x=list(alternating=F)), layout=which.layout, 
                            xlab='Model observed cases'), 
                        split=c(1,2,1,5), more=T)
                print( xyplot(xts(is.infinite(obsratio.xts), index(obsratio.xts))[, which.cities], type=c('l','h','g'), scales=list(x=list(alternating=F)), 
                            xlab='Apparent extinctions', layout=which.layout), 
                        split=c(1,3,1,5), more=T)
                print( xyplot(xts(is.infinite(rescue.xts), index(rescue.xts))[, which.cities], type=c('l','h','g'), scales=list(x=list(alternating=F)), 
                            xlab='Rescues', layout=which.layout), 
                        split=c(1,4,1,5), more=T)
                print( xyplot(xts(is.infinite(rescue.xts) & is.infinite(obsratio.xts), index(rescue.xts))[, which.cities], type=c('l','h','g'), 
                            scales=list(x=list(alternating=F)), xlab='Apparent extinctions and rescues co-occur', layout=which.layout), 
                        split=c(1,5,1,5), more=F)
                dev.off()
            })
        }
        return(list(ratio.df=ratio.df, cases.obs=modrun.xts))
    }, .parallel=.parallel) -> simlist 
#, .progress='text')
    #cat('\n')
    if(do.fpc) {
        ldply(1:nreps, function(nrep) {
            ## make the fpc from the model timeseries
            cases.obs <- simlist[[nrep]][['cases.obs']]
            fpc <- mk.alldist(cases.obs, newdemog, .parallel=.parallel)
            fpc$nrep <- nrep
            fpc
        }, .parallel=.parallel) -> finret$fpc 
    } else { finret$fpc <- NULL}
    if(do.spec) {
        cases.list <- llply(1:nreps, function(nrep) {
            ## pull out all the lists of timeseries
            cases <- simlist[[nrep]][['cases.obs']]
        }) 
        finret$spec <- mk.spec.model(cases.list, spec.names)
    } else { finret$spec <- NULL}
    if(do.ratios) {
        ldply(1:nreps, function(nrep) {
            ## pull out all the ratiodfs, let plyr join them
            ratio <- simlist[[nrep]][['ratio.df']]
            ratio$nrep <- nrep
            ratio
        }) -> finret$ratio.df 
    } else { finret$ratio.df <- NULL}
    if(do.obs) {
        llply(1:nreps, function(nrep) {
            ## pull out all the lists of timeseries
            cases <- simlist[[nrep]][['cases.obs']]
        }) -> finret$cases.list
    } else { finret$cases.list <- NULL}
    return(finret)
}
 


mk.addcols.df = function(mydf, newdf, value=F) {
    ## add all the columns in newdf that are missing to mydf with values from newdf, unless value is specified
    ## NOTE only the first row of newdf is used
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





mk.label.fpc <- function(myfpc) { 
    ## label numeric factors of fpc (see SEIR C++ model for final definition of numeric codes)
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

mk.fpc.avg = function(myfpc) {
        ## should this go in mk.modrun??
        ## pull out the fpc for this sim
        ret = myfpc
        ## for each fpc, average over all repetitions, put in nrep=0
        ## myavg is temp df to write into
        myavg = subset(ret, nrep==1)
        myavg$nrep=-1
        ## at the end, discard any props that haven't changed (e.g. no samples)
        myavg$prop=-1
        for( myplace in unique(ret$placename)) {
            for (mycase in unique(ret$cases)) {
                ## place into here
                myoutrow = with(myavg, which(placename==myplace & cases==mycase & nrep==0))
                ## average over these -- exclude 0 & 1
                myrows = with(ret, placename==myplace & cases==mycase & between(prop,0,1))
                if (length(myoutrow) == 0 || sum(myrows) == 0) next
                myavg[myoutrow,'prop'] = mean(ret[myrows,'prop'])
                myavg[myoutrow,'X1'] = sum(myrows)
            }
        }
        rbind(ret, subset(myavg, prop>=0))
}


## measures

## hardcoded cols, should we fix?
mk.more.breakpoints = function(mydf, aa=-6, bb=-3) {
    ## for each level-comb of mydf, run mk.breakpoints, combine into data.frame
    ret = ddply( mydf, c('rep','phi', 'betaforce', 'R0', 'cum'), function(xx) {
        tmp <- nlm(mk.breakpoints, c(-6, -3), mydf = xx)
        #mk.breakpoints( tmp$estimate, xx, .plot=T)
        #readline(xx[1,])
        ret = as.list(c(tmp$estimate, tmp$minimum))
        names(ret) = c('Lower','Upper','RSS')
        return(data.frame(ret))
    }) #, .progress='text')
    return(ret)
}
