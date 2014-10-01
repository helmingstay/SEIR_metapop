source(
    system.file('examples', 'run.prepmod.R', package='SEIR')
)
##
## setup multicore
library(doMC)
registerDoMC()
## uncomment to use more than 1 cores
#registerDoMC(cores=4)
## run simulation in parallel in inner loop?
in.parallel <- F
## number of simulations, should be >1
nreps <- 2

##  then modify mod.final as a copy of mod.prep
## final simulation parameter specifications
mod.final <- within(mod.prep, {
    ## best parameter values
    US$model.pars <- within(US$model.pars, {
        R0=18  ## R0 must be present in demog, city.rates dfs
        probs=0 
        distmethod=0 ## deprecated??
        betaforce=0.5  ## seasonal forcing
        pop.pow=0 ## deprecated??
        importmethod=c(-1) 
        imports=2e-8 ## per day probability of E import??
        ## schooltype 0 = sin, type 1 = term
        schooltype=c(0)
    })
    `England & Wales`$model.pars <- within(`England & Wales`$model.pars, {
        R0=24
        probs=0 
        ## deprecated?
        distmethod=0
        betaforce=0.5
        pop.pow=0 ## depre
        importmethod=c(-1)
        imports=5e-9
        ## schooltype 0 = sin, type 1 = term
        schooltype=c(1) 
    })
})

##  Get list of vars from US, from run.prepmod.R
## run simulations
ensemble <- llply( mod.final, function(thislist) {
    cat(sprintf("Running simulation with %s reps\n", nreps))
    thislist <- within(thislist, {
        .R0 <- model.pars$R0
        demog <- subset(demog, R0==.R0)
        city.rates <- subset(city.rates, R0==.R0)
    })
    thislist$sims <- with(thislist, {
        ncity <- length(city.order)
        demog$connect <- 1
        ## connection matrix not currently used, 
        ## required by simulation
        dummy.mat = matrix( 1, nrow=ncity, ncol=ncity)
        mk.modrun(model.pars, index.limits, demog, 
            city.rates, cases.obs,
            colnames(thislist$spec),
            .parallel= in.parallel, nreps=nreps, debug=F, 
            do.obs=T, do.spec=F, do.distrib=T, do.ratios=T
        )
    })
    thislist
})

## pull out estimated distributions
ensemble.distrib <- ldply(ensemble, .id='where', function(x) x$sims$distrib)

## summarize simulation results for each country
mod.distrib.df <- ddply(ensemble.distrib,  'placename', function(x){
        ## assumes england and wales doesn't share placenames with us
        ## pull out one row, compute means, place into ret
        ret <- x[1,]
        myvars = c('mean','sd','error') 
        ret[1, myvars] = colMeans( x[, myvars])
        ret
})
mod.distrib.df$model <- factor('Model')

## writed out the results to an object for safekeeping
#save(ensemble, moddist.df, ensemble.distrib, file='cache.ensemble.RData')
## do more processing in run.combine.R

## compare??  use lm to 
#ensemble.fit <- mk.ensemble.fit.all(ensemble.dist, mod.prep) 
#ensemble.xtable <- mk.ensemble.fit.all(ensemble.dist, mod.prep, do.xtable=T) 
#ensemble.xtable <- colname.replace('where', 'Country', ensemble.xtable )
#ensemble.xtable <- subset(ensemble.xtable, select=c(-R0, -betaforce,-imports))
#print(xtable(ensemble.xtable, digits=3), include.rownames=F, file='figs/lmtable.tex')

