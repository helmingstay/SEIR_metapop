source('mk.modrun.R');
# in SOURCE.R
#source('run.prepmod.R')
##
registerDoMC(cores=4)

in.parallel <- T
nreps <- 10

## moved from run.finalensemble.R
## see:
##  prob1 = 0.01; prob2 = 0.90; (subset( resid.measure, value<quantile(value, prob1) &corr>quantile(corr, prob2) & where == 'US'))
mod.final <- mod.prep
mod.final$US$model <- data.frame( R0=18,  ## R0 must be present in demog, city.rates table
                    probs=0, distmethod=c('null'),
                    betaforce=0.5,
                    pop.pow=0,
                    imports=2e-8,
                    ## schooltype 0 = sin, type 1 = term
                    schooltype=c(0), importmethod=c(-1)
)

## see:
## prob1 = 0.01; prob2 = 0.7; (subset( resid.measure, value<quantile(value, prob1) &corr>quantile(corr, prob2) & where != 'US'))
mod.final$`England & Wales`$model <- data.frame( R0=24,
                    probs=0, distmethod=c('null'),
                    betaforce=0.5,
                    pop.pow=0,
                    imports=5e-9,
                    ## schooltype 0 = sin, type 1 = term
                    schooltype=c(1), importmethod=c(-1)
)

##
##  Get list of vars from US, from run.prepmod.R
ensemble <- llply( mod.final, function(thislist) {
    cat(sprintf("Running model with %s reps\n", nreps)) 
    thislist <- within(thislist, {
        .R0 <- model$R0
        newdemog <- subset(newdemog, R0==.R0)
        alldist.df <- subset(alldist.df, R0==.R0)
        city.rates <- subset(city.rates, R0==.R0)
        prev <- mk.prev(cases.obs, newdemog, zero.denom=Inf)
    })
    thislist$sims <- with(thislist, {
        ncity <- length(city.order)
        newdemog$connect <- 1
        ## still required by model, but not used
        dummy.mat = matrix( 1, nrow=ncity, ncol=ncity)
        mk.modrun(model, index.limits, newdemog, 
            city.rates, colnames(thislist$spec),
            .parallel= in.parallel, nreps=nreps, debug=F, 
            do.obs=T, do.spec=F, do.fpc=T, do.ratios=T
        )
    })
    ## measure diff between model ensemble and data
    #thislist$mod.diff <- mk.measure.mod(thislist)
    
    ## add the model results back to the list and return everything together
    #for (ii in names(ret)) {
        #thislist[[ii]] <- ret[[ii]]
    #}
    thislist
})

ensemble.fpc <- ldply(ensemble, function(x) x$sims$fpc)
ensemble.fpc <- colname.replace('.id',  'where', ensemble.fpc)
ensemble.fpc <- wherefix(ensemble.fpc)

moddist.df <- ddply(ensemble.fpc,  'placename', function(x){
        ## assumes england and wales doesn't share placenames with us
        ## pull out one row, compute means, place into ret
        ret <- x[1,]
        myvars = c('mean','sd','CV','error') 
        ret[1, myvars] = colMeans( x[, myvars])
        ret
})
moddist.df$model <- 'Model'
moddist.df$where <- sub('EW', 'England & Wales', moddist.df$where)



## compare??  use lm to 
#ensemble.fit <- mk.ensemble.fit.all(ensemble.fpc, mod.prep) 
#ensemble.xtable <- mk.ensemble.fit.all(ensemble.fpc, mod.prep, do.xtable=T) 
#ensemble.xtable <- colname.replace('where', 'Country', ensemble.xtable )
#ensemble.xtable <- subset(ensemble.xtable, select=c(-R0, -betaforce,-imports))
#print(xtable(ensemble.xtable, digits=3), include.rownames=F, file='figs/lmtable.tex')

## writed out the results to an object for safekeeping
save(ensemble, moddist.df, ensemble.fpc, file='cache.ensemble.RData')
## do more processing in run.combine.R
