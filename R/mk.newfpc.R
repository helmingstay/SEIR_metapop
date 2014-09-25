## these are be deprecated
## in favor of functions in repo extinct
if(F) {
mk.prev <- function(cases.obs, newdemog, inferred=T, prev=T, rm.zeros=T, 
                    zero.denom=2, scale.prev=1) {
    ##
    ## takes a timeseries matrix of observed cases
    ## first exclude NAs
    ## optionally, compute inferred cases from reporting rate
    ## optionally, divide by pop to get inferred prev
    ## 
    ## keep track of whether it's zero
    ## if so use 1/N for 0 instead
    init.zero.denom <- zero.denom
    ## return as list since variable length
    cases.inf <- lapply(colnames(cases.obs), function(mycity) {
            ## pull out this column of data, and reporting rate
            ## omitting nas
            myxts <- coredata(na.omit(cases.obs[,mycity]))
            myrep = subset(newdemog, placename==mycity)$reprate
            ## use zero.denom == 0 as indication to use 1/N
            if (init.zero.denom==0) {zero.denom = 1/myrep}
            ## index for zero and nonzero data
            ii <- myxts==0
            jj <- myxts!=0
            if (inferred) {
                #### convert timeseries to inferred cases
                ## for zeros, set equal to 1/zero.denom reprate, 
                ## which gives a result between zero and 1 observed case reports adjusted by reprate
                ## should we round this??
                myxts[ii] <- 1/(zero.denom*myrep)
                ## otherwise, scale by reprate
                myxts[jj] <- myxts[jj]/myrep
            }
            if (prev) {
                ## divide  cases (inferred or raw) by pop to get prev
                mypop <- subset(newdemog, placename==mycity)$pop
                myxts <- myxts/mypop
            }
            if(rm.zeros) {
                ## at the end if there are any zeros (e.g. zero.demom =  Inf)
                ## should we remove them? 
                ## we can't take mean of log10(0), 
                ## and we get a CDF step function of sum(0 <= thresh.prev)
                ## for arbitrarily small thresh.prev
                ## so this should only be F for testing
                myxts <- myxts[myxts>0]
            }
            ## optionally, rescale into units of incidence per week
            ## e.g. divide by 4.35 for monthly cases
            myxts <- myxts/scale.prev
            return(as.vector(myxts))
    })
    names(cases.inf) <- colnames(cases.obs)
    ## return as a dataframe?
    ## compute ecdf and bind to prev
    ret <- lapply(names(cases.inf), function(y) { x = cases.inf[[y]]; 
                    data.frame(placename=y, prev=x, prop=ecdf(x)(x))})
    ret <- do.call(rbind, ret)
    #colnames(ret)[1] <- 'placename'
    ## reorder by input col, which should be city size order
    ret$placename <- ordered( ret$placename, levels=colnames(cases.obs))
    return(ret)
}

mk.sampledist <- function(prev.list, level, newdemog, do.log=T) {
    ## take prevalence results of mk.prev, return mean
    ## either by pop or for whole pop together.
    mydf.fun = function(x) {
            ## exclude zero from everything else.  
            ## zero prev only seen when zero.denom==Inf
            if(do.log) {
                myprev = log10(x)
            } else {
                ## we're dealing with cases or something
                myprev = x
            }
            mymean=mean(myprev)
            mysd=sd(myprev)
            mycv <- mysd/mymean
            ## equals statistic from ks test
            #myerror <- ks.test(myprev, "pnorm", mymean, mysd,exact=F)$statistic
            myerror <- max(abs(ecdf(myprev)(myprev) - pnorm(myprev, mymean, mysd)))
            data.frame(mean=mymean, sd=mysd, CV=mycv, error=myerror)
    }
    switch(level,
        city = ldply(prev.list, mydf.fun),
        all =   cbind(.id='all', mydf.fun(unlist(prev.list))),
        stop(sprintf('Invalid level %s, must be city or all', level))
    ) -> ret
    colnames(ret) <- sub('.id','placename', colnames(ret))
    ## add pop back in, only meaninful for city
    ret <- merge(newdemog, ret, all.x=F, all.y=T)
    return( cbind(ret, method='Sample', level=level))
}

## convenience function for mk.fpcdist 
## to do the estimation and construct appropriate data.frame
mk.estdf <- function(mydf, init=c(-4.5,1), norm=function(x) max(x)) {
    ## use ddply to do by city
    ## empirical prev
    myprev <- log10(mydf$prev)
    ## from ecdf
    myprop <- mydf$prop
    est = nlm( mk.fitdist, init, x=myprev, y=myprop, norm=norm)$estimate
    mymean <- est[1]
    mysd <- est[2]
    mycv <- mysd/mymean
    ## compute max diff between theoretical (pnorm) and observed probabilities
    ## e.g. KS diff,
    ## identical to minimum above if norm=max
    myerror <- max(abs(myprop - pnorm(myprev, mymean, mysd)))
    #myerror <- max( abs(pnorm(myprev, mymean, mysd) - myprop))
    est.cdf <- data.frame(mean=mymean, sd=mysd, CV=mycv, error=myerror, method='CDF')
}

mk.fpcdist <- function(prop.df, newdemog, init, norm, .parallel=TRUE) {
    ## take the prevalence dataframe and compute the per-city or total normal mean and variance
    ## using nlm/norm to minimize distance between empirical proportion and theoretical 
    ##apply the mk.estdf to either each city individuall, or all together
    ret <- ddply(prop.df, 'placename', mk.estdf, init=init, norm=norm)
    ## change name back after ddply
    #colnames(ret)[1] <- 'placename'
    ## add pop back in
    ret <- merge(newdemog, ret, all.x=F, all.y=T, by='placename')
    ## reorder
    ret
}


#' Infer each city's distribution of inference 
#' From case reports, population, reporting rate
mk.alldist <- function(cases.obs, newdemog, zero.denom=2, 
                    init=c(-4,2), norm=max) {
    ## scale.prev now taken from newdemog
    scale.prev <- unique(newdemog$period)
    if (length(scale.prev) != 1) stop( 'Check weeks in newdemog -- must be single value')
    ## make prevalence list 
    prev.df <- mk.prev(cases.obs, newdemog, zero.denom=zero.denom, 
                        scale.prev=scale.prev)
    ## fit a lognormal to each city
    myfpc <- mk.fpcdist(prev.df, newdemog, init=init, norm=norm)
    newinit <- with(myfpc, c( mean(mean), mean(sd)))
    ## take initial estimates and use them as starting values for another est
    myfpc <- mk.fpcdist(prev.df, newdemog, init=newinit, norm=norm)
}
}
