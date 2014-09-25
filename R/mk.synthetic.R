#' function takes original timeseries and df of normal param estimates
#' and creates synthetic timeseries based on estimates
#' debug?? is this correct??
mk.synthetic <- function(myxts, distdf, prev.seq=10^seq(from=-8, to=0, by=0.1)) {
    mylen = nrow(myxts)
    nbreaks = length(prev.seq)
    ##
    ## vector of evenly spaced probabilities (0,1) to fill in, excluding 0 and 1
    myprobs = head(tail(seq(from=0, to=1, length.out=nbreaks),-1),-1)
    ## places to fill in with which probability
    myindex=quantile(1:mylen, probs=seq(from=0, to=1, length.out=nbreaks+1))
    myindex = cbind(head(myindex, -1), myindex[-1])
    ##
    ## modify data by city
    for (mycity in colnames(myxts))  {
        ## grab the estimated lognormal params for this city
        tmpdf = subset(distdf, placename==mycity, select=c(mean, sd))
        ## lognormal quantiles from probabilities and params
        ## this is prevalence
        mylevel =  10^qnorm(myprobs, mean=tmpdf$mean, sd=tmpdf$sd)
        ## grab pop and reprate
        mypop = subset(newdemog, placename==mycity)$pop
        myreprate = subset(newdemog,  placename ==mycity)$reprate
        ## use pop * per capita prevalence to imitates cases
        mycases = round(mypop*mylevel)#+.001
        ## formulation for negative binomial -- target successes, prob of success
        #myfill = rnbinom(length(mylevel), mycases*myreprate, myreprate)
        ## reporting as a binomial process
        myfill = rbinom(length(mylevel), mycases, myreprate)
        ## pull out the column to modify
        tmp = myxts[,mycity]
        for ( ii in 1:nbreaks ) {
            ## get the index, fill in
            thisindex = myindex[ii,1]:myindex[ii,2]
            tmp[thisindex] = myfill[ii]
        }
        myxts[,mycity] <- tmp
    }
    myxts
}
