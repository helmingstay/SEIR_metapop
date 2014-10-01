## prep for most of the run.*.R files:
## run.comine, run.newfpc, run.ensemble...
## prep objects for model run
if(F) {#!!todo - cleanup
source('mk.newfpc.R')
source('mk.util.R')
source('mk.spec.R')
source('mk.connect.R')
source('mk.cityrates.R')
source('mk.reprates.R')
source('mk.measuremod.R')
source('mk.range.R')
require(reshape)
load('full.RData')
require(plyr)
}

## prereqs - what's the correct way to do this in-package??
data(US.case.reports)
data(EW.case.reports)
data(US.demographics)
data(EW.demographics)
library(reshape2)
library(plyr)
library(xts)

## this list will be used to run simulations 
## one element each for US, EW
mod.prep = list()


## model parameters used for both coutries
mod.template <- list(
    #myR0s = 10:25
    #myR0s = 20
    ## this should include the final values for US, EW
    ## otherwise specify per-list
    myR0s = c(18,20,24),
    model.pars = list(
    ## hardcoded initial births to calc equilib S/E/I/R to start with
        dS=0.02/365,
        ## for model initialization
        dR = 0.0,
        ## incubation
        sigma = 1/8,
        ## infectious period
        gamma = 1/5
    )
)


## ex
mod.prep$US <- within(mod.template, {
    cases.obs <- US.case.reports$Measles$city.zoo["::1945-12",]
    ## weeks per sample
    period <- 1
    ################
    ##  Demographics
    ################
    city.rates <- US.demographics$city.rates
    city.ranks <- US.case.reports$Measles$city.ranks
    where <- 'US'
    ##
})

## same for England and Wales
mod.prep[['England & Wales']] <- within(mod.template, {
    ## to match names in US
    ## there are several administrative changes in 1966, which affects 1965 migration
    ## only use TS (and run model) through end of 1964
    cases.obs <- EW.case.reports["::1964",]
    ## weeks per sample
    period <- 2
    ################
    ##  Demographics
    ################
    city.rates <- EW.demographics$city.rates
    city.ranks <- EW.demographics$city.ranks
    ## London changes in 1966.  Reuse values from 1063 in 1964
    tmp.rates <- subset(city.rates, placename=='London' & year==1963)
    tmp.rates$year <- 1964
    city.rates <- subset(city.rates, !(placename=='London' & year>1963))
    city.rates <- rbind(city.rates, tmp.rates)
    rm(tmp.rates)
    where <- 'England & Wales'
})

## things that are the same for each country
## reprates, limits, 
mod.prep <- llply(mod.prep, function(.list) within(.list,{
    ## max and proportion of zeros for each city
    ## was ranges
    city.stats <- mk.city.stats(list(cases.obs))
    index.limits = paste(range(index(cases.obs)),collapse='::')
    ## grab rates within case report range
    .tmp.rates <- subset(city.rates, year %in% (.indexyear(cases.obs)+1900))
    ## get median pop from available years
    ## median pop?? mean pop??
    demog <- ddply(.tmp.rates, 'placename', function(x) data.frame(pop=median(x$pop)))
    ## city connection matrix, not currently used
    demog <- mk.connect(city.ranks, demog)
    ## make absolutely certain they're ordered by pop
    demog <- demog[ order(-demog$pop),]
    demog$placename <- with(demog, 
            ordered(placename, levels=unique(placename)))
    city.order <- levels(demog$placename)
    ## exclude last year, since demogs are calcd from diffs
    city.rates <- na.omit(city.rates)
    ## ensure order
    city.rates$placename <- with(city.rates, ordered(placename, levels=city.order))
    city.rates <- with(city.rates, city.rates[ order(placename, year),])
    ##
    ## make reprate first, before adding dS and dR to city.rates
    reprates <- mk.reprates(myR0s, city.rates, cases.obs) 
    demog <- merge( demog, reprates)
    ## then make dS and dR by R0
    #
    city.rates <- ldply(myR0s, mk.cityrates, city.rates=city.rates)
}))


## fix!!
if(F) {
## analysis - infer distributions for data
## deprecated?
mod.prep <- llply(mod.prep, function(.list) within(.list,{
    ## get the dist for each R0
    distrib.df <- ldply(myR0s, function(r0) {
                mk.all.distrib(cases.obs, subset( demog, R0==r0))
    })
    distrib.df$where <- where
    distrib.df$model <- factor('Data')
    ## below is deprecated?
    ## make spectrum, biweekly   
    #spec <- mk.spec.mat(gren.cases[1:512,], freq.adjust=26)
    #spec.names <-  colnames(spec)
}))
}
