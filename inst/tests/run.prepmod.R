## prep for most of the run.*.R files:
## run.comine, run.newfpc, run.ensemble...
## prep objects for model run
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

mod.prep = list()
#myR0s = c(18,20)
#myR0s = 20
myR0s = 10:25
my.divide=1
my.which=1

mod.prep$US <- within(list(), {
    cases.obs <- lcasedata$Measles$city.zoo["::1945-12",]
    #!!
    cases.obs <- xts.split(cases.obs, .divide=my.divide, .which=my.which)
    ## weeks per sample
    period <- 1
    ## max and proportion of zeros for each city
    ranges <- mk.range(list(cases.obs))
    index.limits = paste(range(index(lcasedata$Measles$city.zoo)),collapse='::')
    ## for mk.alldist
    ## for US, only get demog cols that are also in gren
    #newdemog = subset(ldemog$city.dec, year == ldemog$order.year,
                               #select=c(placename, pop, year ))
    ##
    city.rates <- ldemog$city.rates
    .tmp.rates <- subset(city.rates, year %in% (.indexyear(cases.obs)+1900))
    newdemog <- ddply(.tmp.rates, 'placename', function(x) data.frame(pop=median(x$pop)))
    ##
    newdemog <- mk.connect(lcasedata$Measles$city.ranks, newdemog)

    ## make absolutely certain they're ordered by pop
    newdemog <- newdemog[ order(-newdemog$pop),]
    newdemog$placename <- with(newdemog, 
            ordered(placename, levels=unique(placename)))
    city.order <- levels(newdemog$placename)
    ## ensure order
    city.rates$placename <- with(city.rates, ordered(placename, levels=city.order))
    city.rates <- with(city.rates, city.rates[ order(placename, year),])
    ## make reprate first, before adding dS and dR to city.rates
    reprates <- mk.reprates(myR0s, city.rates, cases.obs) 
    newdemog <- merge( newdemog, reprates)
    ## then make dS and dR by R0
    city.rates <- ldply(myR0s, mk.cityrates, city.rates=city.rates)
    ## 
    ## get the dist for each R0
    alldist.df <- ldply(myR0s, function(r0) {
                mk.alldist(cases.obs, subset( newdemog, R0==r0))
    })
    alldist.df$where='US'
    alldist.df$model <- factor('Data')
    
    #spec <- mk.spec.mat(lcasedata$Measles$city.approx[1:1024])
    #spec.names <-  colnames(spec)
})

## same for England and Wales
mod.prep[['England & Wales']] <- within(lgrenfell, {
    ## to match names in US
    ## there are several administrative changes in 1966, which affects 1965 migration
    ## only use TS (and run model) through end of 1964
    cases.obs <- gren.cases["::1964",]
    #!!
    cases.obs <- xts.split(cases.obs, .divide=my.divide, .which=my.which)
    ## weeks per sample
    period <- 2

    ## London changes in 1966.  Reuse values from 1063 in 1964
    tmp.rates <- subset(city.rates, placename=='London' & year==1963)
    tmp.rates$year <- 1964
    city.rates <- subset(city.rates, !(placename=='London' & year>1963))
    city.rates <- rbind(city.rates, tmp.rates)
    rm(tmp.rates)
    ## get median pop from available years
    .tmp.rates <- subset(city.rates, year %in% (.indexyear(cases.obs)+1900))
    newdemog <- ddply(.tmp.rates, 'placename', function(x) data.frame(pop=median(x$pop)))
    ## max and proportion of zeros for each city
    ranges <- mk.range(list(cases.obs))
    index.limits = paste(range(index(cases.obs)),collapse='::')
    ## city connectivity
    newdemog <- mk.connect(lgrenfell$city.ranks, newdemog)
    ## make absolutely sure of order
    newdemog <- newdemog[ order(-newdemog$pop),]
    newdemog$placename <- with(newdemog, 
            ordered(placename, levels=unique(placename)))
    city.order <- levels(newdemog$placename)
    ## exclude last year, since demogs are calcd from diffs
    city.rates <- na.omit(city.rates)
    ## ensure order
    city.rates$placename <- with(city.rates, ordered(placename, levels=city.order))
    city.rates <- with(city.rates, city.rates[ order(placename, year),])
    ##
    ## make reprate first, before adding dS and dR to city.rates
    reprates <- mk.reprates(myR0s, city.rates, cases.obs) 
    newdemog <- merge( newdemog, reprates)
    ## then make dS and dR by R0
    ##
    city.rates <- ldply(myR0s, mk.cityrates, city.rates=city.rates)
    ## get the dist for each R0
    alldist.df <- ldply(myR0s, function(r0) {
                mk.alldist(cases.obs, subset( newdemog, R0==r0))
    })
    alldist.df$where='England & Wales'
    alldist.df$model <- factor('Data')
    ## make spectrum, biweekly   
    #spec <- mk.spec.mat(gren.cases[1:512,], freq.adjust=26)
    #spec.names <-  colnames(spec)
})

alldemog <- ldply(mod.prep, function(x) x$newdemog)
alldemog <- colname.replace('.id','where', alldemog)
alldemog <- wherefix(alldemog)

all.city.rates <- ldply(mod.prep, function(x) x$city.rates)
all.city.rates <- colname.replace('.id','where', all.city.rates)
all.city.rates  <- wherefix(all.city.rates)

