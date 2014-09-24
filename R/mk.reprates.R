#' Helper function 
#' get dS and dR from birth, infant, migration, death, and r0
#' return the modified city.rates
mk.cityrates <- function(r0, city.rates) {
    ret <- within(city.rates, {
        ## yearly per cap change in susceptibles
        dS <- (birth*(1-infant) + migrate/r0)
        ## yearly per cap change in recovereds
        dR <- (migrate*(1-1/r0) - death)
    })
    ret$R0 <- r0
    return(ret)
}

#' get reporting rate from case reports,
#' yearly demographic rates, and r0
mk.reprates <- function(r0vec, city.rates, cases.obs) {
    ## actual cases
    cases.a <- adply( cases.obs, 2, sum, na.rm=T)
    colnames(cases.a) <- c('placename', 'Ctotal')
    ## expected cases, prep
    ## make sure we only get births for timeseries years
    city.rates <- subset(city.rates, year %in% (.indexyear(cases.obs)+1900))
    ## for each R0, should deprecate, R0 doesn't do much!!
    cases.e <- ldply(r0vec, function(r0){
        ret <- ddply(city.rates, 'placename', function(x) {
            ## total flow of S (live births + migrating s) x prop of S infected (1-1/r0)
            Syearly = with(x,  (birth*(1-infant) + migrate/r0)*pop)*(1-1/r0)
            ## these should be counts, round to integer
            Stotal = round(sum(Syearly))
            ## 
            data.frame(Stotal)
        })
        ret$R0 <- r0
        ret
    })
    ## combine actual and expected, get reprate
    cases <- merge(cases.a, cases.e)
    cases$reprate <- with(cases, Ctotal/Stotal)
    ## warn if 
    .bad <- subset(cases, reprate>1)
    if(nrow(.bad)>1 ) {
        warning('In mk.reprate: Reporting rate >1.\n')
        print(.bad)
    }
    return( subset(cases, select=c(placename, R0, reprate, Ctotal, Stotal)))
}
