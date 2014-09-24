## from run.prepmod.R
{
    ## pull out demographics 
    alldemog <- ldply(mod.prep, function(x) x$newdemog)
    alldemog <- rename(alldemog, c(.id='where'))
    ## pull out reporting rates for both cities into single df
    all.city.rates <- ldply(mod.prep, function(x) x$city.rates)
    all.city.rates <- rename(all.city.rates, c(.id='where'))
}
