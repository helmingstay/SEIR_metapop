#' make  matrix of connection between cities
#' for use is simulation
#' results not currently used...
mk.connect <- function(ranks, demog, do.norm=T)  {
    ## compute new distance metric
    ## kinda like a ranked gravity model
    ## take log to get a reasonable scale
    ## normalize by max if do.norm=T
    ## merge into demog and return 
    city.connect = ddply(ranks, 'placename', function(x) {
        ## the main distance function
        ## meaning of ^2 for proportion?  Keep proportion-units?
        ret <- data.frame(connect=log10(sum((x$poprank1*x$poprank)/(x$distrank)^2)))
        return(ret)
        ## other things we've tried
        #data.frame(connect=sum(log10(x$pop)+log10(x$pop1)- log10(x$dist)))
        #data.frame(connect=log10(sum(1/(x$poprank1*x$rank*x$poprank))))
        #data.frame(connect=-log10(sum(1/(x$poprank1*x$rank))/unique(x$poprank)))
        #data.frame(connect=log10(sum(1/(x$poprank1*x$rank))))
    })
    ## normalize?
    if ( do.norm )  city.connect$connect <- with(city.connect, connect/max(connect))
    ## negative imports makes bad things happen
    if ( min( city.connect$connect) <= 0 ) stop("Connection must be positive")
    ## order by the final ranking
    ## only keep placename and connect
    city.connect <- 
            city.connect[order(city.connect$connect), c('placename', 'connect')]
    ## finally, merge into demog and return
    demog <- merge( city.connect, demog )
    return(demog)
}
