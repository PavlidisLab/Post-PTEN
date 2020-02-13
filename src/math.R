# description: "Mathematical functions or shortcuts to calculate things quickly."
# author: "Manuel Belmadani"
# date: "23/04/2018"
# updated: "23/04/2018"


stderr <- 
  # Standard error of the mean.
  function(x) sqrt( var(x,na.rm=TRUE)/length(na.omit(x)) )

confidence.intervals <-
  # Calculate the confidence interval error
  function( sample.mean,
            sample.n,
            sample.sd=NULL,
            sample.sem=NULL ) {
    if (is.null(sample.sd) &
        is.null(sample.sem)){
        stop("Error: Need to specify either SD or SEM.")
    } else if ( is.null(sample.sd) ) {
      sample.sd = sample.sem * sqrt(sample.n)  
    }
    
    error <- qt(0.975,df=sample.n-1)*sample.sd/sqrt(sample.n)
    return(error)
}
