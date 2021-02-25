#' Sample size
#'
#' Returns necessary sample size for a given power and dataset
#'
#' @param alpha Numeric alpha.
#' @param margin Numeric error margin.
#' @param dependent.variable The variable being sampled as a vector.
#' @param weights Numeric weights (if present) as a vector.
#' @return A list containing the parameters and a data.frame with the smoothed estimates.
#' @keywords internal
sample.size<-function(alpha,margin,dependent.variable,weights){
	if(missing(weights)){
		if(is.numeric(dependent.variable)){stdev<-stats::sd(dependent.variable);}
		else{
			stdev<-stats::sd(dependent.variable==names(which.max(table(factor(dependent.variable)))));
		}
	}else{
		if(sum(weights)<1){
			return(NA);
		}else{
			if(is.numeric(dependent.variable)){
				stdev<-sqrt(Hmisc::wtd.var(x=dependent.variable,weights=weights));
			}else{
				wtbl<-Hmisc::wtd.table(x=dependent.variable,weights=weights);
				stdev<-sqrt(Hmisc::wtd.var(x=(dependent.variable==wtbl$x[which.max(wtbl$sum.of.weights)]),weights=weights));
			}
		}
	}
	z=stats::qnorm(alpha/2);
	return((stdev*(z/margin))^2);
}
#' Sample size floored variance
#'
#' Returns necessary sample size for a given power and dataset, but treats any standard deviation below 0.1 as 0.1
#'
#' @param alpha Numeric alpha.
#' @param margin Numeric error margin.
#' @param dependent.variable The variable being sampled as a vector.
#' @param weights Numeric weights (if present) as a vector.
#' @return A list containing the parameters and a data.frame with the smoothed estimates.
#' @keywords internal
sample.size.floored.variance<-function(alpha,margin,dependent.variable,weights){
	if(missing(weights)){
		if(is.numeric(dependent.variable)){stdev<-stats::sd(dependent.variable);}
		else{
			stdev<-stats::sd(dependent.variable==names(which.max(table(factor(dependent.variable)))));
		}
	}else{
		if(is.numeric(dependent.variable)){
			stdev<-sqrt(Hmisc::wtd.var(x=dependent.variable,weights=weights));
		}else{
			wtbl<-Hmisc::wtd.table(x=dependent.variable,weights=weights);
			stdev<-sqrt(Hmisc::wtd.var(x=(dependent.variable==wtbl$x[which.max(wtbl$sum.of.weights)]),weights=weights));
		}
	}
	if(is.nan(stdev) | stdev<0.1){stdev=0.1;}
	z=stats::qnorm(alpha/2);
	return((stdev*(z/margin))^2);
}
