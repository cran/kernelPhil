#' Kernel smooth data in space alone
#'
#' This function performs kernel smoothing on a dataset in space alone.
#'
#' @param dataset Dataset to be smoothed as a data.frame.
#' @param dependent.variable String name of the single column in dataset with the factor dependent variable (if data.type=="factor") or a vector of column names with numeric counts (if data.type=="count") (defaults to "dependent.variable").
#' @param x String name of column containing numeric x co-ordinate (defaults to "x").
#' @param y String name of column containing numeric y co-ordinate (defaults to "y").
#' @param weight String name of column in the dataset with numeric weights (defaults to "weight").
#' @param normalise.by String name of column by which data should be normalised (typically factor with document, speaker or writer ids).
#' @param data.type The type of the dependent variable: either "factor", if each row is a token, or "count", if each row is a document, speaker or writer with token counts in separate columns (defaults to "factor").
#' @param alpha Numeric alpha for calculating error margins (defaults to 0.05).
#' @param margin Numeric desired error margin for calculating spatial bandwidths (defaults to 0.1).
#' @param kernel.function The kernel function, one of gaussian.kernel, gaussian.square.kernel, triangular.kernel, square.kernel, or a custom function (defaults to gaussian.kernel).
#' @param adaptive.spatial.bw A boolean indicating whether the spatial bandwidth is adaptive (set to achieve margin at every point) or static (set to the average of bandwidths needed to achieve margin at every point).
#' @param measure.points A data.frame of spatial points at which estimates are to be made, with two columns with the same names as x,y in dataset; if not supplied, estimates are at the same locations as dataset.
#' @param projection The spatial projection as a proj4 string - if given, data will be projected before smoothing and results will be deprojected before returning.
#' @param round.up.low.variance Set to TRUE if there are periods of time with extremely low variance (defaults to TRUE).
#' @param explicit If TRUE, progress will be reported with a progress bar (defaults to TRUE).
#' @return A data.frame with the smoothed estimates.
#' @examples
#' n=400;
#' synthesised.data<-data.frame(x=stats::runif(n),y=stats::runif(n),
#'     year=stats::runif(n,0,sqrt(2)));
#' synthesised.data$dependent.variable<-unlist(lapply(1:nrow(synthesised.data),
#'     function(X){
#'     stats::dist(as.matrix(synthesised.data[c(1,X),1:2]),method =
#'         "euclidean")<synthesised.data$year[X];
#' }))
#' result<-kernelPhil::kernel.smooth.in.space(dataset = synthesised.data);
#' ggplot2::ggplot(result,ggplot2::aes(x=x,y=y,colour=relative_density_TRUE))+
#'     ggplot2::geom_point();
#' @importFrom dplyr %>%
#' @export
kernel.smooth.in.space<-function(dataset,dependent.variable="dependent.variable",x="x",y="y",weight="weight",normalise.by,data.type="factor",alpha=0.05,margin=0.1,kernel.function=gaussian.kernel,adaptive.spatial.bw=TRUE,measure.points,projection=NA,round.up.low.variance=TRUE,explicit=TRUE){
	weight.consistent.name<-dataset.wide<-NULL; # set up variable names for use later
	if(round.up.low.variance){sample.size.func<-sample.size.floored.variance;}else{sample.size.func<-sample.size;}

	# make sure variables present in the right formats
	if(data.type=="factor"){
		dataset[,dependent.variable]<-factor(dataset[,dependent.variable]);
	}else if(data.type=="count"){
		if(!is.vector(dependent.variable)){stop("if data.type is 'count', dependent.variable should be a vector of column names");}
		for(dependent.variable.i in dependent.variable){
			if(class(dataset[,dependent.variable.i])=="factor"){
				dataset[,dependent.variable.i]<-as.numeric(as.character(dataset[,dependent.variable.i]));
			}else if(class(dataset[,dependent.variable.i])=="character"){
				dataset[,dependent.variable.i]<-as.numeric(dataset[,dependent.variable.i]);
			}else if(class(dataset[,dependent.variable.i])!="numeric"){
				stop("for data.type='count', data should be numeric");
			}
		}
	}else{
		stop("data.type should be one of 'factor' or 'count'");
	}
	if(!(weight %in% colnames(dataset))){
		dataset[,weight]<-1;
	}

	# get rid of items with na location
	dataset<-dataset[which(!is.na(dataset[,x]) & !is.na(dataset[,y])),];

	# project dataset if needed
	if(!is.na(projection)){dataset[,c(x,y)]<-rgdal::project(as.matrix(dataset[,c(x,y)]),proj=projection);}

	# if a column name is passed with "normalise.by", then divide weight by count of factor levels (used for normalising weights by text id, speaker id, etc.)
	if(!missing(normalise.by)){
		dataset[,weight]<-unlist(apply(dataset,1,function(X) return(as.numeric(X[weight])/nrow(dataset[which(dataset[,normalise.by]==X[normalise.by]),]))));
	}

	# set up wide dataset with weights per possible value at each location
	if(data.type=="count"){  # if counts in separate columns, first put this into long form
		uncount<-function(df,dep.levs){
			df$temp.id.var<-1:nrow(df);
			non.dep.cols<-colnames(df)[which(!(colnames(df) %in% dep.levs))];
			uncount1<-function(df,dep,non.dep.cols){
					tdf<-df[rep(1:nrow(df),df[,dep]),c(non.dep.cols)];
					tdf$dependent.variable<-dep;
					tdf
			}
			df<-do.call("rbind",lapply(dep.levs,uncount1,df=df,non.dep.cols=non.dep.cols));
			df[,weight]<-df[,weight]/unlist(lapply(df$temp.id.var,function(X) length(df$temp.id.var[which(df$temp.id.var==X)])));
			df$temp.id.var<-NULL;
			df
		}
		dataset<-uncount(dataset,dependent.variable);
		data.type="factor";
		dependent.variable="dependent.variable"
		dataset[,dependent.variable]<-factor(dataset[,dependent.variable]);
	}
	if(data.type=="factor"){
		dataset$weight.consistent.name<-dataset[,weight];
		# reshape data so that we have one row per coordinate pair with counts of each level of the factor multiplied by weights
		suppressMessages(dataset.wide<-dplyr::group_by_(dataset,x,y,dependent.variable) %>% dplyr::summarise(count=dplyr::n(), sum_weight=sum(weight.consistent.name,na.rm=TRUE)) %>% as.data.frame %>% reshape2::melt(id.var=c(x,y,dependent.variable)) %>% reshape2::dcast(sprintf("%s+%s~variable+%s",x,y,dependent.variable)));
		dataset$weight.consistent.name<-NULL;

		# list levels of the factor
		dependent.variable.levels<-levels(dataset[,dependent.variable]);

		# add column with the total weight per point
		dataset.wide$sum_weight<-apply(dataset.wide,1,function(X) sum(X[paste0("sum_weight_",dependent.variable.levels)],na.rm=TRUE));
	}

	# if no separate measure points passed, use the coordinates of the dataset as measure points
	if(missing(measure.points)){
		measure.points=dataset.wide[,c(x,y)];
	}else{
		# project measure points if needed
		if(!is.na(projection)){measure.points[,c(x,y)]<-rgdal::project(as.matrix(measure.points[,c(x,y)]),proj=projection);}
	}

	# calculate distance matrix
	dm<-wordspace::dist.matrix(as.matrix(dataset.wide[,c(x,y)]),as.matrix(measure.points),method="euclidean");
	colnames(dm)<-1:ncol(dm);
	rownames(dm)<-1:nrow(dm);

	# identify minimum sample size needed to achieve desired specificity and confidence (k)
	k=sample.size.func(alpha=alpha,margin=margin,dependent.variable=dataset[,dependent.variable],weights=dataset[,weight]);
	if(k>sum(dataset.wide$sum_weight)){stop("not enough data to achieve this level of specificity and confidence");}
	if(k==0){k=min(dataset.wide$sum_weight[which(dataset.wide$sum_weight>0)]);}


	# identify local bandwidths, depending on kernel
	find.bandwidth<-function(i){
		stats::optimise(function(spatial.bandwidth){
			abs(k-sum(kernel.function(dm[,i],as.numeric(spatial.bandwidth))*dataset.wide$sum_weight));
		},lower=min(dm[which(dm>0)]),upper=max(dm),tol=0.1,maximum = FALSE)$minimum;
	};
	if(explicit){
		message("identifying local bandwidths\n");
		measure.points$bw<-unlist(pbapply::pblapply(1:nrow(dataset.wide),FUN = find.bandwidth));
	}else{
		measure.points$bw<-unlist(lapply(1:nrow(dataset.wide),FUN = find.bandwidth));
	}

	# if static bandwidth, take the mean bandwidth
	if(!adaptive.spatial.bw){
		measure.points$bw<-mean(dataset.wide$bw);
	}

	apply.kernel.smooth.by.location<-function(crow){
		# apply spatial kernel
		dataset.wide$current.kernel.value<<-unlist(lapply(dm[,crow],kernel.function,bandwidth=measure.points$bw[crow]));

		# sum weights per variant across the relevant locations
		cresult<-stats::setNames(data.frame(rbind(unlist(lapply(dependent.variable.levels,function(dependent.variable.level) return(sum(dataset.wide[,paste0("sum_weight_",dependent.variable.level)]*dataset.wide$current.kernel.value,na.rm=TRUE)))))),paste0("sum_weight_",dependent.variable.levels));

		cresult$sum_of_weights_total<-sum(cresult);
		cresult[,paste0("relative_density_",dependent.variable.levels)]<-cresult[,1:ncol(cresult)-1]/cresult$sum_of_weights_total;

		csds<-lapply(dependent.variable.levels,function(dependent.variable.level){return();})

		cresult$effective_sample_size<-sum(dataset.wide$current.kernel.value*dataset.wide$sum_weight);
		cresult$sum_of_weights_total<-NULL;
		return(cresult);
	}

	if(explicit){
		message("calculating kernel smooths\n");
		results<-as.data.frame(do.call(rbind,pbapply::pblapply(1:nrow(measure.points),apply.kernel.smooth.by.location)));
	}else{
		results<-as.data.frame(do.call(rbind,lapply(1:nrow(measure.points),apply.kernel.smooth.by.location)));
	}

	# reassociate with locations
	results<-cbind(measure.points,results[,c(paste0("relative_density_",dependent.variable.levels),"effective_sample_size")]);
	results<-merge(results,dataset.wide[,c(x,y,"sum_weight")],by=c(x,y));
	colnames(results)[length(colnames(results))]<-"weight_at_point";

	results$best<-unlist(apply(results,1,function(X){
	  substr(labels(X)[which(substr(labels(X),0,17)=="relative_density_")][which.max(X[labels(X)[which(substr(labels(X),0,17)=="relative_density_")]])],18,1000)
	}));

	# deproject dataset if needed
	if(!is.na(projection)){results[,c(x,y)]<-rgdal::project(as.matrix(results[,c(x,y)]),proj=projection,inv=TRUE);}

	return(results);
}
#' Kernel smooth data in space alone, returning specific error margins at each point
#'
#' This function performs kernel smoothing on a dataset in space alone. It is the same as kernel.smooth.in.space(), except that the results include the error margins for the estimates at every point. Note that it is *much* slower than kernel.smooth.in.space().
#'
#' @param dataset The dataset to be smoothed as a data.frame.
#' @param dependent.variable String name of the single column in dataset with the factor dependent variable (if data.type=="factor") or a vector of column names with numeric counts (if data.type=="count") (defaults to "dependent.variable").
#' @param x String name of column containing numeric x co-ordinate (defaults to "x").
#' @param y String name of column containing numeric y co-ordinate (defaults to "y").
#' @param weight String name of column in the dataset with numeric weights (defaults to "weight").
#' @param normalise.by String name of column by which data should be normalised (typically factor with document, speaker or writer ids).
#' @param data.type The type of the dependent variable as a string: either "factor", if each row is a token, or "count", if each row is a document, speaker or writer with token counts in separate columns (defaults to "factor").
#' @param alpha Numeric alpha for calculating error margins (defaults to 0.05).
#' @param margin Numeric desired error margin for calculating spatial bandwidths (defaults to 0.1).
#' @param kernel.function The kernel function, one of gaussian.kernel, gaussian.square.kernel, triangular.kernel, square.kernel, or a custom function (defaults to gaussian.kernel).
#' @param adaptive.spatial.bw A boolean indicating whether the spatial bandwidth is adaptive (set to achieve margin at every point) or static (set to the average of bandwidths needed to achieve margin at every point).
#' @param measure.points A data.frame of spatial points at which estimates are to be made, with two columns with the same names as x,y in dataset; if not supplied, estimates are at the same locations as dataset.
#' @param round.up.low.variance Set to TRUE if there are periods of time with extremely low variance (defaults to TRUE).
#' @param projection The spatial projection as a proj4 string - if given, data will be projected before smoothing and results will be deprojected before returning.
#' @param explicit If TRUE, progress will be reported with a progress bar (defaults to TRUE).
#' @return A data.frame with the smoothed estimates.
#' @examples
#' \donttest{n=400;
#' synthesised.data<-data.frame(x=stats::runif(n),y=stats::runif(n),
#'     year=stats::runif(n,0,sqrt(2)));
#' synthesised.data$dependent.variable<-unlist(lapply(1:nrow(synthesised.data),
#'     function(X){
#'     stats::dist(as.matrix(synthesised.data[c(1,X),1:2]),method =
#'         "euclidean")<synthesised.data$year[X];
#' }))
#' result<-kernelPhil::kernel.smooth.in.space.with.margins(dataset = synthesised.data);
#' ggplot2::ggplot(result,ggplot2::aes(x=x,y=y,colour=relative_density_TRUE))+
#'     ggplot2::geom_point();}
#' @export
kernel.smooth.in.space.with.margins<-function(dataset,dependent.variable="dependent.variable",x="x",y="y",weight="weight",normalise.by,data.type="factor",alpha=0.05,margin=0.1,kernel.function=gaussian.kernel,adaptive.spatial.bw=TRUE,measure.points,projection=NA,round.up.low.variance=TRUE,explicit=TRUE){
	dataset.wide<-NULL; # set up variable names for use later
	if(round.up.low.variance){sample.size.func<-sample.size.floored.variance;}else{sample.size.func<-sample.size;}

	# make sure variables present in the right formats
	if(data.type=="factor"){
		dataset[,dependent.variable]<-factor(dataset[,dependent.variable]);
	}else if(data.type=="count"){
		if(!is.vector(dependent.variable)){stop("if data.type is 'count', dependent.variable should be a vector of column names");}
		for(dependent.variable.i in dependent.variable){
			if(class(dataset[,dependent.variable.i])=="factor"){
				dataset[,dependent.variable.i]<-as.numeric(as.character(dataset[,dependent.variable.i]));
			}else if(class(dataset[,dependent.variable.i])=="character"){
				dataset[,dependent.variable.i]<-as.numeric(dataset[,dependent.variable.i]);
			}else if(class(dataset[,dependent.variable.i])!="numeric"){
				stop("for data.type='count', data should be numeric");
			}
		}
	}
	if(!(weight %in% colnames(dataset))){
		dataset$weight<-1;
	}else{
		dataset$weight<-dataset[,weight];
	}

	# get rid of items with na location
	dataset<-dataset[which(!is.na(dataset[,x]) & !is.na(dataset[,y])),];

	# project dataset if needed
	if(!is.na(projection)){dataset[,c(x,y)]<-rgdal::project(as.matrix(dataset[,c(x,y)]),proj=projection);}

	# if a column name is passed with "normalise.by", then divide weight by count of factor levels (used for normalising weights by text id, speaker id, etc.)
	if(!missing(normalise.by)){
		dataset$weight<-unlist(apply(dataset,1,function(X) return(as.numeric(X["weight"])/nrow(dataset[which(dataset[,normalise.by]==X[normalise.by]),]))));
	}

	# set up wide dataset with weights per possible value at each location
	if(data.type=="count"){  # if counts in separate columns, first put this into long form
		uncount<-function(df,dep.levs){
			df$temp.id.var<-1:nrow(df);
		  non.dep.cols<-colnames(df)[which(!(colnames(df) %in% dep.levs))];
		  uncount1<-function(df,dep,non.dep.cols){
		      tdf<-df[rep(1:nrow(df),df[,dep]),c(non.dep.cols)];
		      tdf$dependent.variable<-dep;
		      tdf
		  }
		  df<-do.call("rbind",lapply(dep.levs,uncount1,df=df,non.dep.cols=non.dep.cols));
			df[,weight]<-df[,weight]/unlist(lapply(df$temp.id.var,function(X) length(df$temp.id.var[which(df$temp.id.var==X)])));
			df$temp.id.var<-NULL;
			df
		}
		dataset<-uncount(dataset,dependent.variable);
		data.type="factor";
		dependent.variable="dependent.variable"
		dataset[,dependent.variable]<-factor(dataset[,dependent.variable]);
	}
	if(data.type=="factor"){ # long>wide for factors
		dependent.variable.levels<-levels(dataset[,dependent.variable]);
		dataset.prepped<-dataset[,c(x,y,weight)];
		dataset.prepped[,paste0("count_",dependent.variable.levels)]<-matrix(as.numeric(unlist(lapply(dependent.variable.levels,function(dependent.variable.level){dataset[,dependent.variable]==dependent.variable.level;})),ncol=length(dependent.variable.levels)));
		dataset.prepped[,paste0("sum_weight_",dependent.variable.levels)]<-matrix(as.numeric(unlist(lapply(dependent.variable.levels,function(dependent.variable.level){dataset.prepped[,paste0("count_",dependent.variable.level)]*dataset.prepped[,weight];}))));
		colnames(dataset.prepped)[3]<-"sum_weight";
	}

	# replace NAs with 0s
	for(dependent.variable.level in dependent.variable.levels){
		dataset.prepped[which(is.na(dataset.prepped[,paste0("sum_weight_",dependent.variable.level)])),paste0("sum_weight_",dependent.variable.level)]<-0;
		dataset.prepped[which(is.na(dataset.prepped[,paste0("count_",dependent.variable.level)])),paste0("count_",dependent.variable.level)]<-0;
	}

	# create df of locations in dataset
	location.index.df<-unique(dataset.prepped[,c(x,y)]);

	# if no separate measure points passed, use the coordinates of the dataset as measure points
	if(missing(measure.points)){
		measure.points=location.index.df;
	}else{
		# project measure points if needed
		if(!is.na(projection)){measure.points[,c(x,y)]<-rgdal::project(as.matrix(measure.points[,c(x,y)]),proj=projection);}
	}

	# check available system memory
	if(is.numeric(benchmarkme::get_ram())){
		avmem<-as.numeric(benchmarkme::get_ram())*1073741824;
	}else{
		avmem<-tryCatch({strsplit(system("wmic MemoryChip get Capacity", intern=TRUE), " ")[[2]][1]},error=function(cond){
			warning("couldn't determine system memory; assuming enough memory is available\n");
			return(nrow(dataset.prepped)*nrow(measure.points)*8);
		})
	}
	if(avmem<nrow(dataset.prepped)*nrow(measure.points)*8){
		if(avmem<nrow(location.index.df)*nrow(measure.points)*8){
			stop("not enough memory available even for low-memory version of script");
		}else{
			warning("low available system memory; switching to slower, lower-memory version\n");
			memflag=TRUE;
		}
	}else{
		memflag=FALSE;
	}

	# calculate distance matrix

	dm<-wordspace::dist.matrix(as.matrix(location.index.df),as.matrix(measure.points),method="euclidean");
	colnames(dm)<-1:ncol(dm);
	rownames(dm)<-1:nrow(dm);
	location.index.df$index<-1:nrow(location.index.df);
	dataset.prepped$location_index<- unlist(apply(dataset.prepped,1,function(X){return(location.index.df$index[which(location.index.df[,x] == as.numeric(X[x]) & location.index.df[,y]==as.numeric(X[y]))]);}));

	if(!memflag){
		dm<-dm[dataset.prepped$location_index,];
	}

	# identify minimum sample size needed to achieve desired specificity and confidence (k)
	k=sample.size.func(alpha=alpha,margin=margin,dependent.variable=dataset[,dependent.variable],weights=dataset[,weight]);
	if(k>sum(dataset.prepped$sum_weight)){stop("not enough data to achieve this level of specificity and confidence");}
	if(k==0){k=min(dataset.wide$sum_weight[which(dataset.wide$sum_weight>0)]);}


	# identify local bandwidths, depending on kernel
	if(memflag){ # memory-safe, slow version
		find.bandwidth<-function(i){
			stats::optimise(function(spatial.bandwidth){
				abs(k-sum(kernel.function(dm[dataset.prepped$location_index,i],as.numeric(spatial.bandwidth))*dataset.prepped$sum_weight));
			},lower=min(dm[which(dm>0)]),upper=max(dm),tol=0.1,maximum = FALSE)$minimum;
		};
	}else{ # memory-hungry normal version
		find.bandwidth<-function(i){
			stats::optimise(function(spatial.bandwidth){
				abs(k-sum(kernel.function(dm[,i],as.numeric(spatial.bandwidth))*dataset.prepped$sum_weight));
			},lower=min(dm[which(dm>0)]),upper=max(dm),tol=0.1,maximum = FALSE)$minimum;
		};
	}
	if(explicit){
		message("identifying local bandwidths\n");
		measure.points$bw<-unlist(pbapply::pblapply(1:nrow(measure.points),FUN = find.bandwidth));
	}else{
		measure.points$bw<-unlist(lapply(1:nrow(measure.points),FUN = find.bandwidth));
	}

	# if static bandwidth, take the mean bandwidth
	if(!adaptive.spatial.bw){
		dataset.prepped$bw<-mean(dataset.prepped$bw);
	}

	dataset.prepped[,paste0("prop_",dependent.variable.levels)]<-dataset.prepped[,paste0("sum_weight_",dependent.variable.levels)]/dataset.prepped$sum_weight;

	if(memflag){ # memory-safe, slow version
		apply.kernel.smooth.by.location<-function(crow){
			# apply spatial kernel
			current_kernel<-unlist(lapply(dm[dataset.prepped$location_index,crow],kernel.function,bandwidth=measure.points$bw[crow]));
			current_weights<-current_kernel*dataset.prepped$sum_weight;
			current_n<-sum(current_weights);

			# sum weights, proportions and standard deviations per variant across the relevant locations
			cresult<-as.list(unlist(lapply(dependent.variable.levels,function(dependent.variable.level){
				current_mean<-sum(dataset.prepped[,paste0("prop_",dependent.variable.level)]*current_weights,na.rm=TRUE)/current_n;
				current_var<-Hmisc::wtd.mean(x = (dataset.prepped[,paste0("prop_",dependent.variable.level)]-current_mean)^2,current_weights,normwt=TRUE);
				current_sd<-sqrt(current_var);
				if(current_sd>1){browser();}
				current_confint<- -1*stats::qnorm(alpha/2)*(current_sd/sqrt(current_n));
				return(stats::setNames(c(current_mean,current_sd,current_confint),c(paste0("relative_density_",dependent.variable.level),paste0("sd_",dependent.variable.level),paste0("confint_",dependent.variable.level))));
			})));

			cresult$effective_sample_size<-current_n;

			return(unlist(cresult));
		}
	}else{
		apply.kernel.smooth.by.location<-function(crow){
			# apply spatial kernel
			current_kernel<-unlist(lapply(dm[,crow],kernel.function,bandwidth=measure.points$bw[crow]));
			current_weights<-current_kernel*dataset.prepped$sum_weight;
			current_n<-sum(current_weights);

			# sum weights, proportions and standard deviations per variant across the relevant locations
			cresult<-as.list(unlist(lapply(dependent.variable.levels,function(dependent.variable.level){
				current_mean<-sum(dataset.prepped[,paste0("prop_",dependent.variable.level)]*current_weights,na.rm=TRUE)/current_n;
				current_var<-Hmisc::wtd.mean(x = (dataset.prepped[,paste0("prop_",dependent.variable.level)]-current_mean)^2,current_weights,normwt=TRUE);
				current_sd<-sqrt(current_var);
				if(current_sd>1){browser();}
				current_confint<- -1*stats::qnorm(alpha/2)*(current_sd/sqrt(current_n));
				return(stats::setNames(c(current_mean,current_sd,current_confint),c(paste0("relative_density_",dependent.variable.level),paste0("sd_",dependent.variable.level),paste0("confint_",dependent.variable.level))));
			})));

			cresult$effective_sample_size<-current_n;

			return(unlist(cresult));
		}
	}

	if(explicit){
		message("calculating kernel smooths\n");
		results<-as.data.frame(do.call(rbind,pbapply::pblapply(1:nrow(measure.points),apply.kernel.smooth.by.location)));
	}else{
		results<-as.data.frame(do.call(rbind,lapply(1:nrow(measure.points),apply.kernel.smooth.by.location)));
	}

	# reassociate with locations
	results<-cbind(measure.points,results);
	results<-merge(results,stats::aggregate(stats::as.formula(paste0("sum_weight~",x,"+",y)),dataset.prepped[,c(x,y,"sum_weight")],FUN=sum),by=c(x,y),all.x=TRUE,all.y=FALSE);
	colnames(results)[length(colnames(results))]<-"weight_at_point";

	results$best<-unlist(apply(results,1,function(X){
	  substr(labels(X)[which(substr(labels(X),0,17)=="relative_density_")][which.max(X[labels(X)[which(substr(labels(X),0,17)=="relative_density_")]])],18,1000)
	}));

	# deproject dataset if needed
	if(!is.na(projection)){results[,c(x,y)]<-rgdal::project(as.matrix(results[,c(x,y)]),proj=projection,inv=TRUE);}

	return(results);
}
