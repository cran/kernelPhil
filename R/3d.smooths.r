#' Kernel smooth data in space and time
#'
#' This function performs kernel smoothing on a dataset in time and space. A static temporal kernel is applied first, and then an (optionally) adaptive spatial kernel on this weighted data.
#'
#' @param dataset The dataset to be smoothed as a data.frame.
#' @param dependent.variable String name of the single column in dataset with the factor dependent variable (if data.type=="factor") or a vector of column names with numeric counts (if data.type=="count") (defaults to "dependent.variable").
#' @param x String name of column containing numeric x co-ordinate (defaults to "x").
#' @param y String name of column containing numeric y co-ordinate (defaults to "y").
#' @param time String name of the column in dataset with the time variable (defaults to "year").
#' @param weight String name of column in the dataset with numeric weights (defaults to "weight").
#' @param normalise.by String name of column by which data should be normalised (typically factor with document, speaker or writer ids).
#' @param data.type The type of the dependent variable as a string: either "factor", if each row is a token, or "count", if each row is a document, speaker or writer with token counts in separate columns (defaults to "factor").
#' @param alpha Numeric alpha for calculating error margins (defaults to 0.05).
#' @param margin Numeric desired error margin for calculating spatial bandwidths (defaults to 0.1).
#' @param kernel.function The kernel function, one of gaussian.kernel, gaussian.square.kernel, triangular.kernel, square.kernel, or a custom function (defaults to gaussian.kernel).
#' @param adaptive.spatial.bw Boolean indicating whether the spatial bandwidth is adaptive (set to achieve margin at every point) or static (set to the average of bandwidths needed to achieve margin at every point).
#' @param temporal.bandwidth Numeric bandwidth of the (gaussian) temporal kernel.
#' @param measure.points A data.frame of spatial points at which estimates are to be made, with two columns with the same names as x,y in dataset; if not supplied, estimates are at the same locations as dataset.
#' @param measure.times A numeric vector of specific times at which to make estimates; if not given, will default to seq(from=min(time),to=max(time),length.out=5).
#' @param projection Spatial projection as a proj4 string - if given, data will be projected before smoothing and results will be deprojected before returning.
#' @param explicit If TRUE, progress will be reported with a progress bar (defaults to TRUE).
#' @return A list containing the parameters and a data.frame with the smoothed estimates.
#' @examples
#' n=200;
#' synthesised.data<-data.frame(x=stats::runif(n),y=stats::runif(n),
#'     year=stats::runif(n,0,sqrt(2)));
#' synthesised.data$dependent.variable<-unlist(lapply(1:nrow(synthesised.data),
#'     function(X){
#'     stats::dist(as.matrix(synthesised.data[c(1,X),1:2]),method =
#'         "euclidean")<synthesised.data$year[X];
#' }));
#' result<-kernelPhil::kernel.smooth.in.space.and.time(dataset =
#'     synthesised.data,temporal.bandwidth = 0.25,measure.times =
#'     seq(from=-0.05,to=1.15,length.out=4),alpha = 0.15,margin = 0.2);
#' gridExtra::grid.arrange(ggplot2::ggplot(result$results[[1]],
#'     ggplot2::aes(x=x,y=y,colour=relative_density_TRUE))+
#'     ggplot2::geom_point(),ggplot2::ggplot(result$results[[2]],
#'     ggplot2::aes(x=x,y=y,colour=relative_density_TRUE))+
#'     ggplot2::geom_point(),ggplot2::ggplot(result$results[[3]],
#'     ggplot2::aes(x=x,y=y,colour=relative_density_TRUE))+
#'     ggplot2::geom_point(),ggplot2::ggplot(result$results[[4]],
#'     ggplot2::aes(x=x,y=y,colour=relative_density_TRUE))+
#'     ggplot2::geom_point());
#' @export
kernel.smooth.in.space.and.time<-function(dataset,dependent.variable="dependent.variable",x="x",y="y",time="year",weight="weight",normalise.by,data.type="factor",alpha=0.05,margin=0.1,kernel.function=gaussian.kernel,adaptive.spatial.bw=TRUE,temporal.bandwidth,measure.points,measure.times,projection=NA,explicit=TRUE){
	# make sure variables present in the right formats
	if(!(weight %in% colnames(dataset))){
		dataset[,weight]<-1;
	}

	# if a column name is passed with "normalise.by", then divide weight by count of factor levels (used for normalising weights by text id, speaker id, etc.)
	if(!missing(normalise.by)){
		dataset[,weight]<-apply(dataset,1,function(X) return(as.numeric(X[weight])/nrow(dataset[which(dataset[,normalise.by]==X[normalise.by]),])));
	}

	# if no separate prediction points passed, use the coordinates of the dataset as prediction points
	if(missing(measure.points)){
		measure.points=unique(dataset[,c(x,y)]);
	}

	# if no separate prediction years passed, generate a reasonable sequence of years
	if(missing(measure.times)){
		measure.times<-seq(from=min(dataset[,time]),to=max(dataset[,time]),length.out=5);
	}

	apply.kernel.smooth.in.space.by.time<-function(prediction.time){
		dataset$kernel.weight<-gaussian.kernel(dataset[,time]-prediction.time,bandwidth=temporal.bandwidth)*dataset[,weight];
		return(kernel.smooth.in.space(x=x,y=y,dependent.variable=dependent.variable,dataset=dataset,alpha=alpha,margin=margin,adaptive.spatial.bw=adaptive.spatial.bw,measure.points=measure.points,weight="kernel.weight",data.type=data.type,projection=projection,explicit=FALSE));
	}
	if(explicit){
		message("calculating smooths by point in time\n");
		result<-pbapply::pblapply(measure.times,apply.kernel.smooth.in.space.by.time);
	}else{
		result<-lapply(measure.times,apply.kernel.smooth.in.space.by.time);
	}
	return(list(parameters=list(alpha=alpha,margin=margin,temporal.bandwidth=temporal.bandwidth,measure.times=measure.times),results=result));
}
#' Kernel smooth data in space and time, returning specific error margins at each point
#'
#' This function performs kernel smoothing on a dataset in time and space. A static temporal kernel is applied first, and then an (optionally) adaptive spatial kernel on this weighted data. Note that this is the same as kernel.smooth.in.space.and.time() except that it returns specific error margins with every estimate and is *much* slower.
#'
#' @param dataset The dataset to be smoothed as a data.frame.
#' @param dependent.variable String name of the single column in dataset with the factor dependent variable (if data.type=="factor") or a vector of column names with numeric counts (if data.type=="count") (defaults to "dependent.variable").
#' @param x String name of column containing numeric x co-ordinate (defaults to "x").
#' @param y String name of column containing numeric y co-ordinate (defaults to "y").
#' @param time String name of the column in dataset with the time variable (defaults to "year").
#' @param weight String name of column in the dataset with numeric weights (defaults to "weight").
#' @param normalise.by String name of column by which data should be normalised (typically factor with document, speaker or writer ids).
#' @param data.type The type of the dependent variable as a string: either "factor", if each row is a token, or "count", if each row is a document, speaker or writer with token counts in separate columns (defaults to "factor").
#' @param alpha Numeric alpha for calculating error margins (defaults to 0.05).
#' @param margin Numeric desired error margin for calculating spatial bandwidths.
#' @param kernel.function The kernel function, one of gaussian.kernel, gaussian.square.kernel, triangular.kernel, square.kernel, or a custom function (defaults to gaussian.kernel).
#' @param adaptive.spatial.bw Boolean indicating whether the spatial bandwidth is adaptive (set to achieve margin at every point) or static (set to the average of bandwidths needed to achieve margin at every point).
#' @param temporal.bandwidth Numeric bandwidth of the (gaussian) temporal kernel.
#' @param measure.points A data.frame of spatial points at which estimates are to be made, with two columns with the same names as x,y in dataset; if not supplied, estimates are at the same locations as dataset.
#' @param measure.times A numeric vector of specific times at which to make estimates; if not given, will default to seq(from=min(time),to=max(time),length.out=5).
#' @param projection The spatial projection as a proj4 string - if given, data will be projected before smoothing and results will be deprojected before returning.
#' @param explicit If TRUE, progress will be reported with a progress bar (defaults to TRUE).
#' @return A list containing the parameters and a data.frame with the smoothed estimates.
#' @examples
#' \donttest{n=200;
#' synthesised.data<-data.frame(x=stats::runif(n),y=stats::runif(n),
#'     year=stats::runif(n,0,sqrt(2)));
#' synthesised.data$dependent.variable<-unlist(lapply(1:nrow(synthesised.data),
#'     function(X){
#'     stats::dist(as.matrix(synthesised.data[c(1,X),1:2]),method =
#'         "euclidean")<synthesised.data$year[X];
#' }));
#' result<-kernelPhil::kernel.smooth.in.space.and.time.with.margins(dataset =
#'     synthesised.data,temporal.bandwidth = 0.2,measure.times =
#'     seq(from=0.15,to=0.85,length.out=2),alpha=0.4,margin=0.2);
#' gridExtra::grid.arrange(ggplot2::ggplot(result$results[[1]],
#'     ggplot2::aes(x=x,y=y,colour=relative_density_TRUE))+
#'     ggplot2::geom_point(),ggplot2::ggplot(result$results[[2]],
#'     ggplot2::aes(x=x,y=y,colour=relative_density_TRUE))+
#'     ggplot2::geom_point())}
#' @export
kernel.smooth.in.space.and.time.with.margins<-function(dataset,dependent.variable="dependent.variable",x="x",y="y",time="year",weight="weight",normalise.by,data.type="factor",alpha=0.05,margin=0.1,kernel.function=gaussian.kernel,adaptive.spatial.bw=TRUE,temporal.bandwidth,measure.points,measure.times,projection=NA,explicit=TRUE){
	# make sure variables present in the right formats
	if(!(weight%in%colnames(dataset))){
		dataset$weight<-1;
	}else{
		dataset$weight<-dataset[,weight];
	}

	# if a column name is passed with "normalise.by", then divide weight by count of factor levels (used for normalising weights by text id, speaker id, etc.)
	if(!missing(normalise.by)){
		dataset[,weight]<-apply(dataset,1,function(X)return(as.numeric(X[weight])/nrow(dataset[which(dataset[,normalise.by]==X[normalise.by]),])));
	}

	# if no separate prediction points passed, use the coordinates of the dataset as prediction points
	if(missing(measure.points)){
		measure.points=unique(dataset[,c(x,y)]);
	}

	# if no separate prediction years passed, generate a reasonable sequence of years
	if(missing(measure.times)){
		measure.times<-seq(from=min(dataset[,time]),to=max(dataset[,time]),length.out=5);
	}

	apply.kernel.smooth.in.space.with.margins.by.time<-function(prediction.time){
		dataset$kernel.weight<-gaussian.kernel(dataset[,time]-prediction.time,bandwidth=temporal.bandwidth)*dataset$weight;
		return(kernel.smooth.in.space.with.margins(x=x,y=y,dependent.variable=dependent.variable,dataset=dataset,alpha=alpha,margin=margin,adaptive.spatial.bw=adaptive.spatial.bw,measure.points=measure.points,weight="kernel.weight",data.type=data.type,explicit=FALSE,projection=projection));
	}
	if(explicit){
		message("calculating smooths by point in time \n");
		result<-pbapply::pblapply(measure.times,apply.kernel.smooth.in.space.with.margins.by.time);
	}else{
		result<-lapply(measure.times,apply.kernel.smooth.in.space.with.margins.by.time);
	}

	return(list(parameters=list(alpha=alpha,margin=margin,temporal.bandwidth=temporal.bandwidth,measure.times=measure.times),results=result));
}
#' Save kernel smooths in space and time
#'
#' Saves the output of kernel.smooth.in.space.and.time() or kernel.smooth.in.space.and.time.with.margins() to a directory
#'
#' @param kernel.smooth A list output of kernel.smooth.in.space.and.time() or kernel.smooth.in.space.and.time.with.margins().
#' @param location String location on the disk to save output to.
#' @param variable.name String name of the variable (used in filenames).
#' @return No returned value, called to save data to disk.
#' @export
save.kernel.smooths<-function(kernel.smooth,location,variable.name){
	dir_name=paste0(variable.name,"_smoothed_alpha",kernel.smooth$parameters$alpha,"_margin",kernel.smooth$parameters$margin,"_tbw",kernel.smooth$parameters$temporal.bandwidth,"_",round(min(kernel.smooth$parameters$measure.times),1),"_",round(max(kernel.smooth$parameters$measure.times),1));
	dir.create(file.path(location, dir_name), showWarnings = FALSE);
	location=paste0(location,"/",dir_name);
	for(i in 1:length(kernel.smooth$parameters$measure.times)){

		kernel.smooth$results[[i]]$best<-apply(kernel.smooth$results[[i]],1,function(X){
			return(substr(names(which.max(X[names(X)[substr(names(X),1,17)=="relative_density_"]])),18,100));
		});
		kernel.smooth$results[[i]]$relative_density_best<-apply(kernel.smooth$results[[i]],1,function(X) X[paste0("relative_density_",X["best"])]);

		utils::write.csv(kernel.smooth$results[[i]],paste0(location,"/",variable.name,"_smoothed_alpha",kernel.smooth$parameters$alpha,"_margin",kernel.smooth$parameters$margin,"_tbw",kernel.smooth$parameters$temporal.bandwidth,"_",round(kernel.smooth$parameters$measure.times[i],1),".csv"),fileEncoding="UTF-8");
	}
	kernel.smooth$parameters$measure.times=paste0(kernel.smooth$parameters$measure.times,",");
	utils::write.table(kernel.smooth$parameters,paste0(location,"/",variable.name,"_smoothed_alpha",kernel.smooth$parameters$alpha,"_margin",kernel.smooth$parameters$margin,"_tbw",kernel.smooth$parameters$temporal.bandwidth,"_parameters.txt"));
}
#' Load kernel smooths
#'
#' Loads the output of kernel.smooth.in.space.and.time() or kernel.smooth.in.space.and.time.with.margins() previously saved with save.kernel.smooths()
#'
#' @param location String location to which results were saved.
#' @return A list containing the parameters and a data.frame with the smoothed estimates (same structure as returned by kernel.smooth.in.space.and.time() and kernel.smooth.in.space.and.time.with.margins()).
#' @export
load.kernel.smooths<-function(location){
	kernel.smooth<-list(parameters=list(),results=list());
	filenames<-list.files(location);
	kernel.smooth<-list(parameters=list(alpha=NA,margin=NA,measure.times=NA),results=lapply(filenames[grep("_[0-9\\.]+\\.csv$",filenames)],function(filename) utils::read.csv(paste0(location,"/",filename),header=TRUE,encoding="UTF-8",row.names=1)));
	kernel.smooth$parameters<-as.list(utils::read.table(paste0(location,"/",filenames[grep("_parameters\\.txt",filenames)][1]),header=TRUE));
	kernel.smooth$parameters$measure.times<-as.numeric(strsplit(kernel.smooth$parameters$measure.times,",")[[1]]);
	return(kernel.smooth);
}
