#' Kernel smooth data in time alone
#'
#' This function performs kernel smoothing on a dataset in time alone.
#'
#' @param dataset The dataset to be smoothed as a data.frame.
#' @param dependent.variable String name of the column in dataset with the dependent variable (defaults to "dependent.variable"); this column should be numeric or factor.
#' @param time String name of the column in dataset with the time variable (defaults to "year").
#' @param weight String name of column in the dataset with numeric weights (defaults to "weight").
#' @param bandwidth Numeric bandwidth of the kernel function.
#' @param sample.density.threshold Numeric local density of samples below which no estimates will be returned.
#' @param length.out The number of measure points along the time axis (defaults to 1000).
#' @param alpha Numeric alpha for calculating error margins (defaults to 0.05).
#' @param xlabel String label for the x-axis in returned plot (defaults to "year").
#' @param ylabel String label for the y-axis in returned plot.
#' @param greyscale If TRUE, plot will be in greyscale; if "compatible", plot will use a colour spectrum which also goes light>dark; otherwise, will use a non-greyscale-compatible colour scale.
#' @param save.path String path to save plot to (if not given, plot will not be saved).
#' @param measure.times A numeric vector of specific times at which to make estimates; if given, sample.density.threshold and length.out will be ignored.
#' @param kernel.function The kernel function, one of gaussian.kernel, gaussian.square.kernel, triangular.kernel, square.kernel, or a custom function.
#' @return A list containing a data.frame with the smoothed estimates, and a ggplot grob visualising them.
#' @examples
#' n=1000;
#' synthesised.data<-data.frame(x=stats::runif(n),y=stats::runif(n),
#'     year=stats::runif(n,0,sqrt(2)));
#' synthesised.data$dependent.variable<-unlist(lapply(1:nrow(synthesised.data),
#'     function(X){
#'     stats::dist(as.matrix(synthesised.data[c(1,X),1:2]),method =
#'     "euclidean")<synthesised.data$year[X];
#' }))
#' result<-kernelPhil::kernel.smooth.in.time(dataset = synthesised.data,
#'     bandwidth = 0.05,sample.density.threshold = 100);
#' result$plot;
#' @export
kernel.smooth.in.time<-function(dataset,dependent.variable="dependent.variable",time="year",weight="weight",bandwidth=10,sample.density.threshold=3,length.out=1000,alpha=0.05,xlabel="year",ylabel,greyscale="compatible",save.path="",measure.times,kernel.function=gaussian.kernel){
	x<-conf_int_upper<-conf_int_lower<-variant<-conf_int<-NULL; # set up variable names for use later

	# data formats
	if(typeof(dataset[,dependent.variable])=="logical"){dataset[,dependent.variable]<-factor(dataset[,dependent.variable]);}

	# supply weights=1 if missing
	if(!(weight %in% colnames(dataset))){dataset$weight<-1;}
	else if(weight!="weight"){dataset$weight<-dataset[,weight];}

	suppressWarnings(sample_density<-stats::density(dataset[,time],bandwidth=bandwidth,kernel="gaussian",weights=dataset$weight));
	sample_density<-data.frame(x=sample_density$x,y=sample_density$y);
	sample_density<-sample_density[which(sample_density$y>=sample.density.threshold),];
	if(missing(measure.times)){
		kde_range<-seq(from=min(sample_density$x),to=max(sample_density$x),length.out=length.out);
	}else{
		kde_range<-seq(from=min(measure.times),to=max(measure.times),length.out=length.out);
	}
	suppressWarnings(sample_density<-stats::density(dataset[,time],bandwidth=bandwidth,kernel="gaussian",weights=dataset$weight,n=length.out,from=min(sample_density$x),to=max(sample_density$x)));
	if(!missing(measure.times)){
		dd<-stats::approxfun(sample_density$x,sample_density$y);
		sample_density<-data.frame(x=measure.times,y=unlist(lapply(measure.times,function(X) dd(X))));
		kde_range<-measure.times;
	}else{
		sample_density<-data.frame(x=sample_density$x,y=sample_density$y);
	}

	dep_density<-data.frame(x=kde_range,sample_density=sample_density$y);
	if(is.factor(dataset[,dependent.variable])){
		dep_density<-do.call(rbind,lapply(levels(dataset[,dependent.variable]),function(level){
			dataset$current_dep<-as.numeric(as.character(dataset[,dependent.variable])==level);
			current_dep_density<-dep_density;
			current_dep_density$variant<-level;
			current_dep_density[,c("dep_density","sd","samp_si","conf_int_lower","conf_int_upper")]=t(matrix(unlist(lapply(kde_range,function(year){
				product_weights<-kernel.function(as.numeric(year)-dataset[,time],bandwidth=bandwidth)*dataset$weight;
				samp_si=sum(product_weights);
				cmean=sum(dataset$current_dep*product_weights)/samp_si;
				csd=sum(product_weights*abs(dataset$current_dep-cmean))/samp_si;
				cconf_int=-1*((csd/sqrt(samp_si))*stats::qnorm(alpha/2));
				return(data.frame(dep_density=cmean,dep_sd=csd,samp_si=samp_si,conf_int_lower=cmean-cconf_int,conf_int_upper=cconf_int+cmean));
			})),nrow=5));
			return(current_dep_density);
		}));
		if(missing(ylabel)){ylabel="rate";}
		if(greyscale==TRUE){
			colourlist<-grDevices::gray.colors(n=length(levels(dataset[,dependent.variable])),start = 0,end = 0.95);
		}else if(greyscale=="compatible"){
			#colourlist<-c("#ffeae5","#fffc20","#5de100","#0c8160","#362697","#73005c","#4d0506","#0d0d0d");
			#colourlist<-inferno(length(levels(dataset[,dependent.variable])));
			colourlist<-grDevices::colorRampPalette(c("#000000","#29054d","#a60a67","#db1620","#f58931","#ffc66b","#fffd8c"))(length(levels(dataset[,dependent.variable])));
		}else{
			colourlist<-c("#C50038","#09247e","#FFE100","#94f220","#911CF8","#FF5C00","#20E5bF","#cB8800","#E9a890","#04841d","#000000","#aaaaaa")[1:length(levels(dataset[,dependent.variable]))];
		}
		names(colourlist)<-as.character(levels(dataset[,dependent.variable]));
		g=ggplot2::ggplot(data= dep_density,ggplot2::aes(x=x))+
			ggplot2::geom_ribbon(ggplot2::aes(ymax=conf_int_upper,ymin=conf_int_lower,fill=variant),colour=NA,alpha=0.5)+
			ggplot2::geom_line(ggplot2::aes(y=dep_density,colour=variant))+
			ggplot2::scale_fill_manual(aesthetics=c("colour","fill"),values=colourlist)+
			ggplot2::labs(x=xlabel,y=ylabel)+ggplot2::theme_minimal()+
			ggplot2::theme(text=ggplot2::element_text(family="serif"));
	}else if(is.numeric(dataset[,dependent.variable])){
		dep_density[,c("dep_density","dep_sd","samp_si","conf_int")]=t(matrix(unlist(lapply(kde_range,function(year){
			product_weights<-kernel.function(as.numeric(year)-dataset[,time],bandwidth=bandwidth)*dataset$weight;
			samp_si=sum(product_weights);
			cmean=sum(dataset[,dependent.variable]*product_weights)/samp_si;
			csd=sum(product_weights*abs(dataset[,dependent.variable]-cmean))/samp_si;
			cconf_int=-1*((csd/sqrt(samp_si))*stats::qnorm(alpha/2));
			return(data.frame(dep_density=cmean,dep_sd=csd,samp_si=samp_si,conf_int=cconf_int));
		})),nrow=4));

		if(missing(ylabel)){ylabel="degree of change";}
		if(greyscale==TRUE){areacolour="#646263";}else{areacolour="#C50038"}
		g=ggplot2::ggplot(data= dep_density,ggplot2::aes(x=x,y=dep_density))+
			ggplot2::geom_ribbon(ggplot2::aes(ymax=dep_density-conf_int,ymin=dep_density+conf_int),colour=NA,fill=areacolour,alpha=0.4)+
			ggplot2::geom_line()+
			ggplot2::xlab(xlabel)+
			ggplot2::ylab(ylabel)+
			ggplot2::theme_minimal()+
			ggplot2::theme(text=ggplot2::element_text(family="serif"));
	}else{
		stop("dependent variable should be numeric or factor");
	}
	if(save.path!=""){
		grDevices::png(paste0(save.path,".png"),width=1200,height=800);
		print(g+ggplot2::theme(text = ggplot2::element_text(size=24)));
		grDevices::dev.off();
	}
	return(list(results=dep_density,plot=g));
}
#' Identify nearest point in time to a given estimate value
#'
#' This function takes the output of kernel.smooth.in.time and identifies the point in time when the smoothed estimate comes closest to some specific value. This is useful for tasks like identifying the likely midpoint of a change.
#'
#' @param kernel.smooths A list output by kernel.smooth.in.time().
#' @param density The value of the dependent variable for which a time is to be identified.
#' @param variant If the dependent variable was a factor, which level is being examined (do not give a value if dependent variable was numeric).
#' @param n The number of nearest points to be returned (useful if the estimates cross the relevant threshold multiple times, defauls to 1).
#' @param timerange Numeric vector of length two - used to restrict search to a specific time range within the kernel.smooth.in.time(), in the form c(min,max).
#' @return One or more numeric values.
#' @export
nearest.point<-function(kernel.smooths,density,variant,n=1,timerange){
	if(!missing(timerange)){
		kernel.smooths$results<-kernel.smooths$results[which(kernel.smooths$results$x>=timerange[1] & kernel.smooths$results$x<=timerange[2]),];
	}
	if(missing(variant)){
		kernel.smooths<-kernel.smooths$results[order(abs(kernel.smooths$results$dep_density-density)),];
	}else{
		kernel.smooths<-kernel.smooths$results[order(abs(kernel.smooths$results$dep_density-density)),];
		kernel.smooths<-kernel.smooths[which(kernel.smooths$variant==variant),];
	}
    return(kernel.smooths$x[1:n]);
}
