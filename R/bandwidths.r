#' Calculate temporal bandwidths by spatial resolution
#'
#' This function calculates relationships between temporal bandwidth and possible spatial resolution for a given power and suggests minimum possible temporal bandwidth for a given resolution
#'
#' @param dataset The dataset to be smoothed as a data.frame.
#' @param dependent.variable String name of the column in dataset with the dependent variable (defaults to "dependent.variable"); this column should be numeric or factor.
#' @param time String name of the column in dataset with the time variable (defaults to "year").
#' @param weight String name of column in the dataset with numeric weights (defaults to "weight").
#' @param alpha Numeric alpha for calculating error margins (defaults to 0.05).
#' @param margin Numeric desired error margin for calculating spatial bandwidths (defaults to 0.1).
#' @param measure.times A numeric vector of specific times at which to make estimates; if not given, will default to seq(from=min(time),to=max(time),length.out=5).
#' @param temporal.bandwidth.limits Numeric vector of length 2 specifying minimum and maximum temporal bandwidth to be tested (defaults to the range of time*0.01 to the range of time*2).
#' @param temporal.bandwidth.n.levels Number of distinct levels of temporal bandwidth to be tested (defaults to 200).
#' @param minimum.spatial.resolution Numeric minimum spatial resolution.
#' @param summary.plots If TRUE, plots of smoothed sample density, the dependent variable, and variance are returned along with the plot of resolution by bandwidth.
#' @param kernel.function The kernel function, one of gaussian.kernel, gaussian.square.kernel, triangular.kernel, square.kernel, or a custom function (defaults to gaussian.kernel).
#' @return A plot of spatial resolution by temporal bandwidth, along with other summary plots of the data if summary.plots==TRUE.
#' @export
calculate.bandwidths.by.resolution<-function(dataset,dependent.variable="dependent.variable",time="year",weight="weight",alpha=0.05,margin=0.1,measure.times,temporal.bandwidth.limits,temporal.bandwidth.n.levels=200,minimum.spatial.resolution=5,summary.plots=FALSE,kernel.function=gaussian.kernel){
	bw<-resolution<-change_ks<-variance<-weight_n<-NULL; # define variables for use later

	time.range=range(dataset[,time])[2]-range(dataset[,time])[1];
	if(missing(measure.times)){measure.times=seq(from=stats::quantile(dataset[,time],0.15),to=stats::quantile(dataset[,time],0.85),length.out=200);}
	if(missing(temporal.bandwidth.limits)){temporal.bandwidth.limits=c(time.range/100,time.range*2)}
  if(!(weight %in% colnames(dataset))) {dataset[,weight]<-1;}

	# calculating sample density, change and variance over time
  statistics<-data.frame(t(matrix(unlist(lapply(measure.times,function(cyear,bandwidth){
    dataset$kernel_weight<-kernel.function(dataset[,time]-as.numeric(cyear),bandwidth=bandwidth)*dataset[,weight];
    return(c(time=cyear,change_ks=stats::weighted.mean(dataset[,dependent.variable],dataset$kernel_weight),weight_n=sum(dataset$kernel_weight),variance=Hmisc::wtd.var(x=dataset[,dependent.variable],weights=dataset$kernel_weight)));
  },bandwidth=time.range/10)),nrow=4)));
  colnames(statistics)<-c("time","change_ks","weight_n","variance");
  change_plot<-ggplot2::ggplot(statistics,ggplot2::aes(x=time,y=change_ks))+ggplot2::geom_line(size=1.2)+ggplot2::ylab("change")+ggplot2::theme_minimal()+ggplot2::theme(text=ggplot2::element_text(family="serif",size=14));
  variance_plot<-ggplot2::ggplot(statistics,ggplot2::aes(x=time,y=variance))+ggplot2::geom_line(size=1.2)+ggplot2::ylab("variance")+ggplot2::theme_minimal()+ggplot2::theme(text=ggplot2::element_text(family="serif",size=14));
  n_plot<-ggplot2::ggplot(statistics,ggplot2::aes(x=time,y=weight_n))+ggplot2::geom_line(size=1.2)+ggplot2::ylim(c(min(statistics$weight_n)*0.75,max(statistics$weight_n)*1.25))+ggplot2::ylab("sum weight")+ggplot2::theme_minimal()+ggplot2::theme(text=ggplot2::element_text(family="serif",size=14));

	# calculating resolutions
  resolutions<-expand.grid(time=measure.times,bw=seq(from=temporal.bandwidth.limits[1],to=temporal.bandwidth.limits[2],length.out=temporal.bandwidth.n.levels));
	message("calculating required sample sizes\n");
  resolutions[,c("n","k")]<-data.frame(t(matrix(unlist(pbapply::pbapply(resolutions,1,function(X){
    kernel_weights<-kernel.function(dataset[,time]-as.numeric(X["time"]),bandwidth=as.numeric(X["bw"]))*dataset[,weight];
    return(c(sum(kernel_weights),sample.size(alpha,margin,dataset[,dependent.variable],kernel_weights)));
  })),nrow=2)));
  resolutions$resolution<-resolutions$n/resolutions$k;
	message("calculating best resolutions by time\n");
  bw_at_res_min<-data.frame(time=measure.times,bw=unlist(pbapply::pblapply(levels(factor(as.character(resolutions$time))),function(time){
    current_res<-resolutions[which(as.character(resolutions$time)==time & resolutions$resolution>=minimum.spatial.resolution),];
    current_res<-current_res[order(current_res$bw),];
    return(current_res$bw[1])
  })));
  colourlist<-grDevices::colorRampPalette(c("#000000","#29054d","#a60a67","#db1620","#f58931","#ffc66b","#fffd8c"));
  res_plot<-ggplot2::ggplot(resolutions[which(resolutions$resolution<=stats::quantile(resolutions$resolution,0.99,na.rm=TRUE)),],ggplot2::aes(x=time,y=bw,fill=resolution))+
    ggplot2::geom_tile()+
    ggplot2::geom_vline(xintercept=mean(bw_at_res_min$time[which(bw_at_res_min$bw==max(bw_at_res_min$bw))]),linetype="dotted",colour="grey")+
    ggplot2::geom_hline(yintercept=max(bw_at_res_min$bw),linetype="dotted",colour="grey")+
    ggplot2::scale_fill_gradientn(colours=colourlist(50))+
    ggplot2::theme_minimal()+ggplot2::theme(text=ggplot2::element_text(family="serif",size=14));
  if(is.na(max(bw_at_res_min$bw))){
    stop("it is not possible to achieve this spatial resolution with these settings\n")
  }else{
    message(paste0("temporal bandwidth needed is ",max(bw_at_res_min$bw),", choke time is ",mean(bw_at_res_min$time[which(bw_at_res_min$bw==max(bw_at_res_min$bw))]),"\n"));
  }
  if(!summary.plots){return(res_plot);}
  else{return(gridExtra::grid.arrange(change_plot,n_plot,variance_plot,res_plot));}
}
#' Calculate temporal bandwidths by separated points
#'
#' This function calculates relationships between temporal bandwidth and spatial bandwidths at a series of specified points for a given power and suggests minimum possible temporal bandwidth such that bandwidths at those points are never greater than 2.2365*the distance to the nearest point (for gaussian kernels) or 2*that distance (for other kernels)
#'
#' @param dataset The dataset to be smoothed as a data.frame.
#' @param dependent.variable String name of the column in dataset with the dependent variable (defaults to "dependent.variable"); this column should be numeric or factor.
#' @param x String name of column containing numeric x co-ordinate (defaults to "x").
#' @param y String name of column containing numeric y co-ordinate (defaults to "y").
#' @param time String name of the column in dataset with the time variable (defaults to "year").
#' @param weight String name of column in the dataset with numeric weights (defaults to "weight").
#' @param alpha Numeric alpha for calculating error margins (defaults to 0.05).
#' @param margin Numeric desired error margin for calculating spatial bandwidths (defaults to 0.1).
#' @param separated.points Data.frame containing two columns same names as x,y in dataset with x and y coordinates of points to be kept separate.
#' @param measure.times A numeric vector of specific times at which to make estimates; if not given, will default to seq(from=min(time),to=max(time),length.out=5).
#' @param temporal.bandwidth.limits A numeric vector of length 2 specifying minimum and maximum temporal bandwidth to be tested (defaults to the range of time*0.01 to the range of time*2).
#' @param temporal.bandwidth.n.levels Number of distinct levels of temporal bandwidth to be tested (defaults to 200).
#' @param kernel.function The kernel function, one of gaussian.kernel, gaussian.square.kernel, triangular.kernel, square.kernel, or a custom function (defaults to gaussian.kernel).
#' @param projection A spatial projection as a proj4 string - if given, data will be projected before smoothing and results will be deprojected before returning.
#' @param include.visualisation If TRUE, will return a ggplot visualisation.
#' @param separated.points.labels String vector of the names of the separated points (used in the visualisation).
#' @param round.up.low.variance Set to TRUE if there are periods of time with extremely low variance.
#' @return A list with suggested bandwidth and the choke point and time, plus a visualisation of bandwidths and resolutions if include.visualisation==TRUE.
#' @export
calculate.bandwidths.by.separated.points<-function(dataset,dependent.variable="dependent.variable",x="x",y="y",time="year",weight="weight",alpha=0.05,margin=0.1,separated.points,measure.times,temporal.bandwidth.limits,temporal.bandwidth.n.levels=200,kernel.function=gaussian.kernel,projection=NA,include.visualisation=FALSE,separated.points.labels,round.up.low.variance=FALSE){
	change_ks<-variance<-weight_n<-prediction.time<-temporal.bandwidth<-spatial.bandwidth<-..level..<-line.num<-separated.point<-NULL; # define variables for use later

	time.range=range(dataset[,time])[2]-range(dataset[,time])[1];
	if(missing(measure.times)){measure.times=seq(from=stats::quantile(dataset[,time],0.15),to=stats::quantile(dataset[,time],0.85),length.out=200);}
	if(missing(temporal.bandwidth.limits)){temporal.bandwidth.limits=c(time.range/200,time.range/2)}
	if(temporal.bandwidth.limits[1]<max(diff(dataset[,time][order(dataset[,time])]))){temporal.bandwidth.limits[1]=max(diff(dataset[,time][order(dataset[,time])]));}		# if lower, bottom of temporal bandwidth search range is the maximum temporal distance between any two points (lower than this and we start to get weird artifacts in calculating k because weighted variance will be close to 0 for certain points in time)
	if(!(weight %in% colnames(dataset))){dataset[,weight]<-1;} # if not passed weights, weight everything 1
	if(missing(separated.points.labels)){separated.points.labels=paste0("point ",1:nrow(separated.points));}			# if separated points aren't named, number them
	# project dataset if needed
	if(!is.na(projection)){
		dataset[,c(x,y)]<-rgdal::project(as.matrix(dataset[,c(x,y)]),proj=projection);
		separated.points<-rgdal::project(as.matrix(separated.points[,c(x,y)]),proj=projection);
	}

	# calculating sample density, change and variance over time
	statistics<-data.frame(t(matrix(unlist(lapply(measure.times,function(cyear,bandwidth){
		dataset$kernel_weight<-kernel.function(dataset[,time]-as.numeric(cyear),bandwidth=bandwidth)*dataset[,weight];
		return(c(time=cyear,change_ks=stats::weighted.mean(dataset[,dependent.variable],dataset$kernel_weight),weight_n=sum(dataset$kernel_weight),variance=Hmisc::wtd.var(x=dataset[,dependent.variable],weights=dataset$kernel_weight)));
	},bandwidth=time.range/10)),nrow=4)));
	colnames(statistics)<-c("time","change_ks","weight_n","variance");
	change_plot<-ggplot2::ggplot(statistics,ggplot2::aes(x=time,y=change_ks))+ggplot2::geom_line(size=1.2)+ggplot2::ylab("change")+ggplot2::theme_minimal()+ggplot2::theme(text=ggplot2::element_text(family="serif",size=14));
	variance_plot<-ggplot2::ggplot(statistics,ggplot2::aes(x=time,y=variance))+ggplot2::geom_line(size=1.2)+ggplot2::ylab("variance")+ggplot2::theme_minimal()+ggplot2::theme(text=ggplot2::element_text(family="serif",size=14));
	n_plot<-ggplot2::ggplot(statistics,ggplot2::aes(x=time,y=weight_n))+ggplot2::geom_line(size=1.2)+ggplot2::ylim(c(min(statistics$weight_n)*0.75,max(statistics$weight_n)*1.25))+ggplot2::ylab("sum weight")+ggplot2::theme_minimal()+ggplot2::theme(text=ggplot2::element_text(family="serif",size=14));


	dm<-wordspace::dist.matrix(M2 =  as.matrix(separated.points),M = as.matrix(dataset[,c(x,y)]),method = "euclidean"); 			# construct a distance matrix between the separated points and all other points, for applying the spatial kernel
	separated.points.dm<-wordspace::dist.matrix(as.matrix(separated.points),method = "euclidean"); 								# construct a distance matrix just between the separated points for the purpose of identifying where overlaps start to happen
	separated.points.dm[which(separated.points.dm==0)]<-NA;																# don't include comparisons between points and themselves
	nearest.points<-unlist(apply(separated.points.dm,2,min,na.rm=TRUE)); 													# identify the distance to the nearest neighbour for each separated point
	max_spatial_bandwidths<-nearest.points/ifelse(identical(gaussian.kernel,kernel.function),2.2365,ifelse(identical(gaussian.square.kernel,kernel.function),2.2365,2)); 																		# identify the bw size at which each point would start overlapping with a neighbour
	spatial.bandwidth.search.range<-c(min(max_spatial_bandwidths)/1000,max(max_spatial_bandwidths)*10); 					# search range for spatial bandwidths is 0.1% of max distance to nearest neighbouring separated point to 1000% of max distance to nearest neighbouring separated point

	if(round.up.low.variance){sample.size.func<-sample.size.floored.variance;}else{sample.size.func<-sample.size;}
	spatial.bandwidths.by.time.and.temporal.bandwidth.matrices<-list();CLs<-list();																		# set up variables to contain results for each separated point
	maxima.by.point<-data.frame(separated.point=1:nrow(separated.points),x=separated.points[,x],y=separated.points[,y],maximum.temporal.bandwidth=NA,location.of.maximum=NA);
	maxima.by.point[,c("maximum.temporal.bandwidth","location.of.maximum")]<-t(matrix(unlist(lapply(1:nrow(separated.points),function(i){				# for each separated point...
		message(paste0("calculating minima for ",separated.points.labels[i],"\n"));
		spatial.bandwidths.by.time.and.temporal.bandwidth.matrices[[i]]<<-t(matrix(unlist(pbapply::pblapply(measure.times,function(prediction.time){			# at each point in time...
			unlist(lapply(seq(from=temporal.bandwidth.limits[1],to=temporal.bandwidth.limits[2],length.out=temporal.bandwidth.n.levels),function(current.temporal.bandwidth){	# for every possible temporal bandwidth...
				temporal.kernel.weights<-gaussian.kernel(dataset[,time]-prediction.time,bandwidth=current.temporal.bandwidth)*dataset[,weight];			# apply temporal kernel
				k<-sample.size.func(alpha=alpha,margin=margin,dependent.variable=dataset[,dependent.variable],weights=temporal.kernel.weights);											# calculate k based on these weights
				if(is.na(k)){return(NA);}																												# (if k is na, this indicates that there is 1 or fewer effective observations)
				else if(sum(temporal.kernel.weights)<k){return(NA);}																					# (if sum of weights<k, no spatial bandwidth can possibly achieve the needed sample size)
				else{
					return(stats::optimise(function(spatial.bandwidth){																						# and identify optimal spatial bandwidth such that sum(weights)=k after applying the spatial kernel at this point
						abs(k-sum(kernel.function(dm[,i],as.numeric(spatial.bandwidth))*temporal.kernel.weights));
					},lower=spatial.bandwidth.search.range[1],upper=spatial.bandwidth.search.range[2],tol=0.01,maximum = FALSE)$minimum);
				}
			}));
		})),nrow=temporal.bandwidth.n.levels));
		if(!all(apply(spatial.bandwidths.by.time.and.temporal.bandwidth.matrices[[i]],1,function(X) min(as.numeric(X),na.rm=TRUE)<=max_spatial_bandwidths[i]))){		# check if the necessary spatial bandwidths never fall below the overlap threshold
			stop(paste0("no temporal bandwidth within range achieves desired power for ",separated.points.labels[i],"; rerun with higher temporal bandwidth range or lower power"));
		}else if(max_spatial_bandwidths[i]>max(spatial.bandwidths.by.time.and.temporal.bandwidth.matrices[[i]],na.rm=TRUE)){											# check if the necessary spatial bandwidths are always above all tested values
			warning(paste0("all temporal bandwidths within range achieve desired power for ",separated.points.labels[i]));
			return(c(NA,NA));
		}else{																																			# if everything is fine, then calculate the contour line for the maximum acceptable spatial bandwidth for this point
			CLs[[i]]<<-grDevices::contourLines(x = measure.times,y=seq(from=temporal.bandwidth.limits[1],to=temporal.bandwidth.limits[2],length.out=temporal.bandwidth.n.levels),z = spatial.bandwidths.by.time.and.temporal.bandwidth.matrices[[i]],levels = max_spatial_bandwidths[i]);
			current.combined.CL<-do.call("rbind",lapply(1:length(CLs[[i]]),function(j){
				data.frame(separated.point=separated.points.labels[i],time=CLs[[i]][[j]]$x,temporal.bandwidth=CLs[[i]][[j]]$y);
			}));
			current.combined.CL<-current.combined.CL[which(apply(current.combined.CL,1,function(crow){
				if(nrow(current.combined.CL[which(round(current.combined.CL$time,3)==round(as.numeric(crow["time"]),3)),])==1){
					return(TRUE);
				}else if(round(as.numeric(crow["temporal.bandwidth"]),3)==round(min(current.combined.CL$temporal.bandwidth[which(round(current.combined.CL$time,3)==round(as.numeric(crow["time"]),3))]),3)){
					return(TRUE);
				}else{
					return(FALSE);
				}
			})),];
			return(c(max(current.combined.CL$temporal.bandwidth),current.combined.CL$time[which.max(current.combined.CL$temporal.bandwidth)]));
		}
	})),nrow=2));
	message(paste0("temporal bandwidth needed is ",max(maxima.by.point$maximum.temporal.bandwidth),", choke time is ",maxima.by.point$location.of.maximum[which.max(maxima.by.point$maximum.temporal.bandwidth)]," at point ",separated.points.labels[maxima.by.point$separated.point[which.max(maxima.by.point$maximum.temporal.bandwidth)]],"\n"));
	if(!include.visualisation){
		return(list(suggested.bandwidth=max(maxima.by.point$maximum.temporal.bandwidth),choke.time=maxima.by.point$location.of.maximum[which.max(maxima.by.point$maximum.temporal.bandwidth)],choke.point=separated.points.labels[maxima.by.point$separated.point[which.max(maxima.by.point$maximum.temporal.bandwidth)]]));
	}else{
		colourlist<-grDevices::colorRampPalette(c("#000000","#29054d","#a60a67","#db1620","#f58931","#ffc66b","#fffd8c","#ffffff"));
		plots<-list();
		for(i in 1:nrow(separated.points)){
			cutoff=max_spatial_bandwidths[i]*5;
			colnames(spatial.bandwidths.by.time.and.temporal.bandwidth.matrices[[i]])<-seq(from=temporal.bandwidth.limits[1],to=temporal.bandwidth.limits[2],length.out=temporal.bandwidth.n.levels);
			rownames(spatial.bandwidths.by.time.and.temporal.bandwidth.matrices[[i]])<-measure.times;
			ylims=c(min(unlist(sapply(CLs[[i]],function(CL) CL$y)))*0.9,max(unlist(sapply(CLs[[i]],function(CL) CL$y)))*1.4);
			current_df<-reshape2::melt(spatial.bandwidths.by.time.and.temporal.bandwidth.matrices[[i]]);
			colnames(current_df)<-c("prediction.time","temporal.bandwidth","spatial.bandwidth");
			current_df$spatial.bandwidth[which(current_df$spatial.bandwidth>cutoff)]=cutoff;
			plots[[i]]<-ggplot2::ggplot() +
				ggplot2::geom_raster(data=current_df[which(current_df$temporal.bandwidth>ylims[1] & current_df$temporal.bandwidth<ylims[2]),], ggplot2::aes(x=prediction.time, y=temporal.bandwidth,fill=spatial.bandwidth), show.legend = TRUE) +
				ggplot2::scale_fill_gradientn(limits=range(current_df[which(current_df$temporal.bandwidth>ylims[1] & current_df$temporal.bandwidth<ylims[2]),]$spatial.bandwidth), colours=colourlist(50)[50:1],na.value="black",name="spatial bandwidth") +
				ggplot2::geom_contour(data = current_df[which(current_df$temporal.bandwidth>ylims[1] & current_df$temporal.bandwidth<ylims[2]),],bins=5,ggplot2::aes(z=spatial.bandwidth,x=prediction.time, y=temporal.bandwidth,colour = ..level..)) +
				ggplot2::scale_colour_gradient(guide = 'none',high="#000000",low="#000000");
			for(j in 1:length(CLs[[i]])){
				plots[[i]]<-plots[[i]]+
					ggplot2::geom_path(data=data.frame(x=CLs[[i]][[j]]$x,y=CLs[[i]][[j]]$y),ggplot2::aes(x=x,y=y),linetype="dotted");
			}
			plots[[i]]<-plots[[i]]+
				ggplot2::xlab("time")+ggplot2::ylab("temporal bandwidth")+
				ggplot2::ggtitle(separated.points.labels[[i]])+
				ggplot2::ylim(ylims)+
				ggplot2::theme_minimal()+ggplot2::theme(text=ggplot2::element_text(family="serif",size=14));
			plots[[i]]<-directlabels::direct.label(plots[[i]],list("far.from.others.borders", "calc.boxes", "enlarge.box", hjust = 1, vjust = 1, box.color = NA, fill = "transparent", "draw.rects",colour="black",fontfamily="serif"));
		}
		# generate plot of all contours
		colourlist<-grDevices::colorRampPalette(c("#000000","#29054d","#a60a67","#db1620","#f58931","#ffc66b","#fffd8c"));
		plots[[nrow(separated.points)+1]]<-ggplot2::ggplot(data=do.call("rbind",lapply(1:length(CLs),function(i){
			do.call("rbind",lapply(1:length(CLs[[i]]),function(j){
				data.frame(separated.point=separated.points.labels[i],line.num=paste0(i,".",j),time=CLs[[i]][[j]]$x,temporal.bandwidth=CLs[[i]][[j]]$y);
			}))
		})),ggplot2::aes(x=time,y=temporal.bandwidth,group=line.num,colour=separated.point))+
			ggplot2::scale_colour_manual(values=colourlist(length(CLs)),name="point")+
			ggplot2::geom_path(size=1.2)+
			ggplot2::ggtitle("minimum temporal bandwidth by time and separated point")+
			ggplot2::ylab("temporal bandwidth")+
			ggplot2::theme_minimal()+
			ggplot2::theme(text=ggplot2::element_text(family="serif",size=12));

		layout_matrx<-matrix(if(length(plots)%%2==0) 1:length(plots) else c(1:length(plots),length(plots)),ncol=2);
		return(list(suggested.bandwidth=max(maxima.by.point$maximum.temporal.bandwidth),choke.time=maxima.by.point$location.of.maximum[which.max(maxima.by.point$maximum.temporal.bandwidth)],choke.point=separated.points.labels[maxima.by.point$separated.point[which.max(maxima.by.point$maximum.temporal.bandwidth)]],individual.plots=plots,composite.plot=function(){gridExtra::grid.arrange(grobs=plots,layout_matrix=layout_matrx)}));
	}
}
