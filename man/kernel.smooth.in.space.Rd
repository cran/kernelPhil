% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2d.smooths.r
\name{kernel.smooth.in.space}
\alias{kernel.smooth.in.space}
\title{Kernel smooth data in space alone}
\usage{
kernel.smooth.in.space(
  dataset,
  dependent.variable = "dependent.variable",
  x = "x",
  y = "y",
  weight = "weight",
  normalise.by,
  data.type = "factor",
  alpha = 0.05,
  margin = 0.1,
  kernel.function = gaussian.kernel,
  adaptive.spatial.bw = TRUE,
  measure.points,
  projection = NA,
  round.up.low.variance = TRUE,
  explicit = TRUE
)
}
\arguments{
\item{dataset}{Dataset to be smoothed as a data.frame.}

\item{dependent.variable}{String name of the single column in dataset with the factor dependent variable (if data.type=="factor") or a vector of column names with numeric counts (if data.type=="count") (defaults to "dependent.variable").}

\item{x}{String name of column containing numeric x co-ordinate (defaults to "x").}

\item{y}{String name of column containing numeric y co-ordinate (defaults to "y").}

\item{weight}{String name of column in the dataset with numeric weights (defaults to "weight").}

\item{normalise.by}{String name of column by which data should be normalised (typically factor with document, speaker or writer ids).}

\item{data.type}{The type of the dependent variable: either "factor", if each row is a token, or "count", if each row is a document, speaker or writer with token counts in separate columns (defaults to "factor").}

\item{alpha}{Numeric alpha for calculating error margins (defaults to 0.05).}

\item{margin}{Numeric desired error margin for calculating spatial bandwidths (defaults to 0.1).}

\item{kernel.function}{The kernel function, one of gaussian.kernel, gaussian.square.kernel, triangular.kernel, square.kernel, or a custom function (defaults to gaussian.kernel).}

\item{adaptive.spatial.bw}{A boolean indicating whether the spatial bandwidth is adaptive (set to achieve margin at every point) or static (set to the average of bandwidths needed to achieve margin at every point).}

\item{measure.points}{A data.frame of spatial points at which estimates are to be made, with two columns with the same names as x,y in dataset; if not supplied, estimates are at the same locations as dataset.}

\item{projection}{The spatial projection as a proj4 string - if given, data will be projected before smoothing and results will be deprojected before returning.}

\item{round.up.low.variance}{Set to TRUE if there are periods of time with extremely low variance (defaults to TRUE).}

\item{explicit}{If TRUE, progress will be reported with a progress bar (defaults to TRUE).}
}
\value{
A data.frame with the smoothed estimates.
}
\description{
This function performs kernel smoothing on a dataset in space alone.
}
\examples{
n=400;
synthesised.data<-data.frame(x=stats::runif(n),y=stats::runif(n),
    year=stats::runif(n,0,sqrt(2)));
synthesised.data$dependent.variable<-unlist(lapply(1:nrow(synthesised.data),
    function(X){
    stats::dist(as.matrix(synthesised.data[c(1,X),1:2]),method =
        "euclidean")<synthesised.data$year[X];
}))
result<-kernelPhil::kernel.smooth.in.space(dataset = synthesised.data);
ggplot2::ggplot(result,ggplot2::aes(x=x,y=y,colour=relative_density_TRUE))+
    ggplot2::geom_point();
}
