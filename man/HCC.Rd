\name{HCC}
\alias{HCC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hierarchical Cluster Coherence Method}
\description{
Compute the hierarchical cluster coherence method for a set of time series X.
}
\usage{
HCC(X,Clustfreq=NULL,freq=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Matrix of time series, the series should be located by column.
}
\item{Clustfreq}{
Frequency choosen to perform the clustering. If a interval in provide then clustering is performed based on integrated cluster coherence on the correponding frequency band. If NULL then the disimilarity for all frequencies is returned.
}
  \item{freq}{
Sampling Frequency. Default value is 1.
}
}
\details{
Let X_1,X_2,...,X_N be the signal for each channel with length T and sampling frequency F_s. The procedure starts with N clusters, one for each individual channel.

1) Estimate the coherence matrix C(\omega) at frequency \omega.

2) Compute the initial dissimilarity matrix at band \Omega_{12}.

3) Find the two clusters with the lowest dissimilarity and save this value as a charac- teristic.

4) Merge the signals of the two most similar clusters, reduce the number of clusters by one, i.e., ki = ki−1 − 1, and increase i accordingly , i.e., i = i + 1.

5) Compute the dissimilarity between the new cluster and the existing ones.

6) Repeat steps 2-5 until there is only one cluster left.
}
\value{
A HCC object with the following variables:

Diss.Matrix = Initial dissimilarity matrix.

min.value = trayectory of the minimum value.

Groups = list with the groupping structure at each step.
}
\references{
Euan, C., Sun, Y. and Ombao, H. (2017) "Coherence-based Time Series Clustering for Brain Connectivity Visualization".
}
\author{
Carolina Euan.}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
cutk
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
