\name{HMAlgo}
\alias{HMAlgo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hierarchical merger algorithm}
\description{
Compute the hierarchical merger clustering algorithm for a set of functions.
}
\usage{
HMAlgo(x,fx,parallel=FALSE,normalize=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{Argument of the functions.
}
\item{fx}{
Matrix of functions, the functions should be located by column.
}
}
\details{

}
\value{
An object with the following variables:

Diss.Matrix = Initial dissimilarity matrix.

min.value = trayectory of the minimum value.

Groups = list with the groupping structure at each step.
}
\references{
Euan, C., Ombao, H. and Ortega, J. (2016) "The Hierarchical Spectral Merger algorithm: A New Time Series Clustering Procedure".
}
\author{
Carolina Euan.}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
cutk
HSM
HCC
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
