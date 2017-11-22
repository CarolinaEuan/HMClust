\name{cutk}
\alias{cutk}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
K groups from HSM, HCC or HM.
}
\description{
Returns k groups from a HSM, HCC or HM object.
}
\usage{
cutk(Clust, kg = NA, alpha = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Clust}{
Output from HSM, HCC or HM.
}
  \item{kg}{
Number of groups.
}
  \item{alpha}{
Dissimilarity value before the next clustering step.
}
}
\details{
Returns k groups from the result of the hierarchichal merger clustering method. The value k can be equal to kg or the last number of clusters before the next minimun value where bigger than alpha.
}
\value{
A list with k elements, where the jth element is a vector with the channels that belong to the jth cluster.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
HSM
HCC
HM
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
