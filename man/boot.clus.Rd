\name{boot.clus}
\alias{boot.clus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bootstrap procedure
}
\description{
Executes a bootstrap test to decide the number of clusters when the HSM method is applied. Returns a p-value and if the test rejects K clusters at a level alpha.
}
\usage{
boot.clus(X,Clust,kg0,alpha=c(.01,.05,.1),nboot=1000,parallel=FALSE,freq=1,par.spectrum=c(100,1/(2*dt)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{X}{
Original time series by column.
}
  \item{Clust}{
Output from HSM.
}
  \item{kg0}{
Number of groups under the null hypothesis.
}
  \item{alpha}{
Level of the test.
}
  \item{nboot}{
Number of boostrap samples.
}
  \item{parallel}{
To run in a parallel loop. Library "doParallel" is requiered.
}
  \item{freq}{
Sampling frequency.
}
  \item{par.spectrum}{
Parameters for the estimation method: par.spectrum[1]=Bandwidth value, par.spectrum[2]= maximun evaluated frequency.
}
}
\author{
Carolina Euan.
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
HSM
}
