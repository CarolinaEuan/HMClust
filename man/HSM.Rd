\name{HSM}
\alias{HSM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hierarchical spectral merger algorithm}
\description{
Compute the hierarchical spectral merger clustering algorithm for a set of time series X.
}
\usage{
HSM(X,S=NULL,w=NULL, freq = 1, Merger = 1, par.spectrum = c(100, 1/(2 * dt), 512))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Matrix of time series, the series should be located by column.
}
\item{S}{
Optional. Matrix of spectrum as input.
}
  \item{w}{
Optional. Vector of Frequencies, required if matrix of spectrum is given as input.
}
  \item{freq}{
Sampling Frequency. Default value is 1.
}
  \item{Merger}{
If Merger==1 (default), the algorithm will estimate the new spectral density with the concatenated signals in order to get a better estimation of the original spectral density. If Merger==2 the algorithm will estimate the new spectral density with the mean spectrum using all time series in the cluster.
}
  \item{par.spectrum}{
Parameters for the estimation method: par.spectrum[1]=Bandwidth value, par.spectrum[2]= maximun evaluated frequency, par.spectrum[3]= length of the grid of the frequencies values.
}
}
\details{
Let X_1,X_2,...,X_N be the signal for each channel with length T and sampling frequency F_s. The procedure starts with N clusters, one for each individual channel.

1) Estimate the spectral density for each cluster using the smoothed periodogram and compute the TVD between their spectra.

2) Find the two clusters that have lower TVD, save this value as a characteristic.

3) Merge the signals in the two closest clusters and replace the two clusters by this new one.

4) Repeat steps 1-3 there is only one cluster left.
}
\value{
A HSM object with the following variables:

Diss.Matrix = Initial dissimilarity matrix.

min.value = trayectory of the minimum value.

Groups = list with the groupping structure at each step.
}
\references{
Euan, C., Ombao, H. and Ortega, J. (2016) "The Hierarchical Spectral Merger algo-
rithm: A New Time Series Clustering Procedure".
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
set.seed(9347)
#Series
TT<-1000
ww<-c(2,6,10,21,40);k<-6;ng<-3
# Original Groups (1,3,5) and (2,4,6)
e1<-c(0,1,0,1,0);e2<-c(0,1,1,0,0);A<-rbind(e1,e2)
ZZ<-matrix(NA,ncol=ng*5,nrow=TT)
for(index in 1:(ng*5))ZZ[,index]<-Sim.Ar(TT,ww[1+index\%\%5],1.01,100)
X<-cbind(t(A\%*\%t(ZZ[,1:5])),t(A\%*\%t(ZZ[,6:10])),t(A\%*\%t(ZZ[,11:15])))
#Algorithm
Fs<-100
Clust1<-HSM(X,freq=Fs)
cutk(Clust1,2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
