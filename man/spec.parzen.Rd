\name{spec.parzen}
\alias{spec.parzen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Smoothed periodogram using a Parzen window 
}
\description{
One-side estimated spectrum using a lag window estimador with a parzen window.
}
\usage{
spec.parzen(x, a = 100, dt = 1, w0 = 10^(-5), wn = 1/(2 * dt), nn = 512)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Time series.
}
  \item{a}{
Bandwidth value.
}
  \item{dt}{
Sampling interval. Also, dt=1/Fs where Fs is the sampling frequency. Default value is 1.
}
  \item{w0,wn}{
Range of frequencies of interest. By default (10^{-5},Fs/2), where Fs is the sampling frequency.
}
  \item{nn}{
Number of evaluated frequencies in (w0,wn). 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A matrix of 2 columns and nn rows. Where the first column corresponds to the grid of frequencies and the second column corresponds to the spectrum at those frequencies. 
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Carolina Euan.}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#spec.parzen
x<-cos(2*pi*(10/1000)*(1:1000))+rnorm(1000)/2
fest<-spec.parzen(x,a=300,dt=.1)
par(mfrow=c(2,1))
plot(seq(.1,100,.1),x,type="l",main="Time series",xlab="Time (sec)",ylab="")
plot(fest,type="l",lwd=2,main="Smoothed Periodogram \n Bandwidth = 300",xlab="w (Hz)",ylab="")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
