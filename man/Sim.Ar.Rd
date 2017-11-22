\name{Sim.Ar}
\alias{Sim.Ar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate from a reparametrized AR(2) model
}
\description{
Simulate an AR(2) process with parameters M, and eta. 
}
\usage{
Sim.Ar(Time, eta, M, Fs=1, m_burn = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Time}{
Length of output series. A strictly positive integer.
}
  \item{eta}{
Peak frequency of the spectrum. eta must be less or equal than Fs/2.
}
  \item{M}{
Modulus of the roots of phi(z), where $z_0^{1}=Conj(z_0^{2})$. M must be bigger than 1.
}
  \item{Fs}{
Sampling frequency (Hz). Default value is 1. 
}
  \item{m_burn}{
Length of ‘burn-in’ period.   
}
}
\details{
In time series methods there exist some parametric models such as the autoregressive process of order $p$ (AR($p$)) that have a closed form for the spectral density. In particular the AR(2) model can be parametrized as a function of the norm of the root of its characteristic polynomial (M) and the peak frequency of the spectrum (eta). These parameters will determine the shape of the spectral density (peak and sparseness). }
\value{
A vector of length Time from an AR(2) model. 
}
\examples{
#Sim.Ar
eta<-2;M<-1.01;Fs<-10
Xt<-Sim.Ar(1000,eta,M,Fs)
fest<-spec.parzen(Xt,a=100,dt=1/Fs)
plot(seq(.1,100,.1),Xt,type="l",main="AR(2) Process",xlab="Time (sec)",ylab="")
plot(fest,type="l",lwd=2,main="Smoothed Periodogram",xlab="w (Hz)",ylab="")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
