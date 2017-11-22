\name{HMClust-package}
\alias{HMClust-package}
\alias{HMClust}
\docType{package}
\title{
Hierarchical Merger Clustering Methods
}
\description{
HMClust is a package developed for time series clustering or functional clustering based on hierarchical merger algorithms.}
\details{
\tabular{ll}{
Package: \tab HMClust\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2016-07-01\cr
License: \tab GLP-2\cr
}
This package includes the following functions: TVD, spec.parzen, Sim.Ar, HSM, HCC, C.Coh, HMAlgo, VisClust, cutk, boot.clus and check.clust.
}
\author{
Carolina Euan.
Maintainer: Carolina Euan <carolina.euan.campos@gmail.com>
}
\references{
Euan, C., Ombao, H. and Ortega, J. (2016) "The Hierarchical Spectral Merger algorithm: A New Time Series Clustering Procedure".

Euan, C., Sun, Y. and Ombao, H. (2017) "Coherence-based Time Series Clustering for Brain Connectivity Visualization".

}
\examples{
#TVD
#TVD between two normal densities
w<-seq(0,5,length=1000)
f1<-dnorm(w,2,.5)
f2<-dnorm(w,2.5,.5)
diss<-TVD(w,f1,f2)
plot(w,f1,type="l",lwd=2,col=2,main=paste("TVD =",round(diss,3)),
      xlab="x",ylab="")
lines(w,f2,col=3,lwd=2)

#spec.parzen
x<-cos(2*pi*(10/1000)*(1:1000))+rnorm(1000)/2
fest<-spec.parzen(x,a=300,dt=.1)
par(mfrow=c(2,1))
plot(seq(.1,100,.1),x,type="l",main="Time series",xlab="Time (sec)",ylab="")
plot(fest,type="l",lwd=2,main="Smoothed Periodogram \n Bandwidth = 300",xlab="w (Hz)",ylab="")

#Sim.Ar
eta<-2;M<-1.01;Fs<-10
Xt<-Sim.Ar(1000,eta,M,Fs)
fest<-spec.parzen(Xt,a=100,dt=1/Fs)
plot(seq(.1,100,.1),Xt,type="l",main="AR(2) Process",xlab="Time (sec)",ylab="")
plot(fest,type="l",lwd=2,main="Smoothed Periodogram",xlab="w (Hz)",ylab="")

#HSM
set.seed(9347)
#Series
TT<-1000
ww<-c(2,6,10,21,40) #Frequency involved in the basic signal
k<-6 #Number of Time series
ng<-3
# Original Groups
e1<-c(0,1,0,1,0);e2<-c(0,1,1,0,0);A<-rbind(e1,e2)
ZZ<-matrix(NA,ncol=ng*5,nrow=TT)
for(index in 1:(ng*5))ZZ[,index]<-Sim.Ar(TT,ww[1+index\%\%5],1.01,100)
X<-cbind(t(A\%*\%t(ZZ[,1:5])),t(A\%*\%t(ZZ[,6:10])),t(A\%*\%t(ZZ[,11:15])))
#Algorithm
Fs<-100
Clust1<-HSM(X,freq=Fs)
cutk(Clust1,2)}
