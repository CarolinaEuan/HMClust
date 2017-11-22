\name{TVD}
\alias{TVD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Total variation distance
}
\description{
Computes the total variation distance between f1 and f2 with respect to the values w using the trapezoidal rule.
}
\usage{
TVD(w, f1, f2)
}
\arguments{
  \item{w}{
   Sorted vector of w values.
   %A vector with the grid of the range of values to compare f1 and f2.  
}
  \item{f1,f2}{Numeric vectors with the values of f1(w) and f2(w) which are going to be compare. 
  
  f1,f2 and w must have the same length. f1 and f2 must be normalized functions.
}
}
\details{
The TVD between two probability densities, f1 and f2, is defined as

d_{TV}(f1,f2) = 1 - integral( min(f1(w),f2(w)) dw). 
}
\value{
A numeric value between [0,1].
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Carolina Euan.
}

\examples{
##TVD between two normal densities 
w<-seq(0,5,length=1000)
f1<-dnorm(w,2,.5)
f2<-dnorm(w,2.5,.5)
diss<-TVD(w,f1,f2)
plot(w,f1,type="l",lwd=2,col=2,main=paste("TVD =",round(diss,3)),
      xlab="x",ylab="")
lines(w,f2,col=3,lwd=2)

##TVD between the normalized spectra of two AR2 processes
X1<-Sim.Ar(1000,12,1.01,100)
X2<-Sim.Ar(1000,15,1.01,100)
fest1<-spec.parzen(X1,a=300,dt=1/100)
fest2<-spec.parzen(X2,a=300,dt=1/100)
diss<-TVD(fest1[,1],fest1[,2]/var(X1),fest2[,2]/var(X2))
plot(fest1[,1],fest1[,2]/var(X1),type="l",lwd=2,col=2,main=paste("TVD =",round(diss,3)),xlab="w (Hz)",ylab="",ylim=c(0,max(fest1[,2]/var(X1),fest2[,2]/var(X2))))
lines(fest2[,1],fest2[,2]/var(X2),col=3,lwd=2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
