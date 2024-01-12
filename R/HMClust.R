 ## Distance functions
TVD<-function(w,f1,f2){
  #Compute TV distance using the trapezoidal rule
  if(length(w)!=length(f1)|length(w)!=length(f2))stop("w, f1 and f2 must have the same length")
  n<-length(w)
  int<-(w[2:n]-w[1:(n-1)])%*%(pmin(f1,f2)[2:n]+pmin(f1,f2)[1:(n-1)])/2
  return(as.double(1-int))
}

C.Coh<-function(A,IndCaseL,p=1){
  #Cluster Coherence function
  if(p==1)eigennorm<-function(x) sum(abs(x)) else eigennorm<-function(x) (sum((abs(x))^p))^(1/p)
  eigenA0<-sort(abs(eigen(A)$values),decreasing = TRUE)
  IndCase<-c(IndCaseL[[1]],IndCaseL[[2]])
  eigenA<-eigenA0/eigennorm(eigenA0)
  eigenI<-sort(IndCase,decreasing = TRUE)/eigennorm(IndCase)
  dissaux<-eigennorm(eigenA-eigenI)
  return(list(dissaux,eigenA0))
}

##Spectral estimators
spec.parzen<-function(x,a=100,dt=1,w0=10^(-5),wn=1/(2*dt),nn=512)
{
  # Compute the smoothed periodogram using a Parzen window
  if(a>=(length(x)-1))stop("The bandwidth value is too big")
  kp<-numeric(2*a-1)
  tt<-seq(-1,1,length=2*a-1)
  kp[abs(tt)<.5]<-1-6*abs(tt[abs(tt)<.5])^2+6*(abs(tt[abs(tt)<.5])^3)
  kp[abs(tt)>=.5]<-2*(1-abs(tt[abs(tt)>=.5]))^3
  Cov<-as.vector(acf(x, lag.max=a-1,type = "covariance",plot=F,na.action = na.exclude)$acf)
  CovS<-kp[a:(2*a-1)]*Cov
  time<-seq(0,length(x)*dt,by=dt)[1:a]
  S<-function(ww)return((2*dt/pi)*sum(cos(2*pi*ww*time)[-1]*CovS[-1])+(dt/pi)*CovS[1])
  SS<-c()
  w<-seq(w0,wn,length.out=nn)
  for(i in 1:nn)SS[i]<-S(w[i])
  Result<-matrix(c(w,SS),ncol=2)
  return(Result)
}

Spectral.matrix<-function(X,Fs,op=1,bw=20)
{
  #Computes the spectral matrix using fft and a smoothing is done by the ferej kernel
  Xfft<-mvfft(X)
  N<-nrow(X)
  if(op==1){
    kernel.option<-kernel("fejer",bw,r=2)
    #plot(kernel.option)
    SS <- array(NA, dim = c(N, ncol(X), ncol(X)))
    for (i in 1L:ncol(X)) {
      for (j in 1L:ncol(X)) {
        SS[, i, j] <- Xfft[, i] * Conj(Xfft[, j])/(N*Fs)
        SS[1, i, j] <- 0.5 * (SS[2, i, j] + SS[N, i, j])
        SS[, i, j]<-kernapply(SS[, i, j], kernel.option, circular = TRUE)
      }
    }
  }
  else{
    SS <- array(NA, dim = c(N, ncol(X), ncol(X)))
    for (i in 1L:ncol(X)) {
      for (j in 1L:ncol(X)) {
        SS[, i, j] <- Xfft[, i] * Conj(Xfft[, j])/(N*Fs)
        SS[1, i, j] <- 0.5 * (SS[2, i, j] + SS[N, i, j])
      }
    }
  }
  SS<-SS[1:(N/2),,]
  w<-seq(0,(Fs/2),length=N/2)
  Spec.Matrix<-list(freq=w,f=SS)
  return(Spec.Matrix)
}

##Clustering Algorithms

HSM<-function(X,S=NULL,w=NULL,freq=1,Merger=1,par.spectrum=c(100,1/(2*dt),512),parallel=FALSE){
  ##Compute the HSM algorithm using the TVD
  dt<-1/freq
  if(is.null(S)){
    #Time series as input
    #Initial parameters
    n<-length(X[,1])
    k<-length(X[1,])
    #Time<-n*dt
    np<-length(par.spectrum)
    #if(np==0){a<-100; wn<-1/(2*dt);length.w<-512}
    if(np==1){a<-par.spectrum[1]; wn<-1/(2*dt); length.w<-512}
    if(np==2){a<-par.spectrum[1]; wn<-par.spectrum[2]; length.w<-512}
    if(np==3){a<-par.spectrum[1]; wn<-par.spectrum[2]; length.w<-par.spectrum[3]}
    w<-seq(10^(-5),wn,length.out=length.w)
    #Initial normalized estimated Spectra and Dissimilarity Matrix
    if(parallel){
      library("doParallel")
      cl<-min(k,detectCores()-1)
      registerDoParallel(cl)
      Spec<-foreach(i=1:(k-1),.combine = cbind) %dopar%
      {
        Aux_Spec<-spec.parzen(X[,i],a=a,dt=dt,wn=wn,nn=length.w)
        Aux_Spec[,2]/((w[2:length.w]-w[1:(length.w-1)])%*%(Aux_Spec[2:length.w,2]+Aux_Spec[1:(length.w-1),2])/2)
      }
      Aux_Spec<-spec.parzen(X[,k],a=a,dt=dt,wn=wn,nn=length.w)
      Spec<-cbind(Spec,Aux_Spec[,2]/((w[2:length.w]-w[1:(length.w-1)])%*%(Aux_Spec[2:length.w,2]+Aux_Spec[1:(length.w-1),2])/2))

    }else{
      Spec<-matrix(NA,nrow=length.w,ncol=k)
      for(i in 1:k){
        Aux_Spec<-spec.parzen(X[,i],a=a,dt=dt,wn=wn,nn=length.w)
        w<-Aux_Spec[,1]
        Spec[,i]<-Aux_Spec[,2]/((w[2:length.w]-w[1:(length.w-1)])%*%(Aux_Spec[2:length.w,2]+Aux_Spec[1:(length.w-1),2])/2)
      }}
  }else{
    #Spectrum as input
    Spec<-S
    k<-ncol(S)
    if(is.null(w)){
      np<-length(par.spectrum)
      if(np==1){a<-par.spectrum[1]; wn<-1/(2*dt); length.w<-512}
      if(np==2){a<-par.spectrum[1]; wn<-par.spectrum[2]; length.w<-512}
      if(np==3){a<-par.spectrum[1]; wn<-par.spectrum[2]; length.w<-par.spectrum[3]}
      w<-seq(10^(-5),wn,length.out=nn)
    }
    length.w<-length(w)
    for(i in 1:k){
      Aux_Spec<-Spec[,i]
      Spec[,i]<-Aux_Spec/((w[2:length.w]-w[1:(length.w-1)])%*%(Aux_Spec[2:length.w]+Aux_Spec[1:(length.w-1)])/2)
    }
    if(Merger==1){
      print("Merger has been assigned to 2")
      Merger<-2
    }
  }
  if(parallel){
    library("doParallel")
    cl<-min(k,detectCores()-1)
    registerDoParallel(cl)
    MatDiss<-foreach(i=1:k,.combine = rbind) %dopar%
    {
      Aux.MatDiss<-rep(NA,k)
      for(j in i:k) Aux.MatDiss[j]<-TVD(w,Spec[,i],Spec[,j])
      Aux.MatDiss
    }
  }else{
    MatDiss<-matrix(NA,k,k)
    for(i in 1:k)for(j in i:k) MatDiss[i,j]<-TVD(w,Spec[,i],Spec[,j])
  }
  diag(MatDiss)<-0
  MatDiss[lower.tri(MatDiss)]<-t(MatDiss)[lower.tri(MatDiss)]
  #Dinamic Diss Mat
  MD_Change<-MatDiss# initial MD
  diag(MD_Change)<-rep(Inf,k)
  #Dinamic estimated spectra
  Spec_Change<-Spec
  # initial groups
  g<-as.list(1:k)
  #Save variables
  min.value<-numeric(k-1)
  groups<-list()
  #Version 1
  if(Merger==1){
    Xlist<-list()
    for(j in 1:k)Xlist[[j]]<-X[,j]#initial series
    for(ite in 1:(k-1)){
      #############
      # Identify the closest clusters and the minimun value
      min.value[ite]<-min(MD_Change) ##
      aux1<-which(MD_Change==min(MD_Change),arr.ind = TRUE)[1,]
      g<-c(g[-aux1],list(unlist(g[aux1])))# new groups
      groups[[ite]]<-g ##
      #############
      # Merge Time Series
      Xlist<-c(Xlist[-aux1],list(unlist(Xlist[sort(aux1)])))
      Spec_Change<-Spec_Change[,-aux1]
      Spec_Xnew<-spec.parzen(Xlist[[k-ite]],a=a,dt=dt,wn=wn,nn=length.w)
      Spec_Change<-cbind(Spec_Change,Spec_Xnew[,2]/((w[2:length.w]-w[1:(length.w-1)])%*%(Spec_Xnew[2:length.w,2]+Spec_Xnew[1:(length.w-1),2])/2))
      #############
      #Compute new TVD
      MD_Change<-MD_Change[-aux1,-aux1]
      if(ite<(k-1)){
        aux2<-numeric(k-ite-1)
        for(i in 1:(k-ite-1)){
          aux2[i]<-TVD(w,Spec_Change[,i],Spec_Change[,k-ite])
        }
      MD_Change<-rbind(MD_Change,aux2)
      MD_Change<-cbind(MD_Change,c(aux2,Inf))# new MD
      }else{MD_Change<-0}
    }
  }
  #Version 2
  if(Merger==2){
    for(ite in 1:(k-1)){
      #############
      # Identify the closest clusters and the minimun value
      min.value[ite]<-min(MD_Change)
      aux1<-which(MD_Change==min(MD_Change),arr.ind = TRUE)[1,]
      g<-c(g[-aux1],list(unlist(g[aux1])))# new groups
      groups[[ite]]<-g ##
      #############
      # Spectral Merge
      Spec_Change<-Spec_Change[,-aux1]
      Spec_Change<-cbind(Spec_Change,rowMeans(Spec[,g[[k-ite]]]))
      #############
      #Compute new TVD
      MD_Change<-MD_Change[-aux1,-aux1]
      if(ite<(k-1)){
        aux2<-numeric(k-ite-1)
        for(i in 1:(k-ite-1)){
          aux2[i]<-TVD(w,Spec_Change[,i],Spec_Change[,k-ite])
        }
      MD_Change<-rbind(MD_Change,aux2)
      MD_Change<-cbind(MD_Change,c(aux2,Inf))# new MD
      }else{MD_Change<-0}

    }
  }
  out <- list(MatDiss,min.value, groups)
  names(out) <- c("Diss.Matrix","min.value", "Groups")
  out
}

HCC<-function(X,Clustfreq=NULL,freq=1,dist=1){
  #X must be scaled (mean=0)
  ##Compute the Clustering Akgorithm using coherence matrix in a particular fequency
  # Remark if w is not a fourier frequency then we will compute the clustering at the biggest value <=
  #Initial parameters
  n<-length(X[,1])
  k<-length(X[1,])
  #Initial estimated Spectra and Dissimilarity Matrix
  SSlist<-Spectral.matrix(X,freq)
  w<-SSlist$freq;SS<-SSlist$f
  Coh.Mat<- array(NA,dim=c(n/2,k,k))
  for(i in 1:k)for(j in i:k) {
    Coh.Mat[,i,j]<-Re((Re(SS[,i,j])^2+Im(SS[,i,j])^2)/((SS[,i,i])*(SS[,j,j])))
  }
  out<-list(w,SS,Coh.Mat) #If user does not choose a frequency band then it returns SS and coherence matrix
  basic.names<-c("w","Spectral.Matrix","Coherence")
  names(out) <- basic.names
  if(is.null(Clustfreq)==FALSE){
    ind.band<-findInterval(Clustfreq,w)#
    MatDiss<- Coh.Mat[ind.band[1]:ind.band[length(ind.band)],,]
    #Integrate Coherence in case of band freq
    if(length(ind.band)>1){
      MatDissAux<-matrix(0,nrow=k,ncol=k)
      for(bandcoh in 1:length(ind.band))MatDissAux<-MatDissAux+as.matrix(MatDiss[bandcoh,,])
      MatDiss<-MatDissAux/(length(ind.band)) #(wn-w0)/n partiton size / (wn-w0) normalization constant
    }
    MatDiss[lower.tri(MatDiss)]<-t(MatDiss)[lower.tri(MatDiss)]
    # Dinamic Diss Mat
    MD_Change<-1-MatDiss# initial MD #
    diag(MD_Change)<-rep(Inf,k)
    # initial groups
    g<-as.list(1:k)
    eigen.values<-as.list(rep(1,k))
    #Save variables
    min.value<-numeric(k-1)
    groups<-list()
    groups.eigen.values<-list()
    #Clustering using coherence
    for(ite in 1:(k-1)){
      #############
      # Identify the closest clusters and the maximun value (highest coherence)
      min.value[ite]<-min(MD_Change) ##
      aux1<-which(MD_Change==min(MD_Change),arr.ind = TRUE)[1,]
      eigen.values<-c(eigen.values[-aux1],list(eigen(MatDiss[unlist(g[aux1]),unlist(g[aux1])])$values)) #new eigen values
      g<-c(g[-aux1],list(unlist(g[aux1])))# new groups
      groups[[ite]]<-g ##
      groups.eigen.values[[ite]]<-eigen.values ##
      #############
      #Compute new dissimilarity using eigenvalues
      MD_Change<-MD_Change[-aux1,-aux1]
      aux2<-numeric(k-ite-1)
      if(ite<(k-1)){
        for(i in 1:(k-ite-1)){
          set1<-unlist(g[i])
          set2<-unlist(g[k-ite])
          Diss_aux<-C.Coh(MatDiss[c(set1,set2),c(set1,set2)], #
                          IndCase=list(unlist(eigen.values[i]),unlist(eigen.values[k-ite])),p=dist)#
          aux2[i]<-1-Diss_aux[[1]]
        }
        MD_Change<-rbind(MD_Change,aux2)
        MD_Change<-cbind(MD_Change,c(aux2,Inf))# new MD
      }
    }
    out <- list(1-MatDiss, min.value, groups, groups.eigen.values)
    names(out) <- c("Diss.Matrix","min.value", "Groups","Lambda.Groups")
  }
  out
}

HMAlgo<-function(x,fx,parallel=FALSE,normalize=TRUE,TVD=TRUE){
  n<-length(fx[,1])
  k<-length(fx[1,])
  if(normalize)fxN<-apply(fx,2,normalize,x=x) else fxN<-fx
  if(parallel){
    library("doParallel")
    cl<-min(k,detectCores()-1)
    registerDoParallel(cl)
    if(TVD){
    MatDiss<-foreach(i=1:k,.combine = rbind) %dopar%
    {
      Aux.MatDiss<-rep(NA,k)
      for(j in i:k) Aux.MatDiss[j]<-TVD(x,fxN[,i],fxN[,j])
      Aux.MatDiss
    }
    }else{
      MatDiss<-foreach(i=1:k,.combine = rbind) %dopar%
      {
        Aux.MatDiss<-rep(NA,k)
        for(j in i:k) Aux.MatDiss[j]<-sqrt(sum((fxN[,i]-fxN[,j])^2))
        Aux.MatDiss
      }
    }
  }else{
    MatDiss<-matrix(NA,k,k)
    if(TVD){
      for(i in 1:k)for(j in i:k) MatDiss[i,j]<-TVD(x,fxN[,i],fxN[,j])
    }else{
      for(i in 1:k)for(j in i:k) MatDiss[i,j]<-sqrt(sum((fxN[,i]-fxN[,j])^2))
    }
  }
  diag(MatDiss)<-0
  MatDiss[lower.tri(MatDiss)]<-t(MatDiss)[lower.tri(MatDiss)]
  #Dinamic Diss Mat
  MD_Change<-MatDiss# initial MD
  diag(MD_Change)<-rep(Inf,k)
  #Dinamic estimated spectra
  fx_Change<-fxN
  # initial groups
  g<-as.list(1:k)
  #Save variables
  min.value<-numeric(k-1)
  groups<-list()
  #Average version
  #Version 2
  if(TVD){
    for(ite in 1:(k-1)){
      #############
      # Identify the closest clusters and the minimun value
      min.value[ite]<-min(MD_Change)
      aux1<-which(MD_Change==min(MD_Change),arr.ind = TRUE)[1,]
      g<-c(g[-aux1],list(unlist(g[aux1])))# new groups
      groups[[ite]]<-g ##
      #############
      # Spectral Merge
      fx_Change<-fx_Change[,-aux1]
      fx_Change<-cbind(fx_Change,rowMeans(fxN[,g[[k-ite]]]))
      #############
      #Compute new TVD
      MD_Change<-MD_Change[-aux1,-aux1]
      if(ite<(k-1)){
        aux2<-numeric(k-ite-1)
        for(i in 1:(k-ite-1)){
          aux2[i]<-TVD(x,fx_Change[,i],fx_Change[,k-ite])
        }
        MD_Change<-rbind(MD_Change,aux2)
        MD_Change<-cbind(MD_Change,c(aux2,Inf))# new MD
      }else{MD_Change<-0}
    }
  }else{
    for(ite in 1:(k-1)){
      #############
      # Identify the closest clusters and the minimun value
      min.value[ite]<-min(MD_Change)
      aux1<-which(MD_Change==min(MD_Change),arr.ind = TRUE)[1,]
      g<-c(g[-aux1],list(unlist(g[aux1])))# new groups
      groups[[ite]]<-g ##
      #############
      # Spectral Merge
      fx_Change<-fx_Change[,-aux1]
      fx_Change<-cbind(fx_Change,rowMeans(fxN[,g[[k-ite]]]))
      #############
      #Compute new TVD
      MD_Change<-MD_Change[-aux1,-aux1]
      if(ite<(k-1)){
        aux2<-numeric(k-ite-1)
        for(i in 1:(k-ite-1)){
          aux2[i]<-sqrt(sum((fx_Change[,i]-fx_Change[,k-ite])^2))
        }
        MD_Change<-rbind(MD_Change,aux2)
        MD_Change<-cbind(MD_Change,c(aux2,Inf))# new MD
      }else{MD_Change<-0}
    }
  }
  out <- list(MatDiss,min.value, groups)
  names(out) <- c("Diss.Matrix","min.value", "Groups")
  out
}

HM2D<-function (x,y, fxy, parallel = FALSE, normalize = TRUE){
  k <- length(fxy)
  if (normalize)
    fxyN <- lapply(fxy,normalize, x = x, y = y, D=2)
  else fxyN <- fxy
  if (parallel) {
    library("doParallel")
    cl <- min(k, detectCores() - 1)
    registerDoParallel(cl)
    MatDiss <- foreach(i = 1:k, .combine = rbind) %dopar%
    {
      Aux.MatDiss <- rep(NA, k)
      for (j in i:k) Aux.MatDiss[j] <- TVD2(x, y, fxyN[[i]],fxyN[[j]])
      Aux.MatDiss
    }
  }
  else {
    MatDiss <- matrix(NA, k, k)
    for (i in 1:k) for (j in i:k) MatDiss[i, j] <- TVD2(x, y, fxyN[[i]],fxyN[[j]])
  }
  diag(MatDiss) <- 0
  MatDiss[lower.tri(MatDiss)] <- t(MatDiss)[lower.tri(MatDiss)]
  MD_Change <- MatDiss
  diag(MD_Change) <- rep(Inf, k)
  fxy_Change <- fxyN
  g <- as.list(1:k)
  min.value <- numeric(k - 1)
  groups <- list()
  for (ite in 1:(k - 1)) {
    min.value[ite] <- min(MD_Change)
    aux1 <- which(MD_Change == min(MD_Change), arr.ind = TRUE)[1,                                              ]
    g <- c(g[-aux1], list(unlist(g[aux1])))
    groups[[ite]] <- g
    fxy_Change <- fxy_Change[-aux1]
    fxy_Change[[length(fxy_Change)+1]] <- Reduce("+", fxyN[g[[k - ite]]]) / length(g[[k - ite]])
    MD_Change <- MD_Change[-aux1, -aux1]
    if (ite < (k - 1)) {
      aux2 <- numeric(k - ite - 1)
      for (i in 1:(k - ite - 1)) {
        aux2[i] <- TVD2(x, y, fxy_Change[[i]],fxy_Change[[k-ite]])
      }
      MD_Change <- rbind(MD_Change, aux2)
      MD_Change <- cbind(MD_Change, c(aux2, Inf))
    }
    else {
      MD_Change <- 0
    }
  }
  out <- list(MatDiss, min.value, groups)
  names(out) <- c("Diss.Matrix", "min.value", "Groups")
  out
}

##Auxiliar functions
Sim.Ar<-function(Time,eta,M,Fs=1,m_burn=100)
{
  #Draw from an AR(2) process with parametrization (\eta,M)
  if(eta>Fs/2)stop("Eta can not be bigger that Fs/2")
  if(M<=1)warning("The AR(2) Process is not causal")
  w_0<-2*pi*eta/Fs
  phi1<-(2*cos(w_0))/M
  phi2<-(-1)/(M^2)
  return(as.matrix(arima.sim(list(order = c(2,0,0), ar = c(phi1,phi2)), n = Time,n.star=m_burn)))
}

normalize<-function(x,fx,normfx=FALSE,D=1,y=NULL)
{
  #normalization with trapezoidal rule
  if(D==1){
    length.x<-length(x)
    inte<-((x[2:length.x] - x[1:(length.x-1)]) %*% (fx[2:length.x] + fx[1:(length.x-1)])/2)
  }
  if(D==2){
    length.x<-length(x)
    length.y<-length(y)
    We<-matrix(4,nrow=length.x,ncol=length.y)
    We[c(1,1,length.x,length.x),c(1,length.y,1,length.y)]<-1
    We[2:(length.x-1),1]<-2
    We[2:(length.x-1),length.y]<-2
    We[1,2:(length.y-1)]<-2
    We[length.x,2:(length.y-1)]<-2
    inte<-sum(fx*We)*(x[2]-x[1])*(y[2]-y[1])*(1/4)
  }
  if(normfx){
    return(inte)
  }else{
    if(inte==0){
      return(fx)
    }else{
      return(fx/inte)}
  }
}

TVD2<-function(x,y,fxy1,fxy2){
  length.x<-length(x)
  length.y<-length(y)
  We<-matrix(4,nrow=length.x,ncol=length.y)
  We[c(1,1,length.x,length.x),c(1,length.y,1,length.y)]<-1
  We[2:(length.x-1),1]<-2
  We[2:(length.x-1),length.y]<-2
  We[1,2:(length.y-1)]<-2
  We[length.x,2:(length.y-1)]<-2
  return((1/2)*sum(We*abs(fxy1-fxy2))*(1/4)*(x[2]-x[1])*(y[2]-y[1]))
}

##Clustering Visualization
#Need the ColPal
VisClust<-function(Clust,Order=NULL,kg=NULL,max.kg=min(200,length(Clust$min.value)),cex.xaxis=1.2,namesX=NULL,nplot=NULL){
  opar <- par()      # make a copy of current settings

  colPal<-c("#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
            "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
            "#5A0007", "#809693", "#FEFFE6", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
            "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
            "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
            "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
            "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
            "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
            "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
            "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
            "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
            "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
            "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C",
            "#83AB58", "#001C1E", "#D1F7CE", "#004B28", "#C8D0F6", "#A3A489", "#806C66", "#222800",
            "#BF5650", "#E83000", "#66796D", "#DA007C", "#FF1A59", "#8ADBB4", "#1E0200", "#5B4E51",
            "#C895C5", "#320033", "#FF6832", "#66E1D3", "#CFCDAC", "#D0AC94", "#7ED379", "#012C58",
            "#7A7BFF", "#D68E01", "#353339", "#78AFA1", "#FEB2C6", "#75797C", "#837393", "#943A4D",
            "#B5F4FF", "#D2DCD5", "#9556BD", "#6A714A", "#001325", "#02525F", "#0AA3F7", "#E98176",
            "#DBD5DD", "#5EBCD1", "#3D4F44", "#7E6405", "#02684E", "#962B75", "#8D8546", "#9695C5",
            "#E773CE", "#D86A78", "#3E89BE", "#CA834E", "#518A87", "#5B113C", "#55813B", "#E704C4",
            "#00005F", "#A97399", "#4B8160", "#59738A", "#FF5DA7", "#F7C9BF", "#643127", "#513A01",
            "#6B94AA", "#51A058", "#A45B02", "#1D1702", "#E20027", "#E7AB63", "#4C6001", "#9C6966",
            "#64547B", "#97979E", "#006A66", "#391406", "#F4D749", "#0045D2", "#006C31", "#DDB6D0",
            "#7C6571", "#9FB2A4", "#00D891", "#15A08A", "#BC65E9", "#FFFFFE", "#C6DC99", "#203B3C",
            "#671190", "#6B3A64", "#F5E1FF", "#FFA0F2", "#CCAA35", "#374527", "#8BB400", "#797868",
            "#C6005A", "#3B000A", "#C86240", "#29607C", "#402334", "#7D5A44", "#CCB87C", "#B88183",
            "#AA5199", "#B5D6C3", "#A38469", "#9F94F0", "#A74571", "#B894A6", "#71BB8C", "#00B433",
            "#789EC9", "#6D80BA", "#953F00", "#5EFF03", "#E4FFFC", "#1BE177", "#BCB1E5", "#76912F",
            "#003109", "#0060CD", "#D20096", "#895563", "#29201D", "#5B3213", "#A76F42", "#89412E",
            "#1A3A2A", "#494B5A", "#A88C85", "#F4ABAA", "#A3F3AB", "#00C6C8", "#EA8B66", "#958A9F",
            "#BDC9D2", "#9FA064", "#BE4700", "#658188", "#83A485", "#453C23", "#47675D", "#3A3F00",
            "#061203", "#DFFB71", "#868E7E", "#98D058", "#6C8F7D", "#D7BFC2", "#3C3E6E", "#D83D66",
            "#2F5D9B", "#6C5E46", "#D25B88", "#5B656C", "#00B57F", "#545C46", "#866097", "#365D25",
            "#252F99", "#00CCFF", "#674E60", "#FC009C", "#92896B", "#1E2324", "#DEC9B2", "#9D4948",
            "#85ABB4", "#342142", "#D09685", "#A4ACAC", "#00FFFF", "#AE9C86", "#742A33", "#0E72C5",
            "#AFD8EC", "#C064B9", "#91028C", "#FEEDBF", "#FFB789", "#9CB8E4", "#AFFFD1", "#2A364C",
            "#4F4A43", "#647095", "#34BBFF", "#807781", "#920003", "#B3A5A7", "#018615", "#F1FFC8",
            "#976F5C", "#FF3BC1", "#FF5F6B", "#077D84", "#F56D93", "#5771DA", "#4E1E2A", "#830055",
            "#02D346", "#BE452D", "#00905E", "#BE0028", "#6E96E3", "#007699", "#FEC96D", "#9C6A7D",
            "#3FA1B8", "#893DE3", "#79B4D6", "#7FD4D9", "#6751BB", "#B28D2D", "#E27A05", "#DD9CB8",
            "#AABC7A", "#980034", "#561A02", "#8F7F00", "#635000", "#CD7DAE", "#8A5E2D", "#FFB3E1",
            "#6B6466", "#C6D300", "#0100E2", "#88EC69", "#8FCCBE", "#21001C", "#511F4D", "#E3F6E3",
            "#FF8EB1", "#6B4F29", "#A37F46", "#6A5950", "#1F2A1A", "#04784D", "#101835", "#E6E0D0",
            "#FF74FE", "#00A45F", "#8F5DF8", "#4B0059", "#412F23", "#D8939E", "#DB9D72", "#604143",
            "#B5BACE", "#989EB7", "#D2C4DB", "#A587AF", "#77D796", "#7F8C94", "#FF9B03", "#555196",
            "#31DDAE", "#74B671", "#802647", "#2A373F", "#014A68", "#696628", "#4C7B6D", "#002C27",
            "#7A4522", "#3B5859", "#E5D381", "#FFF3FF", "#679FA0", "#261300", "#2C5742", "#9131AF",
            "#AF5D88", "#C7706A", "#61AB1F", "#8CF2D4", "#C5D9B8", "#9FFFFB", "#BF45CC", "#493941",
            "#863B60", "#B90076", "#003177", "#C582D2", "#C1B394", "#602B70", "#887868", "#BABFB0",
            "#030012", "#D1ACFE", "#7FDEFE", "#4B5C71", "#A3A097", "#E66D53", "#637B5D", "#92BEA5",
            "#00F8B3", "#BEDDFF", "#3DB5A7", "#DD3248", "#B6E4DE", "#427745", "#598C5A", "#B94C59",
            "#8181D5", "#94888B", "#FED6BD", "#536D31", "#6EFF92", "#E4E8FF", "#20E200", "#FFD0F2",
            "#4C83A1", "#BD7322", "#915C4E", "#8C4787", "#025117", "#A2AA45", "#2D1B21", "#A9DDB0",
            "#FF4F78", "#528500", "#009A2E", "#17FCE4", "#71555A", "#525D82", "#00195A", "#967874",
            "#555558", "#0B212C", "#1E202B", "#EFBFC4", "#6F9755", "#6F7586", "#501D1D", "#372D00",
            "#741D16", "#5EB393", "#B5B400", "#DD4A38", "#363DFF", "#AD6552", "#6635AF", "#836BBA",
            "#98AA7F", "#464836", "#322C3E", "#7CB9BA", "#5B6965", "#707D3D", "#7A001D", "#6E4636",
            "#443A38", "#AE81FF", "#489079", "#897334", "#009087", "#DA713C", "#361618", "#FF6F01",
            "#006679", "#370E77", "#4B3A83", "#C9E2E6", "#C44170", "#FF4526", "#73BE54", "#C4DF72",
            "#ADFF60", "#00447D", "#DCCEC9", "#BD9479", "#656E5B", "#EC5200", "#FF6EC2", "#7A617E",
            "#DDAEA2", "#77837F", "#A53327", "#608EFF", "#B599D7", "#A50149", "#4E0025", "#C9B1A9",
            "#03919A", "#1B2A25", "#E500F1", "#982E0B", "#B67180", "#E05859", "#006039", "#578F9B",
            "#305230", "#CE934C", "#B3C2BE", "#C0BAC0", "#B506D3", "#170C10", "#4C534F", "#224451",
            "#3E4141", "#78726D", "#B6602B", "#200441", "#DDB588", "#497200", "#C5AAB6", "#033C61",
            "#71B2F5", "#A9E088", "#4979B0", "#A2C3DF", "#784149", "#2D2B17", "#3E0E2F", "#57344C",
            "#0091BE", "#E451D1", "#4B4B6A", "#5C011A", "#7C8060", "#FF9491", "#4C325D", "#005C8B",
            "#E5FDA4", "#68D1B6", "#032641", "#140023", "#8683A9", "#CFFF00", "#A72C3E", "#34475A",
            "#B1BB9A", "#B4A04F", "#8D918E", "#A168A6", "#813D3A", "#425218", "#DA8386", "#776133",
            "#563930", "#8498AE", "#90C1D3", "#B5666B", "#9B585E", "#856465", "#AD7C90", "#E2BC00",
            "#E3AAE0", "#B2C2FE", "#FD0039", "#009B75", "#FFF46D", "#E87EAC", "#DFE3E6", "#848590",
            "#AA9297", "#83A193", "#577977", "#3E7158", "#C64289", "#EA0072", "#C4A8CB", "#55C899",
            "#E78FCF", "#004547", "#F6E2E3", "#966716", "#378FDB", "#435E6A", "#DA0004", "#1B000F",
            "#5B9C8F", "#6E2B52", "#011115", "#E3E8C4", "#AE3B85", "#EA1CA9", "#FF9E6B", "#457D8B",
            "#92678B", "#00CDBB", "#9CCC04", "#002E38", "#96C57F", "#CFF6B4", "#492818", "#766E52",
            "#20370E", "#E3D19F", "#2E3C30", "#B2EACE", "#F3BDA4", "#A24E3D", "#976FD9", "#8C9FA8",
            "#7C2B73", "#4E5F37", "#5D5462", "#90956F", "#6AA776", "#DBCBF6", "#DA71FF", "#987C95",
            "#52323C", "#BB3C42", "#584D39", "#4FC15F", "#A2B9C1", "#79DB21", "#1D5958", "#BD744E",
            "#160B00", "#20221A", "#6B8295", "#00E0E4", "#102401", "#1B782A", "#DAA9B5", "#B0415D",
            "#859253", "#97A094", "#06E3C4", "#47688C", "#7C6755", "#075C00", "#7560D5", "#7D9F00",
            "#C36D96", "#4D913E", "#5F4276", "#FCE4C8", "#303052", "#4F381B", "#E5A532", "#706690",
            "#AA9A92", "#237363", "#73013E", "#FF9079", "#A79A74", "#029BDB", "#FF0169", "#C7D2E7",
            "#CA8869", "#80FFCD", "#BB1F69", "#90B0AB", "#7D74A9", "#FCC7DB", "#99375B", "#00AB4D",
            "#ABAED1", "#BE9D91", "#E6E5A7", "#332C22", "#DD587B", "#F5FFF7", "#5D3033", "#6D3800",
            "#FF0020", "#B57BB3", "#D7FFE6", "#C535A9", "#260009", "#6A8781", "#A8ABB4", "#D45262",
            "#794B61", "#4621B2", "#8DA4DB", "#C7C890", "#6FE9AD", "#A243A7", "#B2B081", "#181B00",
            "#286154", "#4CA43B", "#6A9573", "#A8441D", "#5C727B", "#738671", "#D0CFCB", "#897B77",
            "#1F3F22", "#4145A7", "#DA9894", "#A1757A", "#63243C", "#ADAAFF", "#00CDE2", "#DDBC62",
            "#698EB1", "#208462", "#00B7E0", "#614A44", "#9BBB57", "#7A5C54", "#857A50", "#766B7E",
            "#014833", "#FF8347", "#7A8EBA", "#274740", "#946444", "#EBD8E6", "#646241", "#373917",
            "#6AD450", "#81817B", "#D499E3", "#979440", "#011A12", "#526554", "#B5885C", "#A499A5",
            "#03AD89", "#B3008B", "#E3C4B5", "#96531F", "#867175", "#74569E", "#617D9F", "#E70452",
            "#067EAF", "#A697B6", "#B787A8", "#9CFF93", "#311D19", "#3A9459", "#6E746E", "#B0C5AE",
            "#84EDF7", "#ED3488", "#754C78", "#384644", "#C7847B", "#00B6C5", "#7FA670", "#C1AF9E",
            "#2A7FFF", "#72A58C", "#FFC07F", "#9DEBDD", "#D97C8E", "#7E7C93", "#62E674", "#B5639E",
            "#FFA861", "#C2A580", "#8D9C83", "#B70546", "#372B2E", "#0098FF", "#985975", "#20204C",
            "#FF6C60", "#445083", "#8502AA", "#72361F", "#9676A3", "#484449", "#CED6C2", "#3B164A",
            "#CCA763", "#2C7F77", "#02227B", "#A37E6F", "#CDE6DC", "#CDFFFB", "#BE811A", "#F77183",
            "#EDE6E2", "#CDC6B4", "#FFE09E", "#3A7271", "#FF7B59", "#4E4E01", "#4AC684", "#8BC891",
            "#BC8A96", "#CF6353", "#DCDE5C", "#5EAADD", "#F6A0AD", "#E269AA", "#A3DAE4", "#436E83",
            "#002E17", "#ECFBFF", "#A1C2B6", "#50003F", "#71695B", "#67C4BB", "#536EFF", "#5D5A48",
            "#890039", "#969381", "#371521", "#5E4665", "#AA62C3", "#8D6F81", "#2C6135", "#410601",
            "#564620", "#E69034", "#6DA6BD", "#E58E56", "#E3A68B", "#48B176", "#D27D67", "#B5B268",
            "#7F8427", "#FF84E6", "#435740", "#EAE408", "#F4F5FF", "#325800", "#4B6BA5", "#ADCEFF",
            "#9B8ACC", "#885138", "#5875C1", "#7E7311", "#FEA5CA", "#9F8B5B", "#A55B54", "#89006A",
            "#AF756F", "#2A2000", "#7499A1", "#FFB550", "#00011E", "#D1511C", "#688151", "#BC908A",
            "#78C8EB", "#8502FF", "#483D30", "#C42221", "#5EA7FF", "#785715", "#0CEA91", "#FFFAED",
            "#B3AF9D", "#3E3D52", "#5A9BC2", "#9C2F90", "#8D5700", "#ADD79C", "#00768B", "#337D00",
            "#C59700", "#3156DC", "#944575", "#ECFFDC", "#D24CB2", "#97703C", "#4C257F", "#9E0366",
            "#88FFEC", "#B56481", "#396D2B", "#56735F", "#988376", "#9BB195", "#A9795C", "#E4C5D3",
            "#9F4F67", "#1E2B39", "#664327", "#AFCE78", "#322EDF", "#86B487", "#C23000", "#ABE86B",
            "#96656D", "#250E35", "#A60019", "#0080CF", "#CAEFFF", "#323F61", "#A449DC", "#6A9D3B",
            "#FF5AE4", "#636A01", "#D16CDA", "#736060", "#FFBAAD", "#D369B4", "#FFDED6", "#6C6D74",
            "#927D5E", "#845D70", "#5B62C1", "#2F4A36", "#E45F35", "#FF3B53", "#AC84DD", "#762988",
            "#70EC98", "#408543", "#2C3533", "#2E182D", "#323925", "#19181B", "#2F2E2C", "#023C32",
            "#9B9EE2", "#58AFAD", "#5C424D", "#7AC5A6", "#685D75", "#B9BCBD", "#834357", "#1A7B42",
            "#2E57AA", "#E55199", "#316E47", "#CD00C5", "#6A004D", "#7FBBEC", "#F35691", "#D7C54A",
            "#62ACB7", "#CBA1BC", "#A28A9A", "#6C3F3B", "#FFE47D", "#DCBAE3", "#5F816D", "#3A404A",
            "#7DBF32", "#E6ECDC", "#852C19", "#285366", "#B8CB9C", "#0E0D00", "#4B5D56", "#6B543F",
            "#E27172", "#0568EC", "#2EB500", "#D21656", "#EFAFFF", "#682021", "#2D2011", "#DA4CFF",
            "#70968E", "#FF7B7D", "#4A1930", "#E8C282", "#E7DBBC", "#A68486", "#1F263C", "#36574E",
            "#52CE79", "#ADAAA9", "#8A9F45", "#6542D2", "#00FB8C", "#5D697B", "#CCD27F", "#94A5A1",
            "#790229", "#E383E6", "#7EA4C1", "#4E4452", "#4B2C00", "#620B70", "#314C1E", "#874AA6",
            "#E30091", "#66460A", "#EB9A8B", "#EAC3A3", "#98EAB3", "#AB9180", "#B8552F", "#1A2B2F",
            "#94DDC5", "#9D8C76", "#9C8333", "#94A9C9", "#392935", "#8C675E", "#CCE93A", "#917100",
            "#01400B", "#449896", "#1CA370", "#E08DA7", "#8B4A4E", "#667776", "#4692AD", "#67BDA8",
            "#69255C", "#D3BFFF", "#4A5132", "#7E9285", "#77733C", "#E7A0CC", "#51A288", "#2C656A",
            "#4D5C5E", "#C9403A", "#DDD7F3", "#005844", "#B4A200", "#488F69", "#858182", "#D4E9B9",
            "#3D7397", "#CAE8CE", "#D60034", "#AA6746", "#9E5585", "#BA6200")
  VisMat<-matrix(1,ncol = length(Clust$Groups), nrow = length(Clust$Groups)+1)
  VisMat[,1]<-1
  for(i in 2:length(Clust$Groups)){
    G1<-cutk(Clust,i-1)
    G2<-cutk(Clust,i)
    G.same<-intersect(G1,G2)
    #G1<-setdiff(G1, G.same)
    G2<-setdiff(G2, G.same)
    small.G2<-which.min(unlist(lapply(G2,length)))
    VisMat[,i]<-VisMat[,i-1]
    VisMat[G2[[small.G2]],i]<-i
  }
  Max.V<-round(rev(Clust$min.value),2)
  if(is.null(Order)==TRUE && is.null(kg)==TRUE) Order<-1:(length(Clust$Groups)+1)
  if(is.null(kg)==FALSE){
    GG<-cutk(Clust,kg)
    lGG<-unlist(lapply(GG, length))
    oGG<-order(lGG,decreasing = TRUE)
    Order<-GG[[oGG[1]]]
    for(j in 2:kg)Order<-c(Order,GG[[oGG[j]]])
  }
  if(is.null(namesX)==TRUE) namesX<-paste0("X",Order) else  namesX<- namesX[Order]
  if(is.null(nplot)==TRUE){
    layout(matrix(c(2,1,1,1,1,1)))
    par(mar=c(4,5,0,2))
    image(x=1:length(Clust$Groups),y=1:(length(Clust$Groups)+1),z=t(VisMat[Order,]),
          col=colPal[1:200],xlab="Number of Clusters", axes=F,
          ylab=" ",cex.lab=1.5, ylim=c(0,(length(Clust$Groups)+4)),xlim=c(0.5,max.kg+0.5))
    axis(1,at=seq(1,length(Clust$Groups),by = max(floor(length(Clust$Groups)/10),1)), cex.axis=1.2,
         labels = seq(1,length(Clust$Groups),by = max(floor(length(Clust$Groups)/10),1)))
    axis(2,at=seq(1,length(Clust$Groups)+1,1), labels = namesX,cex.axis=cex.xaxis,tick=FALSE,las=2)
    par(mar=c(1,5,5,2))
    #(max(Clust$Diss.Matrix)-min(Clust$Diss.Matrix))
    plot(1:length(Clust$Groups),rev(Clust$min.value),type = "p",
         axes=FALSE, main = "Minimum dissimilarity",ylim = c(0,max(1,Clust$min.value)),
         pch=20,col=colorRampPalette(c("blue", "yellow", "red"))(105)[Max.V*100+1],xlab="",ylab=" ",
         xlim=c(1,max.kg))
    segments(1:(length(Clust$Groups)-1),rev(Clust$min.value[-1]),
             x1=2:length(Clust$Groups), y1=rev(Clust$min.value)[-1],
             col=colorRampPalette(c("blue", "yellow", "red"))(105)[Max.V*100+1])
    axis(2)
    axis(1,at=seq(1,length(Clust$Groups),by = max(floor(length(Clust$Groups)/10),1)), cex.axis=1.2,
         labels = seq(1,length(Clust$Groups),by = max(floor(length(Clust$Groups)/10),1)))
  }else{
    if(nplot==1){
      layout(matrix(c(1)))
      par(mar=c(5,5,5,2))
      plot(1:length(Clust$Groups),rev(Clust$min.value),type = "p",
           axes=FALSE, main = "Minimum dissimilarity",ylim = c(0,max(1,Clust$min.value)),
           pch=20,col=colorRampPalette(c("blue", "yellow", "red"))(105)[Max.V*100+1],
           xlim=c(1,max.kg),
           xlab="Number of Clusters",ylab=" ")
      segments(1:(length(Clust$Groups)-1),rev(Clust$min.value[-1]),
               x1=2:length(Clust$Groups), y1=rev(Clust$min.value)[-1],
               col=colorRampPalette(c("blue", "yellow", "red"))(105)[Max.V*100+1])
      axis(2)
      axis(1,at=seq(1,length(Clust$Groups),by = max(floor(length(Clust$Groups)/10),1)), cex.axis=1.2,
           labels = seq(1,length(Clust$Groups),by = max(floor(length(Clust$Groups)/10),1)))
    }else{
      layout(matrix(c(1)))
      par(mar=c(4,5,0,2))
      image(x=1:length(Clust$Groups),y=1:(length(Clust$Groups)+1),z=t(VisMat[Order,]),
            col=colPal[1:200],xlab="Number of Clusters", axes=F,
            ylab=" ",cex.lab=1.5, ylim=c(0,(length(Clust$Groups)+4)),xlim=c(0.5,max.kg+0.5))
      axis(1,at=seq(1,length(Clust$Groups),by = max(floor(length(Clust$Groups)/10),1)), cex.axis=1.2,
           labels = seq(1,length(Clust$Groups),by = max(floor(length(Clust$Groups)/10),1)))
      axis(2,at=seq(1,length(Clust$Groups)+1,1), labels = namesX,cex.axis=cex.xaxis,tick=FALSE,las=2)
    }
  }
  par(opar)          # restore original settings
  }


##Identified Clusters

cutk<-function(Clust,kg=NA,alpha=NA)
{
  #Clust output from HSM
  if(length(na.omit(kg))>0){
    initial_k<-length(Clust$min.value)+1
    GG<-Clust$Groups[[initial_k-kg]]
    return(GG)
  }
  if(length(na.omit(alpha))>0){
    initial_k<-length(Clust$min.value)+1
    kg<-initial_k-(which(Clust$min.value>alpha)[1]-1)
    if(length(na.omit(kg))>0){
      GG<-Clust$Groups[[initial_k-kg]]
      return(GG)
    }
    else{print("One Cluster")}
  }
  if(length(na.omit(alpha))==0 & length(na.omit(kg))==0) print("Give me a value of alpha or kg")
}

boot.clus<-function(X,Clust,kg0,alpha=c(.01,.05,.1),nboot=1000,parallel=FALSE,freq=1,par.spectrum=c(100,1/(2*dt))){
  X<-scale(X,scale=FALSE)
  TV.Observed<-rev(Clust$min.value)[kg0]
  GroupsMembers0<-cutk(Clust,kg0)
  GroupsMembersA<-cutk(Clust,kg0+1)
  testgroupsaux<-matrix(NA,ncol=kg0+1,nrow=1)
  testgroupsaux<-(unlist(lapply(lapply(GroupsMembersA,intersect,y=GroupsMembers0[[kg0]]),length)))/unlist(lapply(GroupsMembersA,length))
  testgroups<- which(testgroupsaux==1,arr.ind=T)
  G1<-as.numeric(na.omit(GroupsMembersA[[testgroups[1]]]))
  G2<-as.numeric(na.omit(GroupsMembersA[[testgroups[2]]]))
  dt<-1/freq
  np<-length(par.spectrum)
  if(np==1){a<-par.spectrum[1]; wn<-1/(2*dt)}
  if(np==2){a<-par.spectrum[1]; wn<-par.spectrum[2]}
  SpecMatrix2<-matrix(0,ncol=length(c(G1,G2)),nrow=length(X[,1])/2)
  for(specj in 1:length(c(G1,G2)))SpecMatrix2[,specj]<-spec.parzen(X[,c(G1,G2)[specj]],a=a,dt=dt,wn=wn,nn=length(X[,1])/2)[,2]
  ww<-spec.parzen(X[,c(G1,G2)[specj]],a=a,dt=dt,wn=wn,nn=length(X[,1])/2)[,1]
  if(length(G1)==1) FG1<-SpecMatrix2[,1]
  if(length(G1)>1) FG1<-rowMeans(SpecMatrix2[,1:length(G1)])
  if(length(G2)==1) FG2<-SpecMatrix2[,-(1:length(G1))]
  if(length(G2)>1) FG2<-rowMeans(SpecMatrix2[,-(1:length(G1))])
  #Bootstrap Test
  f<-(FG1+FG2)/2
  Spec.to.Spec<-function(S,N,dt,ancho,wn){
    Zk<-rnorm(N)
    specZ<-spec.parzen(Zk,a=ancho,dt=dt,w0=(1/N*dt),wn=wn,nn=N/2)
    specsim<-(pi/dt)*S*specZ[,2]
    return(specsim)
  }
  normalize<-function(w,g){
    nor<-((w[2:length(w)]-w[1:(length(w)-1)])%*%(g[2:length(w)]+g[1:(length(w)-1)])/2)
    return(g/nor)
  }

  if(parallel==TRUE){
    library("doParallel")
    cl<-detectCores()-1
    registerDoParallel(cl)
    ResPar<-foreach(icount(nboot),.combine = c) %dopar%
    { #.packages="HSMClust"
      Spec_Boot<-apply(cbind(f,f),2,Spec.to.Spec,N=length(X[,1]),dt=dt,ancho=a,wn=wn)
      Spec_Boot_N<-apply(Spec_Boot,2,normalize,w=ww)
      TVD(ww,Spec_Boot_N[,1],Spec_Boot_N[,2])
    }
    SampleBoot<-ResPar
  }else{
    FMat<-matrix(rep(f,2*nboot),ncol=2*nboot)
    Spec_Boot<-apply(FMat,2,Spec.to.Spec,N=length(X[,1]),dt=dt,ancho=a,wn=wn)
    Spec_Boot_N<-apply(Spec_Boot,2,normalize,w=ww)
    TV_Boot<-rep(NA,nboot)
    for(k in 1:nboot) TV_Boot[k]<-TVD(ww,Spec_Boot_N[,(2*k-1)],Spec_Boot_N[,2*k])
    SampleBoot<-TV_Boot
  }
  Rech<-matrix(NA,ncol=1,nrow=length(alpha))
  for(i in 1:length(alpha))Rech[i,1]<-ifelse(quantile(SampleBoot,probs=1-alpha[i],type=1)<TV.Observed,1,0)
  p.value<-ifelse(length(which(TV.Observed>(sort(SampleBoot))))>0,1-rev(which(TV.Observed>(sort(SampleBoot))))[1]/nboot,0)
  Res<-rbind(Rech,p.value)
  rownames(Res)<-c(paste("Alpha", alpha),"P value")
  return(Res)
}

check.clust<-function(X,Groups,freq=1,Merger=1,par.spectrum=c(100,1/(2*dt),512),iter=FALSE,dif.tol=0)
{
  X<-scale(X,scale=FALSE)
  #Initial parameters
  n<-length(X[,1])
  k<-length(X[1,])
  dt<-1/freq
  ng<- length(Groups)
  #Time<-n*dt
  np<-length(par.spectrum)
  if(np==1){a<-par.spectrum[1]; wn<-1/(2*dt); length.w<-512}
  if(np==2){a<-par.spectrum[1]; wn<-par.spectrum[2]; length.w<-512}
  if(np==3){a<-par.spectrum[1]; wn<-par.spectrum[2]; length.w<-par.spectrum[3]}
  Spec<-matrix(NA,nrow=length.w,ncol=k)
  for(i in 1:k){
    Aux_Spec<-spec.parzen(X[,i],a=a,dt=dt,wn=wn,nn=length.w)
    w<-Aux_Spec[,1]
    Spec[,i]<-Aux_Spec[,2]/((w[2:length.w]-w[1:(length.w-1)])%*%(Aux_Spec[2:length.w,2]+Aux_Spec[1:(length.w-1),2])/2)
  }
  RSpec<-matrix(NA,nrow=length.w,ncol=ng)
  if(Merger==1){
    for(indg in 1:ng){
      RSpec[,indg]<-spec.parzen(c(X[,Groups[[indg]]]),a=a,dt=dt,wn=wn,nn=length.w)[,2]
      RSpec[,indg]<-RSpec[,indg]/((w[2:length.w]-w[1:(length.w-1)])%*%(RSpec[2:length.w,indg]+RSpec[1:(length.w-1),indg])/2)
    }
  }
  if(Merger==2){
    for(indg in 1:ng){
      if(length(Groups[[indg]])==1)RSpec[,indg]<-Spec[,Groups[[indg]]]
      else RSpec[,indg]<-rowMeans(Spec[,Groups[[indg]]])
    }
  }
  TV.to.RSpec.Mat<-matrix(0,ncol=ng,nrow=k)
  for(i in 1:k)for(indg in 1:ng) TV.to.RSpec.Mat[i,indg]<-TVD(w,Spec[,i],RSpec[,indg])
  if(iter==FALSE){
    Group2<-numeric(k)
    for(indg in 1:ng)Group2[Groups[[indg]]]<-indg
    TV.to.RSpec<-list()
    for(i in 1:k){
      TV.to.RSpec[[i]]<-matrix(c(TV.to.RSpec.Mat[i,Group2[i]],TV.to.RSpec.Mat[i,-Group2[i]]),nrow=1)
      colnames(TV.to.RSpec[[i]])<-paste0("Cluster",c(Group2[i],(1:ng)[-Group2[i]]))
    }
    TV.to.RSpec[[k+1]]<-which(apply(TV.to.RSpec.Mat,1,which.min)!=Group2);names(TV.to.RSpec[[k+1]])<-"Wrong Cluster"
    return(TV.to.RSpec)
  }
  if(iter==TRUE){
    Group2<-numeric(k)
    for(indg in 1:ng)Group2[Groups[[indg]]]<-indg
    iter.num<-0
    while(sum(apply(TV.to.RSpec.Mat,1,which.min)!=Group2)>dif.tol){
      iter.num<-1+iter.num
      min.dif<-which(apply(TV.to.RSpec.Mat,1,which.min)!=Group2)
      if(length(min.dif)==1){ new.clust<-which.min(TV.to.RSpec.Mat[min.dif,])
      }else{new.clust<-apply(TV.to.RSpec.Mat[min.dif,],1,which.min)}
      clusters.change<-sort(unique(c(Group2[min.dif],new.clust)))
      Group2[min.dif]<-new.clust
      new.RSpec<-RSpec
      if(Merger==1){
        for(indg in clusters.change){
          new.RSpec[,indg]<-spec.parzen(c(X[,Group2==indg]),a=a,dt=dt,wn=wn,nn=length.w)[,2]
          new.RSpec[,indg]<-new.RSpec[,indg]/((w[2:length.w]-w[1:(length.w-1)])%*%(new.RSpec[2:length.w,indg]+new.RSpec[1:(length.w-1),indg])/2)
        }
      }
      if(Merger==2){
        for(indg in clusters.change){
          if(sum(Group2==indg)==1)new.RSpec[,indg]<-Spec[,Group2==indg]
          else new.RSpec[,indg]<-rowMeans(Spec[,Group2==indg])
        }
      }
      TV.to.RSpec.Mat<-matrix(0,ncol=ng,nrow=k)
      for(i in 1:k)for(indg in 1:ng) TV.to.RSpec.Mat[i,indg]<-TVD(w,Spec[,i],new.RSpec[,indg])
    }
    TV.to.RSpec<-list()
    for(i in 1:k){
      TV.to.RSpec[[i]]<-matrix(c(TV.to.RSpec.Mat[i,Group2[i]],TV.to.RSpec.Mat[i,-Group2[i]]),nrow=1)
      colnames(TV.to.RSpec[[i]])<-paste0("Cluster",c(Group2[i],(1:ng)[-Group2[i]]))
    }
    TV.to.RSpec[["Iteration"]]<-iter.num;
    #names(TV.to.RSpec[[k+1]])<-"Iteration"
    TV.to.RSpec[["New Clusters"]]<-Group2;
    #names(TV.to.RSpec[[k+2]])<-"New Clusters"
    return(TV.to.RSpec)
  }
}

