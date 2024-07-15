##                                  AUXILIAR FUNCTIONS
##  These functions are part of the main function saem_zibr. They are for internal use.

# funczibr1: creates the transformation to determine u_it and p_it

funczibr1<-function(psi,dim,id,X)
{
  psim<-psi[id,dim]
  if(length(dim)==1)
    rt<-psim*X
  else
    rt<-rowSums(psim*X)
  fit<-inv.logit(rt)
  return(fit)
}

# diagJB: auxiliar to work with the variance matrix

diagJB<-function(x)
{
  if(length(x)==1) d<-x
  else d<-diag(x)
  return(d)
}

# func_a: first part of the log-likelihood function

func_a<-function(a,psi.cad,ind.a.aleat,xcovM,id2,idg0M,ideq0M,nxcov,nzcov)
{
  ind<-which(ind.a.aleat==F)
  N<-dim(psi.cad)[1]
  psi.cad[,ind]<-matrix(rep(a,each=N),ncol=length(a),nrow=N)
  p<-funczibr1(psi.cad,1:nxcov,id2,xcovM)
  res<-sum(log(1-p[ideq0M]))+sum(log(p[idg0M]))
  return(-res)
}

# func_b: second part of the log-likelihood function

func_b<-function(p,psi.cad,ind.b.aleat,zcovM,id2,idg0M,yobsM,nxcov,nzcov,nal.b)
{
  v<-p[length(p)]
  N<-dim(psi.cad)[1]
  if(nal.b!=nzcov) b<-p[-length(p)] 
  ind<-which(ind.b.aleat==F)
  if(nal.b!=nzcov) psi.cad[,nxcov+ind]<-matrix(rep(b,each=N),ncol=length(b),nrow=N)
  u<-funczibr1(psi.cad,nxcov+(1:nzcov),id2,zcovM)
  res<-sum(lgamma(v)-
             lgamma(v*u[idg0M])-
             lgamma((1-u[idg0M])*v)+
             v*u[idg0M]*log(yobsM[idg0M])+
             v*(1-u[idg0M])*log(1-yobsM[idg0M]))
  return(-res)
}

# llis.zibr: calculates the approximate log-likelihood by importance sampling

llis.zibr<-function(MU,G,V,yobs,idg0,ideq0,xcov,zcov,nxcov,nzcov,id,psi.mean,psi.var,ind.psi.aleat,naleat,niter=100)
{
  require(boot)
  
  Ginv<-diag(G)^-1
  Gdet<-prod(diag(G))
  
  psiM<-array(rep(psi.mean,niter),dim=c(dim(psi.mean),niter))
  sdM<-array(rep(sqrt(psi.var),niter),dim=dim(psiM))
  r<-array(rt(prod(dim(psiM)),5),dim=dim(psiM))
  
  mat<-psiM+sdM*r
  
  pM<-apply(mat,3,funczibr1,1:nxcov,id,xcov)
  uM<-apply(mat,3,funczibr1,nxcov+(1:nzcov),id,zcov)
  
  aux<-matrix(0,ncol=niter,nrow=length(id))
  aux[ideq0,]<-log(1-pM[ideq0,])
  aux[idg0,]<-log(pM[idg0,])+log(yobs[idg0])*(V*uM[idg0,]-1)+log(1-yobs[idg0])*(V*(1-uM[idg0,])-1)-lbeta(V*uM[idg0,],V*(1-uM[idg0,]))
  P1<-apply(aux,2,function(x) tapply(x,id,sum))
  
  aux<-matrix(rep(MU[ind.psi.aleat],niter),ncol=niter)
  P2<-apply(mat[,ind.psi.aleat,],1,function(x) -0.5*(colSums(Ginv*(x-aux)^2)+log(Gdet)+naleat*log(2*pi)))
  
  aux<-dt(r[,ind.psi.aleat,],df=5,log=T)-log(sdM[,ind.psi.aleat,])
  P3<-apply(aux,3,rowSums)

  mf<-P1+t(P2)-P3
  ll<-sum(log(rowMeans(exp(mf))))
  
  return(ll)
}

# zibr.grad: calculates the gradient of the likelihood function of the ZIBR model

zibr.grad<-function(MU,G,V,psi.cad,ind.psi.aleat,ind.a.aleat,ind.b.aleat,naleat,
                   xcovM,id2,idg0M,ideq0M,nxcov,nzcov,zcovM,yobsM,nal.a,nal.b)
{
  G_diag<-diagJB(G)
  N<-dim(psi.cad)[1]
  A<-MU[1:nxcov]
  B<-MU[nxcov+(1:nzcov)]
  S<-colSums(psi.cad[,ind.psi.aleat])
  S2<-colSums(psi.cad[,ind.psi.aleat]^2)

  mu_g<-rep(0,length(MU))
  mu_g[ind.psi.aleat]<-(S-N*MU[ind.psi.aleat])/G_diag
  G_g<-0.5*((S2-2*MU[ind.psi.aleat]*S+N*(MU[ind.psi.aleat]^2))/G_diag^2-N/G_diag)

  if(nal.a!=nxcov)
  {
    A_g<-(-1)*grad(func_a,A[-ind.a.aleat],psi.cad=psi.cad,ind.a.aleat=ind.a.aleat,xcovM=xcovM,id2=id2,idg0M=idg0M,ideq0M=ideq0M,
              nxcov=nxcov,nzcov=nzcov)
  }
  if(nal.b!=nzcov) paru<-c(B[-ind.b.aleat],V)
  else paru<-V

  paru_g<-(-1)*grad(func_b,paru,psi.cad=psi.cad,ind.b.aleat=ind.b.aleat,zcovM=zcovM,id2=id2,idg0M=idg0M,yobsM=yobsM,
                   nxcov=nxcov,nzcov=nzcov,nal.b=nal.b)
  V_g<-paru_g[length(paru)]

  if(nal.b!=nzcov) B_g<-paru_g[-length(paru)]
  if (naleat!=(nxcov+nzcov)) mu_g[-ind.psi.aleat]<-c(A_g,B_g)
  g<-c(mu_g,V_g,G_g)
  return(g)
}

# zibr.hess: calculates the hessian of the likelihood function of the ZIBR model

zibr.hess<-function(MU,G,V,psi.cad,ind.psi.aleat,ind.a.aleat,ind.b.aleat,naleat,
                      xcovM,id2,idg0M,ideq0M,nxcov,nzcov,zcovM,yobsM,nal.a,nal.b)
{
  G_diag<-diagJB(G)
  N<-dim(psi.cad)[1]
  A<-MU[1:nxcov]
  B<-MU[nxcov+(1:nzcov)]
  S<-colSums(psi.cad[,ind.psi.aleat])
  S2<-colSums(psi.cad[,ind.psi.aleat]^2)
  nhess<-nxcov+nzcov+naleat+1
  H<-matrix(0,ncol=nhess,nrow=nhess)

  diag(H)[ind.psi.aleat]<- -N/G_diag
  diag(H)[(nhess-naleat+1):nhess]<-(1/G_diag^2)*(0.5*N-(S2-2*MU[ind.psi.aleat]*S+N*(MU[ind.psi.aleat])^2)/G_diag)
  diag(H[ind.psi.aleat,(nhess-naleat+1):nhess])<-(-1/G_diag^2)*(S-N*MU[ind.psi.aleat])
  diag(H[(nhess-naleat+1):nhess,ind.psi.aleat])<-(-1/G_diag^2)*(S-N*MU[ind.psi.aleat])

  if(nal.a!=nxcov)
  {
    A_h<-(-1)*hessian(func_a,A[-ind.a.aleat],psi.cad=psi.cad,ind.a.aleat=ind.a.aleat,xcovM=xcovM,id2=id2,idg0M=idg0M,ideq0M=ideq0M,
                       nxcov=nxcov,nzcov=nzcov)
    indx1<-setdiff(1:nxcov,ind.a.aleat)
    H[indx1,indx1]<-A_h
  }
  if(nal.b!=nzcov) paru<-c(B[-ind.b.aleat],V)
  else paru<-V
  paru_h<-(-1)*hessian(func_b,paru,psi.cad=psi.cad,ind.b.aleat=ind.b.aleat,zcovM=zcovM,id2=id2,idg0M=idg0M,yobsM=yobsM,
                        nxcov=nxcov,nzcov=nzcov,nal.b=nal.b)
  indx2<-setdiff(1:(nzcov+1),ind.b.aleat)+nxcov
  H[indx2,indx2]<-paru_h
  return(H)
}

            
##                             MAIN FUNCTION
# saem_zibr: estimate the parameters of the ZIBR model with the SAEM algorithm

saem_zibr<-function(Y,X=NULL,Z=NULL,index,v0,a0,b0,seed,iter,ncad=5,a.fix=NULL,b.fix=NULL)
{
  require(MASS)
  require(boot)
  require(numDeriv)
  set.seed(seed)
  
  Stk<-floor(0.75*iter)
  Gtk<-iter-Stk
  
  # Initialization
  
  yobs<-Y
  nind<-length(unique(index))
  ntot<-length(index)
  xcov<-cbind(rep(1,ntot),X)
  zcov<-cbind(rep(1,ntot),Z)
  nxcov<-ncol(xcov)
  nzcov<-ncol(zcov)
  
  if(is.null(colnames(X)))
    labs.X<-c("Intercept",paste("X",1:(nxcov-1),sep="."))[1:nxcov]
  else
    labs.X<-c("Intercept",colnames(X))
  if(is.null(colnames(Z)))
    labs.Z<-c("Intercept",paste("Z",1:(nxcov-1),sep="."))[1:nzcov]
  else
    labs.Z<-c("Intercept",colnames(Z))
  
  idg0<-(yobs!=0)
  ideq0<-(yobs==0)
  npsi<-length(c(a0,b0))
  ind.a.aleat<-if(is.null(a.fix)) c(T,rep(F,nxcov-1)) else a.fix==0
  ind.b.aleat<-if(is.null(b.fix)) c(T,rep(F,nzcov-1)) else b.fix==0
  ind.psi.aleat<-which(c(ind.a.aleat,ind.b.aleat)==T)
  naleat<-length(ind.psi.aleat)
  nal.a<-sum(ind.a.aleat)
  nal.b<-sum(ind.b.aleat)
  id<-as.numeric(factor(index, levels=unique(index)))
  id1<-rep(id,ncad)
  id2<-id1+(nind*(rep(1:ncad,each=ntot)-1))
    
  graph_str<-NULL

  MU<-c(a0,b0)
  Gfull<-0.5*diag(abs(MU))
  G<-Gfull[ind.psi.aleat,ind.psi.aleat]
  V<-v0
  
  xcovM<-do.call("rbind",replicate(ncad,xcov,F))
  zcovM<-do.call("rbind",replicate(ncad,zcov,F))
  idg0M<-rep(idg0,ncad)
  ideq0M<-rep(ideq0,ncad)
  yobsM<-rep(yobs,ncad)
  PSI_CAD<-matrix(rep(MU,ncad*nind),nrow=ncad*nind,ncol=npsi,byrow = T)
  P_CAD<-funczibr1(PSI_CAD,1:nxcov,id2,xcovM)
  U_CAD<-funczibr1(PSI_CAD,nxcov+(1:nzcov),id2,zcovM)
  
  SD1<-0.5*G
  SD1[SD1<0.5&SD1>0]<-0.5
  SD2<-SD1

  Gk<-rep(0,npsi+naleat+1)
  Dk<-matrix(0,ncol=npsi+naleat+1,nrow=npsi+naleat+1)
  Tk<-matrix(0,ncol=npsi+naleat+1,nrow=npsi+naleat+1)
  Hk<-matrix(0,ncol=npsi+naleat+1,nrow=npsi+naleat+1)
  
  psik<-apply(PSI_CAD,2,function(x){tapply(x,rep(1:nind,ncad),mean)})
  Epsik2<-apply(PSI_CAD^2,2,function(x){tapply(x,rep(1:nind,ncad),mean)})
  Sk1<-nind*MU
  Sk2<-Gfull
  
  for(i in 1:iter)
  {
   if(i<=Stk) gam=1
   else gam=1/(i-Stk)
   
   mum<-matrix(rep(MU,ncad*nind),nrow=ncad*nind,ncol=npsi,byrow = T)
   Ginv<-diagJB(diag(G)^-1)
   dsc<-rep(0,ncad*ntot)
   
   # Kernel 1: A priori distribution
   for(j in 1:4)
   {
       psi_c<-mvrnorm(nind*ncad,MU,Gfull)
       psi_c[,-ind.psi.aleat]<-PSI_CAD[,-ind.psi.aleat]
       p_c<-funczibr1(psi_c,1:nxcov,id2,xcovM)
       u_c<-funczibr1(psi_c,nxcov+(1:nzcov),id2,zcovM)
       dsc[ideq0M]<-log(1-P_CAD[ideq0M])-log(1-p_c[ideq0M])
       dsc[idg0M]<-(lgamma(V*u_c[idg0M])+lgamma(V*(1-u_c[idg0M]))+log(P_CAD[idg0M])-
                      lgamma(V*U_CAD[idg0M])-lgamma(V*(1-U_CAD[idg0M]))-log(p_c[idg0M])+
                      V*(U_CAD[idg0M]-u_c[idg0M])*(log(yobsM[idg0M])-log(1-yobsM[idg0M])))
       DD<-as.vector(tapply(dsc,id2,sum))
       RD<-(DD< -log(runif(nind*ncad)))
       PSI_CAD<-psi_c*RD+PSI_CAD*(!RD)
       P_CAD<-funczibr1(PSI_CAD,1:nxcov,id2,xcovM)
       U_CAD<-funczibr1(PSI_CAD,nxcov+(1:nzcov),id2,zcovM)
   }
   #Kernel 2: Unidimensional random walk in random index
   t2<-n2<-rep(0,naleat)
   for(j in 1:4)
   {
       dpsi<-matrix(0,nrow=nind*ncad,ncol=naleat)
       dpsi[matrix(c(1:(ncad*nind),sample(1:naleat,ncad*nind,replace=T)),nrow=ncad*nind)]<-rnorm(ncad*nind)
       psi_c<-PSI_CAD
       psi_c[,ind.psi.aleat]<-PSI_CAD[,ind.psi.aleat]+dpsi%*%SD1
       p_c<-funczibr1(psi_c,1:nxcov,id2,xcovM)
       u_c<-funczibr1(psi_c,nxcov+(1:nzcov),id2,zcovM)
       dsc[ideq0M]<-log(1-P_CAD[ideq0M])-log(1-p_c[ideq0M])
       dsc[idg0M]<-(lgamma(V*u_c[idg0M])+lgamma(V*(1-u_c[idg0M]))+log(P_CAD[idg0M])-
                      lgamma(V*U_CAD[idg0M])-lgamma(V*(1-U_CAD[idg0M]))-log(p_c[idg0M])+
                      V*(U_CAD[idg0M]-u_c[idg0M])*(log(yobsM[idg0M])-log(1-yobsM[idg0M])))
       d1<-psi_c[,ind.psi.aleat]-mum[,ind.psi.aleat]
       d2<-PSI_CAD[,ind.psi.aleat]-mum[,ind.psi.aleat]
       DD<-as.vector(tapply(dsc,id2,sum))+0.5*(diag(d1%*%Ginv%*%t(d1))-diag(d2%*%Ginv%*%t(d2)))
       RD<-(DD< -log(runif(nind*ncad)))
       PSI_CAD<-psi_c*RD+PSI_CAD*(!RD)
       P_CAD<-funczibr1(PSI_CAD,1:nxcov,id2,xcovM)
       U_CAD<-funczibr1(PSI_CAD,nxcov+(1:nzcov),id2,zcovM)
       t2<-t2+colSums((dpsi*RD)!=0)
       n2<-n2+colSums(dpsi!=0)
   }
   SD1<-(1+0.4*(t2/n2-0.4))*SD1
   # Kernel 3: Multidimensional random walk
   t3<-0
   for(j in 1:4)
   {
       psi_c<-PSI_CAD
       psi_c[,ind.psi.aleat]<-PSI_CAD[,ind.psi.aleat]+matrix(rnorm(ncad*nind*naleat),nrow=ncad*nind,ncol=naleat)%*%SD2
       p_c<-funczibr1(psi_c,1:nxcov,id2,xcovM)
       u_c<-funczibr1(psi_c,nxcov+(1:nzcov),id2,zcovM)
       dsc[ideq0M]<-log(1-P_CAD[ideq0M])-log(1-p_c[ideq0M])
       dsc[idg0M]<-(lgamma(V*u_c[idg0M])+lgamma(V*(1-u_c[idg0M]))+log(P_CAD[idg0M])-
                      lgamma(V*U_CAD[idg0M])-lgamma(V*(1-U_CAD[idg0M]))-log(p_c[idg0M])+
                      V*(U_CAD[idg0M]-u_c[idg0M])*(log(yobsM[idg0M])-log(1-yobsM[idg0M])))
       d1<-psi_c[,ind.psi.aleat]-mum[,ind.psi.aleat]
       d2<-PSI_CAD[,ind.psi.aleat]-mum[,ind.psi.aleat]
       DD<-as.vector(tapply(dsc,id2,sum))+0.5*(diag(d1%*%Ginv%*%t(d1))-diag(d2%*%Ginv%*%t(d2)))
       RD<-(DD< -log(runif(nind*ncad)))
       PSI_CAD<-psi_c*RD+PSI_CAD*(!RD)
       P_CAD<-funczibr1(PSI_CAD,1:nxcov,id2,xcovM)
       U_CAD<-funczibr1(PSI_CAD,nxcov+(1:nzcov),id2,zcovM)
       t3<-t3+sum(RD)
   }
   SD2<-(1+0.4*(t3/(4*ncad*nind)-0.4))*SD2
   
   psik<-psik+gam*(apply(PSI_CAD,2,function(x){tapply(x,rep(1:nind,ncad),mean)})-psik)
   Epsik2<-Epsik2+gam*(apply(PSI_CAD^2,2,function(x){tapply(x,rep(1:nind,ncad),mean)})-Epsik2)
   
   Sk1<-Sk1+gam*((colSums(PSI_CAD)/ncad)-Sk1)
   Sk2<-Sk2+gam*(((t(PSI_CAD)%*%PSI_CAD)/ncad)-Sk2)
   
   if(i>10)
   {
     MU<-Sk1/nind
     Gfull<-Sk2/nind-(Sk1%*%t(Sk1))/(nind^2)
     G<-as.matrix(Gfull[ind.psi.aleat,ind.psi.aleat])
     
     A<-MU[1:nxcov]
     B<-MU[nxcov+(1:nzcov)]
       
     if(nal.a!=nxcov)
     {
       ak<-nlminb(start=A[-ind.a.aleat],objective=func_a,psi.cad=PSI_CAD,ind.a.aleat=ind.a.aleat,xcovM=xcovM,id2=id2,idg0M=idg0M,ideq0M=ideq0M,
                   nxcov=nxcov,nzcov=nzcov)$par
       A[-ind.a.aleat]<-A[-ind.a.aleat]+gam*(ak-A[-ind.a.aleat])
     }
       
     if(nal.b!=nzcov) paru<-c(B[-ind.b.aleat],V)
     else paru<-V
       
     paruk<-nlminb(start=paru,objective=func_b,psi.cad=PSI_CAD,ind.b.aleat=ind.b.aleat,zcovM=zcovM,id2=id2,idg0M=idg0M,yobsM=yobsM,
                    nxcov=nxcov,nzcov=nzcov,nal.b=nal.b,lower=c(rep(-Inf,length(paru)-1),0.0001))$par
     paru<-paru+gam*(paruk-paru)
       
     V<-paru[length(paru)]
     if(nal.b!=nzcov) B[-ind.b.aleat]<-paru[-length(paru)]
       
     MU<-c(A,B)
     Gfull[,-ind.psi.aleat]<-0
     Gfull[-ind.psi.aleat,]<-0
       
     PSI_CAD[,-ind.psi.aleat]<-matrix(rep(MU[-ind.psi.aleat],each=ncad*nind),ncol=npsi-naleat,nrow=ncad*nind)

     ZGk<-zibr.grad(MU,G,V,PSI_CAD,ind.psi.aleat,ind.a.aleat,ind.b.aleat,naleat,
                      xcovM,id2,idg0M,ideq0M,nxcov,nzcov,zcovM,yobsM,nal.a,nal.b)
     ZHk<-zibr.hess(MU,G,V,PSI_CAD,ind.psi.aleat,ind.a.aleat,ind.b.aleat,naleat,
                      xcovM,id2,idg0M,ideq0M,nxcov,nzcov,zcovM,yobsM,nal.a,nal.b)
     ZTk<-matrix(0,ncol=npsi+naleat+1,nrow=npsi+naleat+1)
     for(cont in 1:ncad)
     {
       id3<-nind*(cont-1)+(1:nind)
       g2<-zibr.grad(MU,G,V,PSI_CAD[id3,],ind.psi.aleat,ind.a.aleat,ind.b.aleat,naleat,
                       xcov,id,idg0,ideq0,nxcov,nzcov,zcov,yobs,nal.a,nal.b)
       ZTk<-ZTk+g2%*%t(g2)
     }
       
     Gk<-Gk+gam*(ZGk-Gk)
     Dk<-Dk+gam*(ZHk-Dk)
     Tk<-Tk+gam*(ZTk-Tk)
     Hk<-(1/ncad)*(Dk+Tk)-(1/ncad^2)*(Gk%*%t(Gk))
   }  
   
   graph_str<-rbind(graph_str,c(MU,diag(G),V))
   colnames(graph_str)<-c(paste("A",1:nxcov,sep="."),
                          paste("B",1:nzcov,sep="."),
                          paste("SIGMA",(1:length(ind.psi.aleat)),sep="."),
                          "V")
    
  }
  
  psi.mean<-apply(PSI_CAD,2,function(x){tapply(x,rep(1:nind,ncad),mean)})
  psi.var<-Epsik2-psik^2
  psi.var[,-ind.psi.aleat]<-0
  
  loglik.is<-llis.zibr(MU,G,V,yobs,idg0,ideq0,xcov,zcov,nxcov,nzcov,id,psi.mean,psi.var,ind.psi.aleat,naleat,niter=500)
  
  l<-list(MU=MU,G=G,V=V,psi.mean=psi.mean,psi.var=psi.var,loglik=loglik.is,
          graph=graph_str,labs.X=labs.X,labs.Z=labs.Z,nxcov=nxcov,nzcov=nzcov,ind.a.aleat=ind.a.aleat,ind.b.aleat=ind.b.aleat,
          FIM.stoch=Hk,var.stoch=diag(solve(-Hk)))
  class(l)<-c("list","SAEM_ZIBR_result")
    
  return(l)
}

## plot.SAEM_ZIBR: plot method
            
plot.SAEM_ZIBR_result<-function(l)
{
  graph_str<-l$graph
  iter<-nrow(graph_str)
  Stk<-floor(0.75*iter)
  g1<-ncol(graph_str)
  g2<-ceiling(g1/3)
  
  par(mfrow=c(g2,3))
  for(i in 1:g1)
  {
    plot(1:iter,graph_str[,i],xlab=colnames(graph_str)[i],ylab="Value",type="l")
    abline(v=Stk)
  }
  par(mfrow=c(1,1))
}

## print.SAEM_ZIBR: print method

print.SAEM_ZIBR_result<-function(l)
{
  nxcov<-l$nxcov
  nzcov<-l$nzcov
  A<-l$MU[1:nxcov]
  B<-l$MU[nxcov+(1:nzcov)]
  Varest<-diag(l$G)
  ind.a.aleat<-l$ind.a.aleat
  ind.b.aleat<-l$ind.b.aleat
  na.al<-sum(ind.a.aleat)
  nb.al<-sum(ind.b.aleat)
  labs.X<-l$labs.X
  labs.Z<-l$labs.Z
  var.stoch<-l$var.stoch
  
  mlog<-data.frame(Estimate=A,"Std.Error"=0,"p.value"=0,Type=ifelse(ind.a.aleat,"Random","Fixed"),Variance=0,"sqrt(Var)"=0,
                   row.names = labs.X)
  suppressWarnings({mlog[,"Std.Error"]=ifelse(var.stoch[1:nxcov]>0,sqrt(var.stoch[1:nxcov]),NA)}) 
  mlog[,"p.value"]=ifelse(var.stoch[1:nxcov]>0,pchisq((A^2)/var.stoch[1:nxcov],1,lower.tail=F),NA)
  mlog[ind.a.aleat,"Variance"]=Varest[1:na.al]
  mlog[,6]=sqrt(mlog$Variance)
  mbeta<-data.frame(Estimate=B,"Std.Error"=0,"p.value"=0,Type=ifelse(ind.b.aleat,"Random","Fixed"),Variance=0,"sqrt(Var)"=0,
                    row.names = labs.Z)
  suppressWarnings({mbeta[,"Std.Error"]=ifelse(var.stoch[nxcov+(1:nzcov)]>0,sqrt(var.stoch[nxcov+(1:nzcov)]),NA)}) 
  mbeta[,"p.value"]=ifelse(var.stoch[nxcov+(1:nzcov)]>0,pchisq((B^2)/var.stoch[nxcov+(1:nzcov)],1,lower.tail=F),NA)
  mbeta[ind.b.aleat,"Variance"]=Varest[na.al+(1:nb.al)]
  mbeta[,6]=sqrt(mbeta$Variance)
  
  cat("===== ESTIMATION RESULTS SAEM-ZIBR =====",'\n')
  cat("== Logistic part ==",'\n')
  print(mlog[,1:4])
  cat("== Beta part ==",'\n')
  print(mbeta[,1:4])
  cat("=== Variance of random effects ===",'\n')
  cat("== Logistic part ==",'\n')
  print(mlog[ind.a.aleat,5:6])
  cat("== Beta part ==",'\n')
  print(mbeta[ind.b.aleat,5:6])
  cat("=== Phi estimate (beta part): ",l$V,'\n')
  cat("=== Log-likelihood (importance sampling): ",l$loglik,'\n')
  return()
}

## EXAMPLES

#data.sim <- simulate_zero_inflated_beta_random_effect_data(
#  subject_n=100,time_n=5,
#  X = as.matrix(c(rep(0,50*5),rep(1,50*5))),
#  Z = as.matrix(c(rep(0,50*5),rep(1,50*5))),
#  alpha = as.matrix(c(-0.5,0.5)),
#  beta = as.matrix(c(-0.5,0.5)),
#  s1 = 0.43,s2 = 0.76,
#  v = 17.2,
#  sim_seed=332)
#

#data.zibr<-data.frame(ID=data.sim$subject_ind,
#                      Y=data.sim$Y,
#                      X=as.vector(data.sim$X),
#                      Z=data.sim$Z)
#            
# RES1<-saem_zibr(data.zibr$Y,data.zibr$X,data.zibr$Z,data.zibr$ID,15,c(-0.3,0.5),c(-0.2,0.8),ncad=10,232,500)
# RES2<-saem_zibr(data.zibr$Y,data.zibr$X,data.zibr$Z,data.zibr$ID,15,c(-0.4,0.9),c(-0.2,0.8),ncad=10,232,500,a.fix=c(0,0),b.fix = c(0,0))
