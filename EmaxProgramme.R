#global parameters
#install.packages("DoseFinding")
source('MCPModSims.R', echo=F)
GP=list(G=20000,gamma1=1,gamma2=1,K=5,deltaCRD=-999,deltaCRDMax=999,S=c(1,0.9,0.8,0.75,0.6),
        mu_alpha=c(0,1),mu_beta=c(1),Nabla_alpha=diag(c(0.3,0.75)),Nabla_beta=0.5,doses=c(0,0.05,0.2,0.6,1),sigma=3,
        alpha=0.025,N1=c(10,30,60,80,100,120,200),N2=seq(from=0,to=2000,by=200),nSim_from_post=500,nSim_dec1_mc=100)

######################### FUNCTIONS ################
#Calculate Expected Gain given I1
ExpGain_func<-function(mu,istar,nPhII,n2,prej,GP){
  theta=mu[istar]-mu[1]
  GP$G*zeta(theta,istar,GP$deltaCRD,GP$deltaCRDMax)*prej-GP$gamma1*sum(nPhII)-2*GP$gamma2*n2
}  
 
#Specify zeta function
zeta<-function(theta_istar,istar,deltaCRD,deltaCRDMax){
  if(theta_istar<=deltaCRD){z=0}
  if(theta_istar>deltaCRD && theta_istar<deltaCRDMax){z=1}
  if(theta_istar>=deltaCRDMax){z=deltaCRDMax/deltaCRD}
  S=GP$S[istar]
  z*S
}

#simulate model parameters from their prior distributions
sim_gam_from_prior=function(GP){
  c(rmvnorm(1,mean=GP$mu_alpha,sigma=GP$Nabla_alpha),rtruncnorm(1,a=0,b=Inf,mean=GP$mu_beta,sd=GP$Nabla_beta))
}

#plot the prior distribution of responses
plot_prior<-function(seed){
  set.seed(seed)
  gam=sim_gam_from_prior(GP)
  plot(GP$doses,sapply(GP$doses,model_resp,gamma=c(gam[1],gam[2]-gam[1],gam[3]),"emax"),type="l",ylim=c(-2,2),pch="*",cex=.5,xlab="dose",ylab="Response",main=c("emax prior"))
  for(i in 1:100){
    gam=sim_gam_from_prior(GP)
    points(GP$doses,sapply(GP$doses,model_resp,gamma=c(gam[1],gam[2]-gam[1],gam[3]),"emax"),type="l",ylim=c(-1,2),pch="*",cex=.5,xlab="dose",ylab="Response",main=c("emax prior"))
  }
}

#plot the posterior distribution of responses
plot_post<-function(data,post_resp_mat){
  plot(GP$doses,post_resp_mat[1,],type="p",ylim=c(-1,2),pch="*",cex=.5,xlab="dose",ylab="Response",main=c("emax post"))
  for(i in sample(1:length(post_resp_mat[,1]),100,F)){
    points(GP$doses,post_resp_mat[i,],type="p",ylim=c(-1,2),pch="*",cex=.5,xlab="dose",ylab="Response",main=c("emax post"))
  }
  points(data$d,data$Ybar,col=2,pch=17)
  points(data$d,mu,col=3,pch=18)
  abline(h=GP$deltaCRD,lty=2,col=8)
  abline(h=GP$deltaCRDMax,lty=3,col=8)
  abline(h=0,lty=4,col=8)
}

#create data
create_data_func=function(mu,nPhII,seed,GP){
  set.seed(seed)
  noise=rnorm(length(GP$doses))
  data=list(d=GP$doses,Ybar=mu+(GP$sigma/sqrt(nPhII))*noise)
  return(data)
}

#Generate posterior responses samples using rejection sampler within Gibbs sampler.
source('RejSampler.R', echo=F)

#expected gain given mu, n2 and istar. return expected gain and prej
exp_gain_giventheta<-function(n2,nPhII,mu,istar,GP){
  prej=prej_giventheta(n2,mu,istar,GP)
  exp_gain=ExpGain_func(mu,istar,nPhII,n2,prej,GP)
  return(list(exp_gain=exp_gain,prej=prej))
}

#prej given theta
prej_giventheta<-function(n2,mu,istar,GP){
  if(n2==0){return(0)}else{
    theta=mu[istar]-mu[1]
    return(pnorm(theta*sqrt((n2)/(2*GP$sigma^2))-qnorm(1-GP$alpha)))}
}

#decision 2: find maximising n2
find_n2<-function(nPhII,post_resp_mat,GP,clever_search=T){
  eNPV_mat<-matrix(NA,ncol=GP$K,nrow=length(GP$N2));prej_mat=matrix(NA,ncol=GP$K,nrow=length(GP$N2))
  #clever search which only calculates enough to find optimal n2.
  if(clever_search==T && length(GP$N2)>2){
    for(j in 2:GP$K){
      #initial search
      init=floor(length(GP$N2)/2)
      for(i in c(init,init+1)){
        comp=exp_gain_givenPhII(GP$N2[i],j,nPhII,post_resp_mat,GP)
        eNPV_mat[i,j]=comp$exp_gain
        prej_mat[i,j]=comp$prej
      }
      #choose direction
      if(eNPV_mat[init,j]<eNPV_mat[init+1,j]){increment=1}else{increment=-1;i=i-1}
      #keep going until change in increment
      stop=F
      while(stop==F){
        i=i+increment
        comp=exp_gain_givenPhII(GP$N2[i],j,nPhII,post_resp_mat,GP)
        eNPV_mat[i,j]=comp$exp
        prej_mat[i,j]=comp$prej
        if(eNPV_mat[i,j]<eNPV_mat[i-increment,j]){stop=T}
        if(i+increment>length(GP$N2) | i+increment<1){stop=T}
      }
    }
    
  }else{
    for(i in 1:length(GP$N2)){
      for(j in 2:GP$K){
        comp=exp_gain_givenPhII(GP$N2[i],j,nPhII,post_resp_mat,GP)
        eNPV_mat[i,j]=comp$exp_gain
        prej_mat[i,j]=comp$prej
      }
    }
  }
  colnames(eNPV_mat)=as.character(paste("istar=",1:GP$K,sep=""))
  colnames(prej_mat)=as.character(paste("istar=",1:GP$K,sep=""))
  rownames(eNPV_mat)=as.character(paste("n2=",GP$N2,sep=""))
  rownames(prej_mat)=as.character(paste("n2=",GP$N2,sep=""))
  
  opt_indices=as.vector(which(eNPV_mat==max(eNPV_mat,na.rm=T),arr.ind=T))
  opt_n2=GP$N2[opt_indices[1]];opt_istar=opt_indices[2]
  return(list(eNPV_mat=eNPV_mat,prej_mat=prej_mat,opt_n2=opt_n2,opt_istar=opt_istar,opt_eNPV=max(eNPV_mat,na.rm=T),opt_prej=prej_mat[opt_indices[1],opt_indices[2]],opt_indices=opt_indices))
}

#find_n2_pairwise<-function(nPhII,data,GP){
#  thetahat=sapply(setdiff(GP$doses,0),function(dose_var){mean(data[which(data[,1]==dose_var),2])-mean(data[which(data[,1]==0),2])})
#  istar<-which(thetahat==max(thetahat))
  #?????
  #need the prior specified as a mvn, not paramters of a emax model.
#  }

#mc evaluation of expected gain given posterior responses from post_resp_mat, n2 and istar.
exp_gain_givenPhII<-function(n2,istar,nPhII,post_resp_mat,GP){
  no_mcsims=nrow(post_resp_mat)
  exp_gain_vec<-vector()
  prej_vec<-vector()
  for(mcsim in 1:no_mcsims){
    compute=exp_gain_giventheta(n2,nPhII,post_resp_mat[mcsim,],istar,GP)
    exp_gain_vec[mcsim]=compute$exp_gain
    prej_vec[mcsim]=compute$prej
  }
  return(list(exp_gain=mean(exp_gain_vec),prej=mean(prej_vec)))
}

prog_sim<-function(mu,nPhII,seed,GP,model=T,plot=F){
  
  data=create_data_func(mu,nPhII,seed,GP)
  post_resp_mat=Gen_posterior_samples(data,GP,seed=seed,nPhII,no_column_grid_pts = 100,Nsamples=GP$nSim_from_post,plot=T)
  
  if(plot==T){par(mfrow=c(1,2))
    plot_prior(seed);plot_post(data,post_resp_mat)
    par(mfrow=c(1,1))
  }
  
  PhIIIdesign=find_n2(nPhII,post_resp_mat,GP)
  
  compute_actual_theta=exp_gain_giventheta(PhIIIdesign$opt_n2,nPhII,mu,PhIIIdesign$opt_istar,GP)
  opt_eNPV_actual_theta=compute_actual_theta$exp_gain
  opt_prej_actual=compute_actual_theta$prej
  
  return(list(data=data,mu=mu,Ybaristar=data$Ybar[PhIIIdesign$opt_istar],post_resp_mat=post_resp_mat,opt_eNPV=PhIIIdesign$opt_eNPV,opt_prej=PhIIIdesign$opt_prej,opt_eNPV_actual_theta=opt_eNPV_actual_theta,opt_prej_actual=opt_prej_actual,opt_n2=PhIIIdesign$opt_n2,opt_istar=PhIIIdesign$opt_istar,opt_indices=PhIIIdesign$opt_indices,
              eNPV_mat=PhIIIdesign$eNPV_mat,prej_mat=PhIIIdesign$prej_mat,seed=seed))
}

PhII_sample_size_rule<-function(n1,GP){
  rep(n1,GP$K)
}

dec1_mc<-function(GP,seeds,parallel=T){
  ResultsMat=data.frame()
  eNPV_vec=vector()
  
  library(parallel)
  sim_theta_prog_sim=function(GP,i_sim,seeds){
    seed=seeds[i_sim]
    set.seed(seed)
    gam=sim_gam_from_prior(GP)
    mu=gam[1] + (gam[2]-gam[1])*(GP$doses/(gam[3]+GP$doses))
    nPhII=PhII_sample_size_rule(GP$N1[i_n1],GP)
    sim=prog_sim(mu,nPhII,seed=seeds[i_sim],GP,plot=F)
    ResultsMatAdd=data.frame(n1=GP$N1[i_n1],mu1=round(mu[1],2),mu2=round(mu[2],2),mu3=round(mu[3],2),mu4=round(mu[4],2),mu5=round(mu[5],2),muhat_istar=round(sim$Ybaristar,2),opt_eNPV=round(sim$opt_eNPV,5),opt_n2=sim$opt_n2,eNPV_true=round(sim$opt_eNPV_actual_theta,5),prej_true=round(sim$opt_prej_actual,3),opt_istar=sim$opt_istar,seed=seed)
  return(ResultsMatAdd)
    }
  
  
  
  for(i_n1 in 1:length(GP$N1)){print(paste("Trying n1t= ",GP$N1[i_n1]))
    
    if(parallel==T){
    
    cl=makeCluster(detectCores())
    ex <- Filter(function(x) is.function(get(x, .GlobalEnv)), ls(.GlobalEnv))
    clusterExport(cl,c("GP","seeds","i_n1",ex),envir=environment())
    clusterEvalQ(cl, c(library(mvtnorm),library(Bolstad2),library(DoseFinding),library(truncnorm))) #add all the packages you need
    ResultsMatAdd<-t(parSapply(cl=cl,1:GP$nSim_dec1_mc,sim_theta_prog_sim,GP=GP,seeds=seeds))
    stopCluster(cl)
    ResultsMat=rbind(ResultsMat,ResultsMatAdd)
    
    }else{
    
    for(i_sim in 1:GP$nSim_dec1_mc){print(i_sim)
      seed=seeds[i_sim]
      set.seed(seed)
      gam=sim_gam_from_prior(GP)
      mu=model_resp(gam,GP$doses,"emax")
      nPhII=PhII_sample_size_rule(GP$N1[i_n1],GP)
      sim=prog_sim(mu,nPhII,seed=seeds[i_sim],GP,plot=F)
      ResultsMatAdd=data.frame(n1=GP$N1[i_n1],mu1=round(mu[1],2),mu2=round(mu[2],2),mu3=round(mu[3],2),mu4=round(mu[4],2),mu5=round(mu[5],2),opt_eNPV=round(sim$opt_eNPV,5),opt_n2=sim$opt_n2,eNPV_true=round(sim$opt_eNPV_actual_theta,5),prej_true=round(sim$opt_prej_actual,3),opt_istar=sim$opt_istar,seed=seed)
      ResultsMat=rbind(ResultsMat,ResultsMatAdd)
      }
    } 
    }
  eNPV=sapply(GP$N1,function(n1_arg){mean(as.numeric(subset(ResultsMat,n1==n1_arg)$eNPV_true))})
  sd=sapply(GP$N1,function(n1_arg){sd(as.numeric(subset(ResultsMat,n1==n1_arg)$eNPV_true))})
  plot(GP$N1,eNPV)
  
  pairwise_se_mat=matrix(NA,ncol=length(GP$N1),nrow=length(GP$N1))
  for(n1_1 in 1:length(GP$N1)){
    for(n1_2 in 1:n1_1){
  pairwise_se_mat[n1_1,n1_2]=sqrt(var( as.numeric(subset(ResultsMat,n1==GP$N1[n1_1])$eNPV_true)-
                                         as.numeric(subset(ResultsMat,n1==GP$N1[n1_2])$eNPV_true)  )/GP$nSim_dec1_mc)
    }
  }
  savelist=list(ResultsMat=ResultsMat,eNPV=eNPV,sd=sd,pairwise_se_mat=pairwise_se_mat)
  save(savelist,file="Sims.Rdata")
  return(savelist)
}

##############################################

seeds=sample(1:100000,10000,replace=T)
dec=dec1_mc(GP,seeds,T)
dec


#emax parameter priors to mvn prior on theta
if(22==33){
set.seed(seed)
theta_mat=matrix(NA,ncol=length(GP$doses),nrow=100)
gam_mat<-matrix(NA,ncol=3,nrow=100)
for(i in 1:100){
  gam=sim_gam_from_prior(GP)
  gam_mat[i,]=gam
  theta_mat[i,]=model_resp(gam,GP$doses,"emax")
}
mu=colMeans(theta_mat)
Sigma=cov(theta_mat)
mvnorm_resp=rmvnorm(100,mean=mu,sigma=Sigma)

plot(GP$doses,theta_mat[1,],ylim=c(-2,2),type="l")
for(i in 2:100){lines(GP$doses,theta_mat[i,],ylim=c(-2,2),type="l",col=i)}
plot(GP$doses,mvnorm_resp[1,],ylim=c(-2,2),type="l")
for(i in 2:100){lines(GP$doses,mvnorm_resp[i,],ylim=c(-2,2),type="l",col=i)}
}
