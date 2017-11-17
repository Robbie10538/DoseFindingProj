#DATA
seed=32232
GP=list(G=20000,gamma1=1,gamma2=1,K=5,deltaCRD=-999,deltaCRDMax=999,S=c(1,0.9,0.8,0.75,0.6),
        mu_alpha=c(0,1),mu_beta=c(0),Lambda_alpha=diag(c(0.3,0.3)),Lambda_beta=2,doses=c(0.00,0.05,0.20,0.60,1.00),sigma=3,
        alpha=0.025,N1=c(10,30,60,80,100,120,200),N2=seq(from=0,to=2000,by=200),nSim_from_post=1000,nSim_dec1_mc=500)
sim_gam_from_prior=function(GP){
  c(rmvnorm(1,mean=GP$mu_alpha,sigma=GP$Lambda_alpha),rlnorm(1,meanlog=GP$mu_beta,sdlog=GP$Lambda_beta))
}
nresp_per_dose=5
doses=rep(c(0.00,0.05,0.20,0.60,1.00),each=nresp_per_dose)
set.seed(seed)
gam=sim_gam_from_prior(GP)
mu=model_resp(gam,GP$doses,"emax")
biom=list(d=doses,resp=rep(mu,each=nresp_per_dose)+rnorm(nresp_per_dose*length(c(0.00,0.05,0.20,0.60,1.00)),mean=0,sd=GP$sigma))

data=list(d=sort(unique(biom$d)),Ybar=sapply(sort(unique(biom$d)),function(d){mean(biom$resp[which(biom$d==d)])}))
doses=sort(unique(biom$d))
mu


#prior
prior_mat=matrix(NA,nrow=5000,ncol=length(GP$doses))
for(i in 1:5000){
  gam=sim_gam_from_prior(GP)
  prior_mat[i,]=model_resp(gam,GP$doses,"emax")
}
colMeans(prior_mat)
sapply(1:ncol(prior_mat),function(d){sd(prior_mat[,d])})

#BORNKAMP WAY
modnames=names(Mods(emax=Inf,doses=doses))
prior <- list(list(norm = c(0, .3), norm = c(1,sqrt(2)*0.3), lnorm=c(0,2)))

anMod <- lm(resp~factor(biom$d)-1, data=biom);drFit <- coef(anMod);S <- vcov(anMod)
gsample <- bFitMod(doses, drFit, S, model = modnames[1],
                   MCMCcontrol=list(adapt=T,thin=2),start = rep(0.5,length(prior[[1]])),
                   nSim = 500000, prior = prior[[1]])
gsample
head(gsample$samples)
colMeans((gsample$samples))
predictive_response_at_dose=function(){
  response_mat=matrix(NA,ncol=length(doses),nrow=length(gsample$samples[,1]))
  for(i in 1:length(gsample$samples[,1])){
    response_mat[i,]=model_resp(gamma=as.numeric(gsample$samples[i,]),doses,model=modnames[1])}
  response_mat}
response_mat=predictive_response_at_dose()
colMeans(response_mat)
sapply(1:ncol(response_mat),function(d){sd(response_mat[,d])})
find_n2(sapply(sort(unique(biom$d)),function(d){sum(biom$d==d)}),response_mat,GP)

#MY WAY
Gen_posterior_samples<-function(data,GP,seed=50,nPhII=rep(50,length(data$d)),no_column_grid_pts=5500,Nsamples=50000,plot=T,post_para=F){
  
  #BUILD MARGINAL POSTERIOR MEANS AND COVARIANCES
  Sigma=(GP$sigma^2)*diag(nPhII^-1)
  z1_beta_vec=function(beta){(data$d/(beta+data$d))}
  Xbeta=function(beta){matrix(c(1-z1_beta_vec(beta),z1_beta_vec(beta)),ncol=2,byrow=F)}
  xi=function(beta){t(Xbeta(beta))%*%solve(Sigma)%*%data$Ybar +solve(GP$Lambda_alpha)%*%GP$mu_alpha}
  Delta=function(beta){solve(t(Xbeta(beta))%*%solve(Sigma)%*%Xbeta(beta)+solve(GP$Lambda_alpha))}
  
  #FUNCTION TO RETURN RESPONSE USING DOSE RESPONSE MODEL
  dr_model=function(d,alpha,beta){
    alpha[1] + (alpha[2]-alpha[1])*(d/(beta+d))
  }
  
  #FUNCTION TO RETURN VALUE PROPORTIONAL TO POSTERIOR DISTRIBUTION OF BETA
  beta_post_func=function(beta,data,sigma){
    as.numeric(((det(Delta(beta)))^0.5)*exp(0.5*t(xi(beta))%*%Delta(beta)%*%xi(beta))*dlnorm(beta,meanlog=GP$mu_beta,sdlog=GP$Lambda_beta))
  }
  
  #GENERATE GRID POINTS
  lowerbnd=GP$mu_beta-4*GP$Lambda_beta
  if(lowerbnd>0){grid_pts=seq(from=0.00000001,to=500,length.out=no_column_grid_pts)}else{
    grid_pts=seq(from=0.00000001,to=500,length.out=no_column_grid_pts)
  }
  eval_grid_pt=sapply(grid_pts,beta_post_func,data=data,sigma=GP$sigma)
  heights_columns=sapply(1:(length(eval_grid_pt)-1),function(i){max(eval_grid_pt[i],eval_grid_pt[i+1])*1.1})
  
  #PLOT GRID POINTS: grid points [k] and [k+1] correspond to column [k] and heights_columns[k]
  if(plot==T){
    plot(grid_pts,eval_grid_pt,col=1,pch=18,ylim=c(0,max(max(eval_grid_pt)*1.1,0.1)),xlab="beta",ylab="Posterior density of beta|data")
    invisible(sapply(1:(length(grid_pts)-1),function(i){
      lines(c(grid_pts[-length(grid_pts)][i],grid_pts[-1][i]),c(heights_columns[i],heights_columns[i]),col=2,lty=2);
      lines(c(grid_pts[-length(grid_pts)][i],grid_pts[-length(grid_pts)][i]),c(0,heights_columns[i]),col=2,lty=2);
      lines(c(grid_pts[-1][i],grid_pts[-1][i]),c(0,heights_columns[i]),col=2,lty=2)}))
  }
  
  #SAMPLING SCHEME: Importance sampling (of beta) within gibbs (alpha,beta)
  Post_alphabeta_Mat<-matrix(NA,ncol=3,nrow=Nsamples)
  set.seed(seed)
  for(i_sample in 1:Nsamples){
    accept=0
    while(accept==0){
      #Sample a column proportionally with the height of each column, and a beta uniformly in this column.
      column_id=which(rmultinom(1,1,heights_columns)==1)
      beta_sample=runif(1,grid_pts[column_id],grid_pts[column_id+1])
      #accept/reject
      prob=min(beta_post_func(beta_sample,data,GP$sigma)/heights_columns[column_id],1)
      accept=rbinom(1,1,prob)
      if(prob>1){print("Y value goes above height of column")}
    }
    beta=beta_sample
    #sample E0,Emax
    alpha=as.vector(rmvnorm(1,mean=Delta(beta)%*%xi(beta),sigma=Delta(beta)))
    Post_alphabeta_Mat[i_sample,]=c(alpha,beta)
  }
  
  if(post_para==T){
    return(Post_alphabeta_Mat)
  }else{
    Post_responses_Mat<-matrix(NA,ncol=length(data$d),nrow=Nsamples)
    for(i in 1:Nsamples){
      Post_responses_Mat[i,]=dr_model(data$d,c(Post_alphabeta_Mat[i,c(1,2)]),Post_alphabeta_Mat[i,3])
    }
    return(Post_responses_Mat)
  }
}


post_resp_mat=Gen_posterior_samples(data,GP,seed=100,sapply(sort(unique(biom$d)),function(d){sum(biom$d==d)}),
                                    no_column_grid_pts =15000,Nsamples=10000,plot=T,post_para=F)
colMeans(post_resp_mat)
sapply(1:ncol(post_resp_mat),function(d){sd(post_resp_mat[,d])})
find_n2(sapply(sort(unique(biom$d)),function(d){sum(biom$d==d)}),post_resp_mat,GP)
