library(mvtnorm)
library(truncnorm)
#data
data=list(d=c(0,0.2,0.5,0.8,1),Ybar=c(0,0.3,0.6,0.9,1))
#dr(dose)= alpha[1]+ (alpha[2]-alpha[1])(d/(beta+d))
#alpha~ N(mu_alpha,Nabla_alpha)
#beta ~ N+(mu_beta,Nabla_beta)


Gen_posterior_samples<-function(data,GP,seed=50,nPhII=rep(50,length(data$d)),no_column_grid_pts=100,Nsamples=500,plot=T,post_para=F){
  
#BUILD MARGINAL POSTERIOR MEANS AND COVARIANCES
Sigma=(GP$sigma^2)*diag(nPhII^-1)
z1_beta_vec=function(beta){(data$d/(beta+data$d))}
Xbeta=function(beta){matrix(c(1-z1_beta_vec(beta),z1_beta_vec(beta)),ncol=2,byrow=F)}
xi=function(beta){t(Xbeta(beta))%*%solve(Sigma)%*%data$Ybar +solve(GP$Nabla_alpha)%*%GP$mu_alpha}
Delta=function(beta){solve(t(Xbeta(beta))%*%solve(Sigma)%*%Xbeta(beta)+solve(GP$Nabla_alpha))}

#FUNCTION TO RETURN RESPONSE USING DOSE RESPONSE MODEL
dr_model=function(d,alpha,beta){
  alpha[1] + (alpha[2]-alpha[1])*(d/(beta+d))
}

#FUNCTION TO RETURN VALUE PROPORTIONAL TO POSTERIOR DISTRIBUTION OF BETA
beta_post_func=function(beta,data,sigma){
  as.numeric(((det(Delta(beta)))^0.5)*exp(0.5*t(xi(beta))%*%Delta(beta)%*%xi(beta))*dtruncnorm(beta,a=0,b=Inf,mean=GP$mu_beta,sd=GP$Nabla_beta))}

#GENERATE GRID POINTS
lowerbnd=GP$mu_beta-4*GP$Nabla_beta
if(lowerbnd>0){grid_pts=seq(from=lowerbnd,to=GP$mu_beta+4*GP$Nabla_beta,length.out=no_column_grid_pts)}else{
  grid_pts=seq(from=0.0001,to=GP$mu_beta+5*GP$Nabla_beta,length.out=no_column_grid_pts)
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