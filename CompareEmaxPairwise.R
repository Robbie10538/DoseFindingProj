#Safety S fixed.
#Priors on alpha beta. Sample from alpha,beta. sample data.
#Emax: Observe data, update posterior on alpha beta, choose n2,istar
#Pairwise: Convert priors on alpha beta to mvn prior on theta. Force Sigma to be diagnonal.
#..........Update posterior on theta pairwise. Choose n2,istar.
#Compare E(Gain) from Emax and Pairwise approaches.

Prog3_Sim=function(GP){
  dec1_mc(GP,seeds,parallel=T)
}