#library(R2jags)
library(mcmcse)
library(abind)
func1<- function(){
  for (j in 1:J) {
    z[j, 1:2] ~ dmnorm(mu[], var1[, ])
  }
  for (k in 1:K) {
    for(j in 1:J){
      temp1[j, k] <- pow(zt[k,1] - z[j,1], 2)
      temp2[j,k] <- pow(zt[k,2] - z[j,2],2)
      temp3[j,k] <- temp1[j,k]+temp2[j,k]
      temp4[j,k] <- pow(temp3[j,k],.5)
      temp5[j,k] <- exp(-temp4[j,k])
    }
  }
  for (k in 1:K) {
    for(v in 1:V){
      bb[k,v]~dgamma(phi[v],1)
      # weird
      beta[k, v] <- bb[k,v]/sum(bb[k,1:V])
    }
    zt[k, 1:2] ~ dmnorm(mu[], var1[, ])
    for(j in 1:J){
      theta[j,k] <- temp5[j,k]/(temp5[j,1]+temp5[j,2]+temp5[j,3])
    }
  }
  
  for(j in 1:J){
    for(v in 1:V){
      for(k in 1:K){
        temp6[j,k,v]<-theta[j,k]*beta[k,v]
      }
      klbeta[j,v]<-sum(temp6[j,,v])
    }
  }
  
  for (j in 1:J) {
    for (h in 1:J) {
      for(v in 1:V){
        tempkl[j,v,h]<-klbeta[j,v]*log(klbeta[j,v]/klbeta[h,v])
      }
      kl[j,h]<-sum(tempkl[j,,h])
      kl1[j,h]<-kl[j,h]+kl[h,j]
    }
  }
  for (j in 1:J) {
    for (y in 1:J) {
      p[j, y] <-a+b*kl1[j,y]-pow(pow((z[j, 1] - z[y, 1]),2)+pow((z[j, 2] - z[y, 2]),2),.5)
      padj[j, y] <- exp(p[j, y])/(1 + exp(p[j, y]))
      adj[j, y] ~ dbern(padj[j, y])
    }
  }
  for (j in 1:J) {
    for (n in 1:N[j]) {
      top[j, n] ~ dcat(theta[j,1:K])
      word[j, n] ~ dcat(beta[top[j, n], 1:V])
    }
  }
  a ~ dunif(lowera, uppera)
  b ~ dunif(lowerb, upperb)
}

docscombined <- read.csv('~/1014_B_wordmat_num.csv')
adj1 <- read.csv("~/1014_sim_B_v3_adjm.csv")
adj1 <- adj1[, 2:length(adj1)]

doclengthscomb <- matrix(NA,100,1)
for(i in 1:100){
  sum(!is.na(docscombined[i,]))->doclengthscomb[i,]
}
doclengthscomb<-as.vector(doclengthscomb)

data1 <- list(adj=adj1, word=docscombined, phi= c(rep(0.01,636)), mu=c(0,0), var1=matrix(c(.25,0,0,.25),2,2), N=doclengthscomb, V=636,K=3, J=100, lowera=.1 , uppera=5, lowerb=-20 , upperb=20)
parameters1<-c("theta","beta","zt","z")

library(R2WinBUGS)
testfile<-file.path(tempdir(), "testfile.txt")
write.model(func1, testfile)

library(mvtnorm)
Num <- 100
u <- rbinom(Num,1,.5)
mix1 <- rmvnorm(Num,c(-1,-1),sigma=matrix(c(.5,.1,.1,.5), nrow=2))
mix2 <- rmvnorm(Num,c(3,3),sigma=matrix(c(.6,.2,.2,.6), nrow=2))
points1 <- u*mix1+(1-u)*mix2
#plot(points1, ylab="Y Location", xlab="X Location", main="Social Space Locations", xlim=c(-2,4), ylim=c(-2,4))


topics<-matrix(NA,3,2)
topics[1,]<-c(-1.3,-1.3)
topics[2,]<-c(1.5,1.5)
#topics[2,]<-c(5.5,5)
topics[3,]<-c(3.5,3)
#topics[2,]<-c(1,.5)
#topics[3,]<-c(2.5,2.5)
#points(topics, col=c(2,3,4), pch=17)
rownames(topics)<-c("t1","t2","t3")
t(as.matrix(dist(rbind(topics,points1)))[1:3,-(1:3)])->dist_char_top
#rownames(dist_char_top)<-c("c1","c2","c3","c4","c5","c6","c7","c8","c9","c10")

inits1<- list( zt=topics, z=points1)

library(rjags)
jags.model(testfile, data1, inits=inits1)->jags1
update(jags1,n.iter=500)
jags.samples(jags1, parameters1, n.iter=5000)->jags2
jags2$theta->thetajags
ess<-multiESS(t(apply(thetajags[,1:2,,1], 3, matrix, ncol=1)))
if(ess<400){
	jags.samples(jags1, parameters1, n.iter=1000)->jags2a
	abind(jags2$beta,jags2a$beta, along=3)->jags_beta
	abind(jags2$theta,jags2a$theta, along=3)->jags_theta
	abind(jags2$z,jags2a$z, along=3)->jags_z
	abind(jags2$zt,jags2a$zt, along=3)->jags_zt
}else{
	jags2$beta-> jags_beta
	jags2$theta-> jags_theta
	jags2$z-> jags_z
	jags2$zt-> jags_zt
}
jagsfinal<-list(beta=jags_beta, theta=jags_theta, z=jags_z,zt=jags_zt)
	

#save(jags2, file="jagssamplesnew.Rda")