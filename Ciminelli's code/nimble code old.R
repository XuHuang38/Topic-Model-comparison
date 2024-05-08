library(methods)
library(nimble)
testcode<-nimbleCode({
	for (j in 1:J) {
		z[j,1:2] ~ dmnorm(mu[1:2], var1[1:2,1:2])
	}
	for (k in 1:K) {
		for(j in 1:J){
	        temp1[j,k] <- pow(zt[k,1]-z[j,1],2)
			temp2[j,k] <- pow(zt[k,2]-z[j,2],2)
			temp3[j,k] <- temp1[j,k]+temp2[j,k]
			temp4[j,k] <- pow(temp3[j,k],0.5)
			temp5[j,k] <- exp(-temp4[j,k])
		}
	}
    for (k in 1:K) {
    	for(v in 1:V){
    		bb[k,v]~dgamma(phi[v],1)
    		beta[k,v] <- bb[k,v]/sum(bb[k,1:V])
		}
		zt[k,1:2] ~ dmnorm(mu[1:2], var1[1:2,1:2])
		for(j in 1:J){
			theta[j,k] <- temp5[j,k]/(temp5[j,1]+temp5[j,2]+temp5[j,3])
		}
    }
	for(j in 1:J){
		for(v in 1:V){
			for(k in 1:K){
				temp6[j,k,v]<-theta[j,k]*beta[k,v]
			}
			klbeta[j,v]<-sum(temp6[j,1:K,v])
		}
	}

	for (j in 1:J) {
    	for (h in 1:J) {
			for(v in 1:V){
				tempkl[j,v,h]<-klbeta[j,v]*log(klbeta[j,v]/klbeta[h,v])
			}
			kl[j,h]<-sum(tempkl[j,1:V,h])
			kl1[j,h]<-kl[j,h]+kl[h,j]
		}
    }
    for (j in 1:J) {
    	for (y in 1:J) {
    		p[j,y] <-a+b*kl1[j,y]-pow(pow((z[j,1]-z[y,1]),2)+pow((z[j,2]-z[y,2]),2),0.5)
    		padj[j,y] <- exp(p[j,y])/(1+exp(p[j,y]))
    		adj[j,y] ~ dbern(padj[j, y])
        }
    }
    for (j in 1:J) {
        for (n in 1:N[j]) {
            top[j,n,1:K] ~ dmulti(theta[j,1:K],1)
            for(k in 1:K){
            		for(v in 1:V){
            			betatemp[k,v] <- pow(beta[k,v],top[j,n,k])
            		}	
            }
            for(v in 1:V){
            		betatemp2[v] <- prod(betatemp[1:K,v])
            }
            word[j,n] ~ dcat(betatemp2[1:V])
        }
    }
    a ~ dunif(lowera, uppera)
    b ~ dunif(lowerb, upperb)
}
)

source("~/Dropbox/Univ/research/prof Love/network-old (fixed seed)-3.R")
data1 <- list(adj=adj1, word=docscombined)

consts <- list(phi= c(rep(0.01,60)), mu=c(0,0), var1=matrix(c(.25,0,0,.25),2,2),V=60,K=3, J=Num, lowera=.1 , uppera=5, lowerb=-20 , upperb=20,  N=doclengthscomb)

model1 <- nimbleModel(testcode, data=data1, constants=consts)
Cmodel1 <- compileNimble(model1)


model1Conf <- configureMCMC(model1,useConjugacy=TRUE, onlySlice=TRUE, print=T)
#model1Conf$addSampler(target="top[1, 1, 1:3]", type= 'slice')
model1Conf$addMonitors(c("theta","beta","z","zt","top"))
model1MCMC <- buildMCMC(model1Conf)
Cmodel1MCMC <- compileNimble(model1MCMC, project=model1)
niter <- 500
Cmodel1MCMC$run(niter)

samples <- as.matrix(Cmodel1MCMC$mvSamples)
vars <- Cmodel1MCMC$mvSamples$getVarNames()
col.names <- list()
for(col in colnames(samples)){
  if(grepl("z", col, fixed = TRUE)){  
    col.names<- append(col.names, col)
  }
}
results <- matrix(NA, ncol = 26, nrow = 1)
colnames(results)<- col.names
for(col in colnames(samples)){
  if(grepl("z", col, fixed = TRUE)){
    val <- mean(samples[,col])
    results[1, col] <- val
  }
}



results_vt <- matrix(NA, ncol = 6, nrow = 1)
col.names_vt <- list()
for(col in colnames(samples)){
  if(grepl("zt", col, fixed = TRUE)){  
    col.names_vt<- append(col.names_vt, col)
  }
}
colnames(results_vt)<- col.names_vt
for(col in colnames(samples)){
  if(grepl("zt", col, fixed = TRUE)){
    val <- mean(samples[,col])
    results_vt[1, col] <- val
  }
}

pos.vec <- results[, 1:20]
pos.vec_vt <- results_vt[, 1:6]
pos.vec_vt <- matrix(pos.vec_vt, ncol = 2)
pos.vec <- matrix(pos.vec, ncol = 2)
plot(pos.vec[,1],pos.vec[,2], xlim=c(-2,2),ylim=c(-2, 2),xlab="X Location",ylab="Y Location",main="Scenario B")#, main="Simulated Social Space")
points(pos.vec_vt, col=c("red","green","blue"),pch=17,cex=1.5)

for(i in 1:10){
  for(j in (i):10){
    if(adj1[i,j]==1){
      segments(pos.vec[i,1],pos.vec[i,2],pos.vec[j,1],pos.vec[j,2])
    }
  }
}


write.csv(results, '1124nimble_result.csv')
