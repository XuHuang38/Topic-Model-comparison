library(methods)
library(nimble)
library(tidyr)
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
    		bb[k,v]~dgamma(phi[v],0.01)
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

adj1 <- read.csv("~/small_sim_A_0102_adjm.csv")
adj1 = adj1[,-1]
adj1[is.na(adj1)] <- 0

docscombined <- read.csv("~/0102_smA_wordmat_num.csv")
Num = 15
doclengthscomb <- matrix(NA,Num,1)
for(i in 1:Num){
  sum(!is.na(docscombined[i,]))->doclengthscomb[i,]
}
doclengthscomb<-as.vector(doclengthscomb)


data1 <- list(adj=adj1, word=docscombined)

consts <- list(phi= c(rep(0.1,82)), mu=c(2,2), var1=matrix(c(.5,0,0,.5),2,2),V=82,K=4, J=Num, lowera=.1 , uppera=5, lowerb=-20 , upperb=20,  N=doclengthscomb)
#290 for 60 words

model1 <- nimbleModel(testcode, data=data1, constants=consts)

Cmodel1 <- compileNimble(model1)


model1Conf <- configureMCMC(model1,monitors = c('theta', 'beta', 'z', 'zt'),useConjugacy=TRUE, onlySlice=TRUE, print=T)
#
model1MCMC <- buildMCMC(model1Conf)
Cmodel1MCMC <- compileNimble(model1MCMC, project=model1)
niter <- 15000
Cmodel1MCMC$run(niter)




samples <- as.matrix(Cmodel1MCMC$mvSamples)
vars <- Cmodel1MCMC$mvSamples$getVarNames()
col.names <- list()
for(col in colnames(samples)){
  if(grepl("z", col, fixed = TRUE)){  
    col.names<- append(col.names, col)
  }
}
results <- matrix(NA, ncol = 38, nrow = 1)
colnames(results)<- col.names
for(col in colnames(samples)){
  if(grepl("", col, fixed = TRUE)){
    val <- mean(samples[,col])
    results[1, col] <- val
  }
}



write.table(pos.vec, file="0323_sim_A_result_pos.csv",row.names=FALSE,col.names=FALSE, sep=",")
write.table(results, file="0323_sim_A_result.csv",row.names=FALSE,col.names=FALSE, sep=",")


#0202_smsim_A_result.csv
# 0220_sim_A_result.csv


# analysis

## analysis on bb & top <- here top is saved only temporarily.

clus.res <- c(3, 3, 1, 3, 3, 3, 3, 1, 2, 3, 3, 1, 2, 3, 3)
adj1 <- read.csv("small_sim_A_0102_adjm.csv")
adj1 = adj1[,-1]
adj1[is.na(adj1)] <- 0
pos.vec <- read.csv('0220_sim_A_result.csv', header = FALSE)
bbsummary <- array(NA, dim = c(82,4,6))


for(col in colnames(samples)){
  if(grepl("bb", col, fixed = TRUE)){
    print(col)
    matches <- regmatches(col, gregexpr("[[:digit:]]+", col))
    word_top <- as.numeric(unlist(matches))
    val <- (samples[,col])
    bbsummary[word_top[2],word_top[1] , ] <- summary(val)
  }
}

par(mfrow = c(2, 2))
plot(samples[, "bb[1, 1]"],type = "l")
plot(samples[, "bb[2, 1]"],type = "l" )
plot(samples[, "bb[3, 1]"],type = "l")
plot(samples[, "bb[4, 1]"],type = "l")


top_mat <- as.matrix(Cmodel1MCMC$mvSaved)
topic_subm<- array(NA, dim = c(15, 25, 4))
for(col in colnames(top_mat)){
  if(grepl("top[", col, fixed = TRUE)){
    print(col)
    matches <- regmatches(col, gregexpr("[[:digit:]]+", col))
    node_word_top <- as.numeric(unlist(matches))
    val <- (top_mat[,col])
    topic_subm[node_word_top[1],node_word_top[2] , node_word_top[3]] <- val
  }
}


as.matrix(Cmodel1MCMC$mvSaved)[, "top[1, 16, 4]"]




# beta larger <- more to the topic


## analysis and summary on location parameter z.
pos.vec <- matrix(pos.vec, ncol = 2)
plot(pos.vec[,1],pos.vec[,2], xlim=c(1.5,2.5),ylim=c(1.5, 2.5),xlab="X Location",ylab="Y Location",main="Scenario A", col= c('blue','blue','red','blue','blue','blue','blue','red','green','blue','blue','red','green','blue','blue', 'black', 'black', 'black', 'black'))
for(i in 1:15){
  for(j in (i):15){
    if(adj1[i,j] == 1){
      segments(pos.vec[i,1][[1]],pos.vec[i,2][[1]], pos.vec[j,1][[1]],pos.vec[j,2][[1]])
    }
  }
}

loc_sam <- array(NA, dim = c(15, 2, 30000))
for(col in colnames(samples)){
  if(grepl("z[", col, fixed = TRUE)){
    print(col)
    matches <- regmatches(col, gregexpr("[[:digit:]]+", col))
    node_loc <- as.numeric(unlist(matches))
    val <- (samples[,col])
    loc_sam[node_loc[1],node_loc[2], 1:30000] <- val
  }
}

library(rjags)

traceplot(as.mcmc(loc_sam[1,2,]))


## analysis again on beta and theta (representing top)
# 15000 iterations in total
# beta[k,v] <- bb[k,v]/sum(bb[k,1:V]) which makes us want to look w.r.t word, not group.
# i.e. we can have 1 word representing different groups.

betasum_mat <- array(NA, dim = c(82,4,6))
for(col in colnames(samples)){
  if(grepl("beta", col, fixed = TRUE)){
    matches <- regmatches(col, gregexpr("[[:digit:]]+", col))
    word_top <- as.numeric(unlist(matches))
    val <- (samples[,col])
    betasum_mat[word_top[2],word_top[1] , ] <- summary(val)
  }
}  


thetasum_mat <- array(NA, dim = c(15, 4, 6))
for(col in colnames(samples)){
    if(grepl("theta", col, fixed = TRUE)){
      matches <- regmatches(col, gregexpr("[[:digit:]]+", col))
      top <- as.numeric(unlist(matches))
      val <- (samples[,col])
      thetasum_mat[top[1],top[2] , ] <- summary(val)
    }
}
