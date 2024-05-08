# set.seed(2541133) #Used for figured in proposal set 1
set.seed(125678) #Used for figured in proposal set 2
#set.seed(1468) #Used for figured in proposal set 2
#set.seed(841242)
t1.s <- readLines('~/Dropbox/Univ/research/prof Love/msgA_f.txt')
t2.s <- readLines('~/Dropbox/Univ/research/prof Love/msgB_f.txt')
t3.s <- readLines('~/Dropbox/Univ/research/prof Love/msgC_f.txt')
t4.s <- readLines('~/Dropbox/Univ/research/prof Love/msgD_f.txt')
word_pack_1 <- t1.s[sample(1:length(t1.s), size = 20, replace = FALSE)]
word_pack_2 <- t2.s[sample(1:length(t2.s), size = 20, replace = FALSE)]
word_pack_3 <- t3.s[sample(1:length(t3.s), size = 20, replace = FALSE)]
word_pack <- c(word_pack_1,word_pack_2,word_pack_3)

#set.seed(2525451)
library(mvtnorm)
Num <- 10
u <- rbinom(Num,1,.5)
mix1 <- rmvnorm(Num,c(-1,-1),sigma=matrix(c(.5,.1,.1,.5), nrow=2))
mix2 <- rmvnorm(Num,c(3,3),sigma=matrix(c(.6,.2,.2,.6), nrow=2))
points1 <- u*mix1+(1-u)*mix2
#plot(points1, ylab="Y Location", xlab="X Location", main="Social Space Locations", xlim=c(-2,4), ylim=c(-2,4))


dist1<-dist(points1)
 names1<-NULL
 for(i in 1:(Num-1)){
	for(j in (i+1):Num){
		 paste(i,"_",j,sep="")->names
		 names1<-c(names1,names)
	 }
 }
as.matrix(as.vector(dist1))->dist1
rownames(dist1)<-names1

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

theta<-matrix(NA,Num,3)
for(i in 1:Num){
	for(j in 1:3){
		theta[i,j]<-exp(-dist_char_top[i,j])/(sum(exp(-dist_char_top[i,])))
	}	
}

V <- 60
beta <- matrix(NA, 3, V)
library(MCMCpack)
beta[1,] <- rdirichlet(1,c(rep(300,20),rep(.01,20), rep(.01,20)))
beta[2,] <- rdirichlet(1,c(rep(3,20),rep(3,20), rep(.02,20)))
beta[3,] <- rdirichlet(1,c(rep(.03,20),rep(.01,20), rep(3,20)))


dpc<-4
char<-NULL
for(i in 1:Num){
	rep(i,dpc)->l
	char<-c(char,l)
}

sample(seq(20,30,1),length(char),T)->doclengths
matrix(NA, length(char),max(doclengths))->docs
for(j in 1:length(doclengths)){
	for(i in 1:doclengths[j]){
		sample(1:3,1,prob=theta[char[j],])->top
		sample(seq(1:60),1,prob=beta[top,])->word
		docs[j,i]<-word
	}
}

rownames(docs)<-char

docscombined<-matrix(NA, Num,V*dpc)
for(j in 1:Num){
	temp1<-NULL
	for(i in which(rownames(docs)==j)){
		docs[i,which(!is.na(docs[i,]))]->temp
		temp1<-c(temp1,temp)
	}
	temp1->docscombined[j,1:length(temp1)]
}

docscombined[,which(apply(docscombined,2,function(x) sum(is.na(x)))!=Num)]->docscombined

doclengthscomb <- matrix(NA,Num,1)
for(i in 1:Num){
	sum(!is.na(docscombined[i,]))->doclengthscomb[i,]
}
doclengthscomb<-as.vector(doclengthscomb)

odtm1<-matrix(0,Num,max(docscombined,na.rm=T))
for(i in 1:length(docscombined[,1])){
		a<-unique(docscombined[i,])
	for(j in 1:(length(a)-1)){
		if(a[j]>0){
			odtm1[i,a[j]]<-sum(docscombined[i,]==a[j],na.rm=T)
		}
	}
}

klbeta<-matrix(NA,length(theta[,1]),length(beta[1,]))
for(i in 1:length(theta[,1])){
	theta[i,]%*%beta->klbeta[i,]
}

kl<-matrix(NA, length(klbeta[,1]),length(klbeta[,1]))
for(i in 1:length(klbeta[,1])){
	for(j in 1:length(klbeta[,1])){
		if(i==j){
		}else{
			sum(klbeta[i,]*log(klbeta[i,]/klbeta[j,]))->kl[i,j]
		}
	}
}

(kl+t(kl))/2->kl



kllist<-NULL
for(i in 1:(Num-1)){
	for(j in (i+1):Num){
		kl[i,j]->kllist1
		kllist<-c(kllist,kllist1)
	}
}
as.matrix(kllist)->kllist
rownames(kllist)<-names1


alpha<-3
beta1<-10
ad_prob<-exp(alpha-beta1*kllist-dist1)/(1+exp(alpha-beta1*kllist-dist1))
ad_prob
cbind(1-ad_prob,ad_prob)->adprob
edge<-matrix(NA,length(adprob[,1]))
for(i in 1:length(adprob[,1])){
	sample(0:1,prob=adprob[i,],1)->edge[i,]
}

adj<-matrix(0,Num,Num)
l<-1
for(i in 1:(Num-1)){
	for(j in (i+1):Num){
		adj[i,j]<-edge[l]
		l+1->l
	}
}

t(adj)+adj->adj
adj->adj1

col1<-matrix(NA,length(theta[,1]),1)
for(i in 1:length(theta[,1])){
	which(max(theta[i,])==theta[i,])->col
	if(col==1){
		col1[i,]<-"red"
	}else{
		if(col==2){
			col1[i,]<-"green"
		}else{
			if(col==3){
				col1[i,]<-"blue"
			}
		}
	}
}


#plot(points1[1,1],points1[1,2], xlim=c(-2,4),ylim=c(-2,4),col=col1[1,], xlab="X Location",ylab="Y Location", main="Simulated Social Space")
plot(points1[1,1],points1[1,2], xlim=c(-2.5,6),ylim=c(-2.2,5),col=col1[1,], xlab="X Location",ylab="Y Location",main="Scenario 3")#, main="Simulated Social Space")

text(points1[1,1],points1[1,2],label=1,pos=2)
for(i in 2:Num){
	points(points1[i,1],points1[i,2],col=col1[i,])
	text(points1[i,1],points1[i,2],label=i,pos=2)
}


for(i in 1:(Num-1)){
	for(j in (i+1):Num){
		if(adj[i,j]==1){
			segments(points1[i,1],points1[i,2],points1[j,1],points1[j,2])
		}
	}
}

points(topics, col=c("red","green","blue"),pch=17,cex=1.5)
#legend(-2,4, legend=c("Topic 1", "Topic 2", "Topic 3"), pch=17, col=c(2,3,4))
legend(-2,5, legend=c("Topic 1", "Topic 2", "Topic 3"), col=c(2,3,4),pch=17)

for(j in 1:length(kl)){
	if(is.na(kl[j])){
		kl[j]<-0
	}
}

doclengthscomb <- matrix(NA,Num,1)
for(i in 1:Num){
  sum(!is.na(docscombined[i,]))->doclengthscomb[i,]
}
doclengthscomb<-as.vector(doclengthscomb)

word_assign_group <- matrix(NA, ncol = dim(docscombined)[2], nrow = dim(docscombined)[1])
for(i in 1:Num){
  total_group <- sum(adj1[i,])
  for(j in 1: doclengthscomb[i]){
    word_assign_group[i,j] <- sample(total_group, 1, replace = TRUE)
  }
}

result <- matrix(NA, nrow = sum(adj1)/2, ncol = 3)
count <- 1

for(i in 1:Num){
  text_count <- 1
  for (j in i:Num) {
    if(adj1[i,j] == 1){
      result[count, 1] = i
      result[count, 2] = j
      words <- word_assign_group[i, ]
      temp <- which(words == text_count)
      toSave_no <- docscombined[i, temp]
      toSave <- word_pack[toSave_no]
      result[count, 3] = paste(toSave,collapse = ' ')
      count = count + 1
      text_count = text_count + 1
  }
}
}
write.table(result, file="1023_sim_Cimin_v1.csv",row.names=FALSE,col.names=FALSE, sep=",")
write.table(word_pack, file = "1023_sim_cimin_word.csv", row.names = FALSE, col.names=FALSE, sep = ',')



# library(latentnet)
# network(adj,directed=F)->net1
# ergmm(net1~edgecov(kl)+euclidean(d=2),control=ergmm.control(sample.size=1,burnin=0))->ergmm1
# plot(ergmm1)