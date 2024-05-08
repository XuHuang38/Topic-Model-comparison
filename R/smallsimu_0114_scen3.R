library('Linkage')
#BBC articles. Each topic is associated with a given article. 
#Knowing that a word is in topic say j, the word is sampled uniformly at random from article j. 
#100-150 words for each existing pair of vertices. 

# Topics & Words generation

t1.s <- readLines('~/Dropbox/Univ/research/prof Love/msgA_f.txt')
t2.s <- readLines('~/Dropbox/Univ/research/prof Love/msgB_f.txt')
t3.s <- readLines('~/Dropbox/Univ/research/prof Love/msgC_f.txt')
t4.s <- readLines('~/Dropbox/Univ/research/prof Love/msgD_f.txt')

# Simulation of C
N <- 20
K <- 1:3
Q <- 1:4

# Adj list:
adj_m <- matrix(NA, nrow = N, ncol = N)
group_l <- rep(0, N)
for (i in 1:N) {
  group_i <- sample(Q, size = 1)
  group_l[i] <- group_i
}

prob.con <- c(0.25, 0.01)
for(i in 1:N){
  connect_l <- rep(0, N+1-i)
  group_no <-group_l[i]
  # print("group number of i is: ")
  # print(group_no)
  for(j in i:N){
    group_no_j <- group_l[j]
    # print("group number of j is: ")
    # print(group_no_j)
    u <- runif(1)
    # print("u is: ")
    # print(u)
    if(group_i == group_no_j){
      if(u <= prob.con[1]){
        connect_l[j-i+1] = 1
      }
    }
    else{
      if(u<=prob.con[2]){
        connect_l[j-i+1] = 1
      }
    }
  }
  print(connect_l)
  adj_m[i, i:N] <- connect_l
}
# adj_m[lower.tri(adj_m)] <- t(adj_m)[lower.tri(adj_m)]
to_del = adj_m
to_del[which(is.na(adj_m))] = 0
result_l <- matrix(NA, nrow = sum(to_del), ncol = 3)
count = 1
for(i in 1:N){
  group_i <- group_l[i]
  for(j in i:N){
    if(adj_m[i,j] == 1){
      result_l[count, 1] = i
      result_l[count, 2] = j
      group_j <- group_l[j]
      if((group_i == group_j)  & (group_j == 1)){
        w_loc <- sample(1:length(t1.s), size = 5, replace = TRUE)
        words = t1.s[w_loc]
        result_l[count, 3] = paste(words,collapse = ' ')
      }
      else if((group_i == group_j) & (group_j == 2)){
        w_loc<-sample(1:length(t2.s), size = 5, replace = TRUE)
        words = t2.s[w_loc]
        result_l[count, 3] = paste(words,collapse = ' ')
      }
      else if((group_i == group_j) & (group_j == 3)){
        w_loc<-sample(1:length(t1.s), size = 5, replace = TRUE)
        words = t1.s[w_loc]
        result_l[count, 3] = paste(words,collapse = ' ')
      }
      else if((group_i == group_j) & (group_j == 4)){
        w_loc<-sample(1:length(t2.s), size = 5, replace = TRUE)
        words = t2.s[w_loc]
        result_l[count, 3] = paste(words,collapse = ' ')
      }
      else{
        w_loc <- sample(1:length(t3.s), size = 5, replace = TRUE)
        words = t3.s[w_loc]
        result_l[count, 3] = paste(words,collapse = ' ')
      }
      count = count + 1
    }
  }
}

result_l

write.table(result_l, file="small_sim_C_0114.csv",row.names=FALSE,col.names=FALSE, sep=",")


adj_m[lower.tri(adj_m)] <- t(adj_m)[lower.tri(adj_m)]
write.csv(adj_m, file="small_sim_C_0114_adjm.csv", )


library(tm)
library(topicmodels)
docmat <- matrix(NA, nrow = N, ncol = 1)
for(i in 1:N){
  toAdd <- list(NA)
  for(j in 1:nrow(result_l)){
    #print(j)
    if((as.numeric(result_l[j, 1]) == i) | (as.numeric(result_l[j, 2]) == i)){
      #print(TRUE)
      toAdd <- append(toAdd, values= result_l[j,3])
    }
  }
  #print(length(toAdd))
  toAdd <- toAdd[2:length(toAdd)]
  docmat[i, 1] = paste(toAdd,collapse = ' ')
}
library(stringr)
maxlength <- 0
for(i in 1:N){
  maxlength2 <- length(str_split(docmat[i,1], pattern = " ")[[1]])
  if(maxlength< maxlength2){
    maxlength = maxlength2
  }
}

wordmat <- matrix(NA, nrow = N, ncol = maxlength)
for(i in 1:N){
  maxlength2 <- length(str_split(docmat[i,1], pattern = " ")[[1]])
  wordmat[i,1:maxlength2] <- str_split(docmat[i,1], pattern = " ")[[1]]
}

write.table(wordmat, file="0114_smsim_C_wordmat.csv",row.names=FALSE,col.names=FALSE, sep=",")

result_l

