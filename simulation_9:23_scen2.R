library('Linkage')
#BBC articles. Each topic is associated with a given article. 
#Knowing that a word is in topic say j, the word is sampled uniformly at random from article j. 
#100-150 words for each existing pair of vertices. 

# Topics & Words generation

t1.s <- readLines('msgA_f.txt')
t2.s <- readLines('msgB_f.txt')
t3.s <- readLines('msgC_f.txt')
t4.s <- readLines('msgD_f.txt')

# Simulation of B
N <- 100
K <- 1:3
Q <- 1:2

# Adj list:
adj_m <- matrix(NA, nrow = N, ncol = N)
group_l <- rep(0, 100)
for (i in 1:100) {
  group_i <- sample(Q, size = 1)
  group_l[i] <- group_i
}

prob.con <- 0.25
for(i in 1:100){
  connect_l <- rep(0, 101-i)
  group_no <-group_l[i]
  # print("group number of i is: ")
  # print(group_no)
  for(j in i:100){
    group_no_j <- group_l[j]
    # print("group number of j is: ")
    # print(group_no_j)
    u <- runif(1)
    # print("u is: ")
    # print(u)
    if(group_no == 1){
      if(group_no_j == 1){
        if(u <= prob.con){
          connect_l[j-i+1] = 1
        }
      }
      else{
        if(u <= prob.con){
          connect_l[j-i+1] = 1
        }
      }
    }
    else{
      if(group_no_j == 1){
        if(u <= prob.con){
          connect_l[j-i+1] = 1
        }
      }
      else{
        if(u <= prob.con){
          connect_l[j-i+1] = 1
        }
      }
    }
  }
  print(connect_l)
  adj_m[i, i:100] <- connect_l
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
      if((group_i == group_j) & (group_j == 1)){
        w_loc <- sample(1:length(t2.s), size = 150, replace = TRUE)
        words = t2.s[w_loc]
        result_l[count, 3] = paste(words,collapse = ' ')
      }
      else if((group_i == group_j) & (group_j == 2)){
        w_loc<-sample(1:length(t3.s), size = 150, replace = TRUE)
        words = t3.s[w_loc]
        result_l[count, 3] = paste(words,collapse = ' ')
      }
      else{
        w_loc <- sample(1:length(t4.s), size = 150, replace = TRUE)
        words = t4.s[w_loc]
        result_l[count, 3] = paste(words,collapse = ' ')
      }
      count = count + 1
    }
  }
}

result_l

write.table(result_l, file="925_sim_B_v3.csv",row.names=FALSE,col.names=FALSE, sep=",")
