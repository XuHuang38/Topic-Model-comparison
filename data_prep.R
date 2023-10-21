library(tm)
library(topicmodels)
t1 <- readLines('~/Dropbox/Univ/research/prof Love/msgA.txt')
t2 <- readLines('~/Dropbox/Univ/research/prof Love/msgB.txt')
t3 <- readLines('~/Dropbox/Univ/research/prof Love/msgC.txt')
t4 <- readLines('~/Dropbox/Univ/research/prof Love/msgD.txt')

data_prep <- function(toWork, name){
  toWork <- strsplit(toWork, " ")[[1]]
  toWork <- removeWords(toWork, stopwords('en'))
  toWork<-gsub("<ef>","", toWork)
  toWork<-gsub("<cd>","", toWork)
  toWork<-gsub("<cf>","", toWork)
  toWork<-gsub("<d5>","", toWork)
  toWork <- gsub("[[:space:]]+",' ',toWork)
  toWork <- gsub("[[:punct:]]+","",toWork)
  toWork <- sapply(toWork, function(row) iconv(row, "latin1","ASCII", sub=""))
  toWork <- sapply(toWork, function(row) iconv(row, "UTF-8","ASCII", sub=""))
  toWork <- gsub("[\r\n]","",toWork)
  toWork <- tolower(toWork)
  toWork <- gsub("[[:digit:]]","",toWork)
  toWork <- removeWords(toWork, stopwords('en'))
  toWork <- removeWords(toWork, c("a","an","and","im","cant","said","can","may", "upon","will","come","came","like","thing","now","frodo","bilbo","gandalf","merri","sam","pippin"))
  toWork <- removeWords(toWork, letters)
  toWork <- gsub("[[:space:]]+",' ',toWork)
  toWork <- toWork[which(toWork != "")]
  fileConn<-file(name)
  writeLines(toWork, fileConn)
  close(fileConn)
  return(toWork)
}

t2.f <- data_prep(t2, "msgB_f.txt")
t1.f <- data_prep(t1, "msgA_f.txt")
t3.f <- data_prep(t3, "msgC_f.txt")
t4.f <- data_prep(t4, "msgD_f.txt")




