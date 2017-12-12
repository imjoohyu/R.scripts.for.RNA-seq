#Make a word cloud out of functional assignments for core genes
#Dec 11th, 2017
#Joo Hyun Im (ji72)

#Based on https://www.r-bloggers.com/word-clouds-using-text-mining/

rm(list=ls(all=TRUE))
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/")

upregulated = readLines("upregulated_core_functional_catergories.txt")

library(tm)
myCorpus = Corpus(VectorSource(upregulated))

myCorpus = tm_map(myCorpus, tolower)
myCorpus = tm_map(myCorpus, removePunctuation)
myCorpus = tm_map(myCorpus, removeNumbers)
myCorpus = tm_map(myCorpus, removeWords, stopwords("english"))

myDTM = TermDocumentMatrix(myCorpus, control = list(minWordLength = 1))

m = as.matrix(myDTM)

v = sort(rowSums(m), decreasing = TRUE)

library(wordcloud)
set.seed(4363)
wordcloud(names(v), v, min.freq = 50)
