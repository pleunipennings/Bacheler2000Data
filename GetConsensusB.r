#load the libraries ape and seqinr 
library(ape)
library(seqinr)
#read the fasta file 
consensusfasta<-read.dna("HIV1_CON_2004_POL_DNA.fasta", format = "fasta")
#create a matrix with genetic distances 
dist.dna(consensusfasta)->ConsD
#make a neighbourjoining tree
NJ<-nj(ConsD)
#plot the neighbourjoining tree
#plot(NJ)

consensusB<-consensusfasta[6, 1:984]
