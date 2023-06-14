##cite https://doi.org/10.1093/sysbio/syu022
## Integrating Incomplete Fossils by Isolating Conflicting Signal in Saturated and Non-Independent Morphological Characters
##Liliana M. Dávalos, Paúl M. Velazco, Omar M. Warsi, Peter D. Smits, Nancy B. Simmons 
##Systematic Biology, Volume 63, Issue 4, July 2014, Pages 582–600

## data simulation 

##libraries
library(cluster)
library(coin)
library(ggplot2)
library(reldist)
library(reshape2)
library(ape)
library(readxl)


##simulate data
##get phylogenies based on morph characters
##from same run as the rate calculation (post burnin)
tree<-read.nexus("noct_morp32.nex.run1.t")

##scaling of tree using rate of change
##This meant that the branch lengths of each phylogeny used as simulation input had to be multiplied by the rate of character evolution so that the rate of change in the corresponding transition matrix equaled 1. 
##this is just an example to simplify but multiple runs can be used
p<-read.delim("noct_morp32.nex.run1.p",comment.char="[")

##get rid of burnin
tree<-tree[1002:20001]
p<-p[1002:20001,]

##sample rows from posterior (e.g., 100, 500)
##using 10 because life is short, paper used 500 and from 6 different runs
N<-10
rows<-sample(length(tree), N)

##rescale and sample trees
##rate is on column 5
##without rescaling the branch lengths are very short and simulated chars become quite uniform
for (i in rows){tree[[i]]$edge.length<-tree[[i]]$edge.length*p[i,5]}

##sample trees only the ones rescaled
tre1<-tree[rows]

##get number of char states per column
##the file morp_char_desc.xlsx contains the breakdown of character ordering and models
##the vector k comes from the char states (no NA) observed in the data, but accounting for how they were coded not just the simplified table morp
desc<-read_excel("morp_char_desc.xlsx")
k<-desc$'no. char. states'

##read in the model, see morp_char_desc.xlsx for building these matrices depending on ordering
##model here depends on the ordering of characters
##each row is the state x state matrix that defines transitions between states
##this can be replaced with rTraitDisc ARD etc. depending on how characters states are ordered or not
model<-read.csv(file="model.csv")

##blank list
simmor<-list()
sidat<-list()

##generate chars based on the same tree no need to reorder them if generated with same tree
##these data can also be used in RAxML to calculate per char log likelihoods of simulation vs. observed
##this step takes the most time when generating lots of data sets
for (j in 1:length(tre1)) {
for (i in 1:length(k)) {simmor[[i]]<-rTraitDisc(tre1[[j]], model =  matrix(as.numeric(as.character(model[i,1:k[i]^2])), k[i]), k=k[i])}

#merge all chars together and convert into dataframe #
sidat[[j]]<-data.frame(matrix(unlist(simmor), nrow=length(simmor[[i]]), ncol = length(k)))

#make rownames from the order in which they were generated#
rownames(sidat[[j]])<-names(simmor[[i]])
}

##save workspace
save.image("simdata.RData")
