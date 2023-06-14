##cite https://doi.org/10.1093/sysbio/syu022
## Integrating Incomplete Fossils by Isolating Conflicting Signal in Saturated and Non-Independent Morphological Characters
##Liliana M. Dávalos, Paúl M. Velazco, Omar M. Warsi, Peter D. Smits, Nancy B. Simmons 
##Systematic Biology, Volume 63, Issue 4, July 2014, Pages 582–600

##dissimilarity analysis including data simulation and relative frequency distribution comparison between observed and simulated dissimilarities

##libraries
library(cluster)
library(coin)
library(ggplot2)
library(reldist)
library(reshape2)
library(ape)
library(readxl)

##load simulated data
load("simdata.RData")

##this is a character state matrix formatted so states can be read
##no multistate, missing data as blank
morp<-read.csv("noct_morp.csv")

##take species names as rownames
rownames(morp)<-morp$species

##delete species names from matrix
morp$species<-NULL

##calculate distances between species based on each character separately
eachdist<-apply(morp, 2, dist, simplify=F)

##turn distance matrices between species into vectors in a dataframe
##characters are now comparable across one another
##insight from Nicole Gerardo
newchar<-sapply(eachdist,function(eachdist) as.vector(eachdist))

##transpose so characters are in rows and make into a dataframe
##distances are calculated between characters
newchar<-as.data.frame(t(newchar))

##calculate Gower distances for real data
newdist<-daisy(newchar, metric= "gower")

##make list for many distances from simulation
sim2com<- list()

##break data down into characters
##initializes a list of lists
sim2com<-sapply(sidat, function(sidat) for(i in 1:dim(sidat)[2]){sidat[,i]})

##generates list of lists with characters
for(k in 1:length(sidat)){for(i in 1:dim(sidat[[1]])[2]){sim2com[[k]][[i]] <- sidat[[1]][,i]}}

##initialize another list of lists
mandist<- sim2com
##translate each char into a distance, this makes chars comparable across one another
for(k in 1:length(sim2com)){for(i in 1:dim(sidat[[k]])[2]) {mandist[[k]][[i]] <- dist(as.numeric(as.factor(sim2com[[k]][[i]])))} }

##initialize list of lists to get characters back
manchar<-sim2com
##make matrices of distance into vectors
for(k in 1:length(mandist)) { for(i in 1:length(mandist[[k]]) ) { manchar[[k]][[i]] <- as.vector(mandist[[k]][[i]])  }}

##now get those vectors and turn back into matrix like datasets
##as many rows as there are chars in simulated/original data which matches the arrangement of newchar
mandat<-sapply(manchar, function(manchar) data.frame(matrix(unlist(manchar), nrow = dim(sidat[[1]])[2], byrow=T)), simplify=F )

##get Gower distance for all these data sets
man2dis<-sapply(mandat, function(mandat) daisy(mandat, metric= "gower"), simplify=F)

##pull together Gower distance values 
##create a vector of labels for observed and simulated datasets
label<-c(rep("observed",length(newdist)),rep("simulated",length(man2dis)*length(man2dis[[1]]) ))
##make a vector with all the distances
gowdist<-c(as.vector(newdist),as.vector(unlist(man2dis)))
##make a dataframe of both vectors
disp<-as.data.frame(cbind(label,gowdist))

##get rid of missing data
disp<-na.omit(disp)

##make Gower distances numeric (issue after cbind)
disp$gowdist<-as.numeric(as.character(disp$gowdist))

##plot histogram in facets
ggplot(aes(gowdist), data=disp) + geom_histogram(binwidth=0.025, aes(y = after_stat(width*density))) + facet_grid(label~.,scales="free_y") + theme_bw() + xlab("Dissimilarity between morphological characters") + ylab("Frequency") + scale_y_continuous(labels = scales::percent)

##save plot
ggsave("Gower_distances_morphology.pdf")

##get the plot and the statistical analysis comparing relative freqs
pdf("relative_density_morphology.pdf")
##this estimates and plots the relative frequency distributions of observed and simulated
reldist(y=newdist,yo=unlist(man2dis),
                      yowgt=rep(1,length(man2dis) * length(man2dis[[1]])),ywgt=rep(1,length(newdist)),
                      cdfplot=F, ci=T,
                      method="gam",
                      cex=0.8,bar=T,
                      ylab="Relative density of dissimilarity between observed morphological characters",
                      xlab="Dissimilarity between simulated morphological characters")
dev.off()
                      
##now generate a printable dataset 
##make a massive matrix                  
temp<-matrix(unlist(man2dis), nrow = length(man2dis[[1]]), ncol = length(man2dis), byrow = T)

##Empirical Cumulative Distribution Function of each simulated dataset by rows
emp2cdf<-apply(temp, 1, ecdf)

##make new list of ranks
rank<-list()

##Empirical Cumulative Distribution Function relative to simulated data
for(i in 1:length(emp2cdf)) { rank[[i]]<-emp2cdf[[i]](newdist[[i]]) }

##observed distances
newmat <- as.matrix(newdist)

##melt into 3 columns: characters and distance between them
m2 <- melt(newmat)[melt(upper.tri(newmat))$value,]

##rename columns
colnames(m2)[1:3]<-c("char1", "char2", "Observed Gower")

##order as they correspond in matrix
m2<-m2[order(m2$char1, m2$char2),]

##insert ranks
m2$rank<-unlist(rank)

##print out for sharing
write.csv(m2,"table_Gower_distances_morphology.csv", row.names=F)

##save workspace
save.image("Gower_distance_morphology.RData")
