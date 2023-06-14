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
newchar<-sapply(eachdist,function(eachdist) as.vector(eachdist))

##transpose so characters are in rows and make into a dataframe
##distances are calculated between characters
newchar<-as.data.frame(t(newchar))

##calculate Gower distances for real data
newdist<-daisy(newchar, metric= "gower")

##read in variable molecular data
##this has been formatted to be numbers
matm <-read.csv("variable.csv")

##remove the species names
matm1<-matm[,2:4602]

##make into factors
##deprecated indicate in daisy command
#matm1<-as.data.frame(apply(matm1,2,function(matm1) as.factor(matm1)))

##transpose matrix
matm1<-as.data.frame(t(matm1))

##calculate distance from factors
moldist<-daisy(matm1, metric = "gower", type = list(factor = seq(1: dim(matm1)[2] )))

##create a vector of labels
label<-c(rep("molecular",length(moldist)),rep("morphology",length(newdist)))

##pull the distances together
gowdist<-c(moldist, newdist)

##create dataframe for plotting
dist4plot<-as.data.frame(cbind(label,gowdist))

##make gowdist numeric
dist4plot$gowdist<-as.numeric(dist4plot$gowdist)

##order factor in label
dist4plot$label <- ordered(dist4plot$label, levels = c("morphology", "molecular"))

##plot histogram in facets
ggplot(aes(gowdist), data=dist4plot) + geom_histogram(binwidth=0.025, aes(y = after_stat(width*density))) + facet_grid(label~.,scales="free_y") + theme_bw() + xlab("Dissimilarity between characters") + ylab("Frequency") + scale_y_continuous(labels = scales::percent)

##save plot
ggsave("Gower_distances_combined.pdf")

##get the plot and the statistical analysis comparing relative freqs
pdf("relative_density_combined.pdf")
##this estimates and plots the relative frequency distributions of observed and simulated
reldist(y=newdist,yo=moldist,
                      yowgt=rep(1,length(moldist)), ywgt=rep(1,length(newdist)),
                      cdfplot=F, ci=T,
                      method="gam",
                      cex=0.8,bar=T,
                      ylab="Relative density of dissimilarity between observed morphological characters",
                      xlab="Dissimilarity between variable molecular characters")
dev.off()
                      
##save workspace
save.image("Gower_distance_combined.RData")
