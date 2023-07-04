##cite https://doi.org/10.1111/j.1469-185X.2012.00240.x
## Understanding phylogenetic incongruence: lessons from phyllostomid bats.
## Liliana M. Dávalos, Andrea L. Cirranello*, Jonathan H. Geisler, Nancy B. Simmons 
## Biological Reviews 87(4), 991-1041

##state to step analysis and plot

##libraries
library(ggplot2)
library(segmented)

##get data
##step category 1 = is end step a new state
## an example of how to extract this from parsimony mapping is in sat_11ascending.xls
ss<-as.data.frame(read.csv("statstep.csv"))

##make into factor
ss$step.category<-as.factor(ss$step.category)

##two alternative functions to linear
##rarefaction or saturation
rare<- function(step,a) {
	a*(1-((1-1/a)^step))
	}

##hypergeometric	
hype <- function(step,b,c) {
	c*(step^b)
	}

##make a linear model # state is a function of # steps
linear<-lm(state~step,data=ss)

##fit a segmented regression to find break
fitseg<-segmented(linear,seg.Z=~step,data=ss)

##parameter corresponding to break, in this case step 879
snap<-as.integer(round(fitseg$psi[2]))

##fit rarefaction function after break	
fitrar<-nls(state~rare(step,a), data=subset(ss,step>snap),start =c(a=100))

##fit hypergeometric function after break
fithyp<-nls(state~hype(step,b,c), data=subset(ss,step>snap),start =c(b=.1,c=1))

##make comparable linear fit
fitlin<-lm(state~step,data=subset(ss,step>snap))

##print out to check best fit
sink("fit_comparisons.txt")
print("AIC linear")
print(AIC(fitlin))
print("AIC rarefaction")
print(AIC(fitrar))
print("AIC hypergeometric")
print(AIC(fithyp))
sink()

##rank test for step category being equal rank after linear phase
mwt<-wilcox.test(step~step.category,data=subset(ss,step>snap))

##print rank test
sink("rank.txt")
print(mwt)
sink()

##plot results	
##predict  values for non linear phase
ss$predrare[(snap+1): dim(ss)[1]]<-predict(fitrar)
ss$predhype[(snap+1): dim(ss)[1]]<-predict(fithyp)

##pull together observed and modeled
y<-c(ss$state,ss$predrare,ss$predhype)
x<-rep(ss$step,3)

##make into dataframe
allp<- as.data.frame(cbind(x,y))
allp$data<- c(rep("observed", dim(ss)[1]),rep("finite", dim(ss)[1]),rep("ordered", dim(ss)[1]))

##plot
ggplot(data=allp, aes(x = x, y = y), na.rm = TRUE) + geom_point( aes(colour = data)) + theme_bw() + scale_colour_manual(values = c("grey50","black", "grey70")) + xlab("No. of steps in MP phylogeny") + ylab("No. of states observed")

##save plot
ggsave("state_by_step.pdf")

##save analysis
save.image("saturation.RData")
