#########################################
#                                       #
#   Title: Fit Poisson models to genet  # 
#   and ramet count data from 10        #
#   Long Island and NYC, NY Populations #
#   Distributions for MLL               #
#   Size Distributions                  # 
#   Author: Matthew B. Hamilton         #
#   Date last edited: 26 June 2025      #
#   Manuscript: S. patens landscape     # 
#   genetics - Long Island & NYC        #
#                                       #
#########################################

# --------------------------------------

# Load rstudioapi package to be able to get directory of source file
library("rstudioapi")

setwd(dirname(getActiveDocumentContext()$path)) # Set working directory to source file location

getwd() # Check updated working directory

dat <- read.csv("4.Long_Island_MLLs.csv")

# subset of data
tmp1 <- dat[,c("pop","ramet.count","obs.genets")] 

max.ramet.count <- max(tmp1$ramet.count)

n.pops <- max(tmp1$pop)

# create data frame with max.ramet.count ramets for every population
tmp2 <- data.frame(pop=rep(1:n.pops,rep(max.ramet.count,n.pops)), ramet.count=rep(1:max.ramet.count,n.pops))

# merge data sets
dat2 <- merge(x=tmp2,y=tmp1,by=c("pop","ramet.count"),all.x=TRUE)

# fill in the NA's with 0's
dat2$obs.genets[is.na(dat2$obs.genets)] <- 0

# ANOTHER WAY to VIEW DATA 
# =========================
M <- as.table(matrix(dat2$obs.genets,nrow=n.pops,ncol=max.ramet.count,byrow=TRUE))
dimnames(M) <- list(pop=1:n.pops,ramet.count=1:max.ramet.count)

Xsq <- chisq.test(M) # -->> leads to error message
#fisher.test(M,simulate.p.value=TRUE,B=1000000) # p=0.02393 => reject the null hypothesis of independence
fisher.test(M,simulate.p.value=TRUE,B=100000) # can run smaller reps just to see

simplerate.mat <- matrix(NA, n.pops, 2)
yhat.simple <- matrix(NA, n.pops, max.ramet.count)
yhat.ols <- matrix(NA, n.pops, max.ramet.count)

# pool observations in tails when 2 or fewer observations (see Boos and Stephanski 2004 p. 22)
# goodness-of-fit statistics => by defn, the squared observed minus expected, divided by the expected
# ----------------------------------------------------------------------------------------------------
# 1st column : all the columns/cells
# 2nd column : collapse tail cells starting with the one whose value is 2 or fewer
#
g2.statistics <- matrix(NA, n.pops, 4) 
dimnames(g2.statistics) <- list(c(paste("pop",1:n.pops,sep=" ")),c("actual","collapsed", "no.cells", "p.val"))

# maximum number of ramets per genet in pooled tails data
max.pooled.ramet.count <- 5

data.collapsed.tails <- data.frame(pop=rep(1:n.pops,rep(max.pooled.ramet.count,n.pops)), ramet.count=rep(1:max.pooled.ramet.count,n.pops), obs.genets=(0), pred.genets=(0))

for (J in 1:n.pops) {
  m <- glm(obs.genets ~ ramet.count,family=poisson(link="log"),data=dat2,subset=pop==J)
  
  m.ols <- lm(log(obs.genets+1) ~ log(ramet.count),data=dat2,subset=pop==J)
  
  #m <- glm(obs.genets ~ ramet.count,family=quasipoisson(link="log"),data=dat2,subset=pop==J)
  
  simplerate.mat[J,] <- m$coefficients
  yhat.simple[J,] <- predict(m,newdata=NULL,type="response",se.fit=FALSE)
  yhat.ols[J,] <- predict(m.ols,newdata=NULL,type="response",se.fit=FALSE)
  
  observed.counts <- dat2$obs.genets[dat2$pop==J]
  predicted.values <- yhat.simple[J,]
  #cat("predicted.values = ", predicted.values,"\n")
  
  g2.poisson <- sum((observed.counts-predicted.values)^2/predicted.values)
  
  # find the ramet count categories that are >2
  tmp <- dat2[(dat2$pop==J) & (observed.counts>2),"ramet.count"]
  #print(tmp)
  #stop("testing\n")
  
  # store the maximum ramet count category that is >2
  last.cell <- max(tmp)
  
  # make index vector for obs and exp with pooled tail
  observed.cells <- 1:last.cell
  #cat("observed.cells indices = ", observed.cells,"\n")
  
  # make index vector
  collapsed.cells <- (last.cell+1):max.ramet.count
  #cat("collaposed.cells indices = ", collapsed.cells,"\n")
  
  observed.collapsed <- sum(observed.counts[collapsed.cells])

  predicted.collapsed <- sum(predicted.values[collapsed.cells])
  #cat("obs and exp collapsed cells = ",c(observed.collapsed,predicted.collapsed),"\n")
  
  #print(observed.cells)
  #print(observed.counts[observed.cells])
  
  # some steps to save observed and predicted values for later plotting
  # number of ramet categories after pooling tails
  length.collapsed <- length(observed.counts[observed.cells])
  

  # store observed ramets per genet counts for plotting later, except pooled
  data.collapsed.tails$obs.genets[((J*max.pooled.ramet.count)-(max.pooled.ramet.count-1)):((J*max.pooled.ramet.count)-(max.pooled.ramet.count-last.cell))] <- observed.counts[observed.cells]
  
  # store pooled ramets per genet counts for plotting later
  data.collapsed.tails$obs.genets[((((J-1)*max.pooled.ramet.count))+(last.cell + 1))] <- observed.collapsed
  
  
  # store predicted ramets per genet counts other than the pooled cell for plotting later
  data.collapsed.tails$pred.genets[((J*max.pooled.ramet.count)-(max.pooled.ramet.count-1)):((J*max.pooled.ramet.count)-(max.pooled.ramet.count-last.cell))] <- predicted.values[observed.cells]
  
  # store predicted ramets per genet counts from pooled cells for plotting later
  data.collapsed.tails$pred.genets[((((J-1)*max.pooled.ramet.count))+(last.cell + 1))] <- predicted.collapsed
  
  
  g2.poisson.approx <- 0
  # sum up over categories with observed ramets >2 first
  g2.possion.approx <- sum((observed.counts[observed.cells]-predicted.values[observed.cells])^2/predicted.values[observed.cells])
  # then add pooled tail categories of observed ramets
  g2.poisson.approx <- g2.poisson.approx + (observed.collapsed-predicted.collapsed)^2/predicted.collapsed
  
  #g2.ols <- sum((log(dat2$obs.genets[dat2$pop==J]+1)-yhat.ols[J,])^2/yhat.ols[J,])
  g2.statistics[J,1:2] <- c(g2.poisson,g2.poisson.approx)
  g2.statistics[J,3] <- last.cell + 1
  
  # degrees of freedom for conservative GOF test => no.of.cells minus 1
  g2.statistics[J,4] <- 1-pchisq(q=g2.poisson.approx,df=last.cell,ncp=0,lower.tail=TRUE,log.p=FALSE)
}

g2.statistics

# m.saturated => model for different rate and different intercept for each population
m.saturated <- glm(obs.genets ~ as.factor(pop)+ramet.count+as.factor(pop)*ramet.count,family=poisson(link="log"),data=dat2)
summary(m.saturated)

# m.reduced.1 => model for same rate parameter for all pops, but different intercepts for each pops
m.reduced.1 <- glm(obs.genets ~ as.factor(pop)+ramet.count,family=poisson(link="log"),data=dat2)
summary(m.reduced.1)

# m.reduced.2 => model for same rate parameter and same intercept parameter for all pops
m.reduced.2 <- glm(obs.genets ~ ramet.count,family=poisson(link="log"),data=dat2)
summary(m.reduced.2)


### ----------------------------------------------------------------------------------------------------------

# define plot window size
dev.new(width = 500, height = 300, unit = "px")

### Plot of the observed genet counts as a function of the ramet counts, plus a line through the fitted values
plot(dat2$ramet.count[dat2$pop==1],dat2$obs.genets[dat2$pop==1],type="n", ylim=c(0,70), xlim=c(1,17),
     xlab="Ramets per MLL", ylab="MLL count")
#, main="Distributions of ramets per genet\nin each population"

axis(1,at=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17))

for (J in 1:10) {
  points(dat2$ramet.count[dat2$pop==J],dat2$obs.genets[dat2$pop==J],pch=J,cex=0.5,col=J)
  
  lines(dat2$ramet.count[dat2$pop==J],yhat.simple[J,],col=J)
}

legend("topright", inset=.02, 
       legend=c("Pelham Bay (PB)", "Caumsett SP (CS)","Ambro Preserve (JA)","Wading River (WR)",
                "Sammys Beach (SB)","William T Davis (WT)","Alley Creek (AC)","Jones Beach (JB)", 
                "Cedar Beach (GI)","Cupsogue Beach (CB)"), 
       pch = c(1:10), lty= c(1:10), col = c(1:10), cex=0.9, pt.cex = 1
       )



# define plot window size
#dev.new(width = 500, height = 300, unit = "px")

### using ramets per genet with pooled tails
### Plot of observed genet counts as a function of the collapsed ramet counts, plus a line through the predicted values

plot(data.collapsed.tails$ramet.count[data.collapsed.tails$pop==1],data.collapsed.tails$obs.genets[data.collapsed.tails$pop==1],type="n", ylim=c(0,70), xlim=c(1,5), xaxt="n",
     xlab="Ramets per MLL with pooled tails for counts ≤2", ylab="MLL count")
#, main="Distributions of ramets per MLL genet\nin each population with pooled tails"
axis(1,at=c(1,2,3,4,5))

for (J in 1:10) {
  points(data.collapsed.tails$ramet.count[data.collapsed.tails$pop==J],data.collapsed.tails$obs.genets[data.collapsed.tails$pop==J],pch=J,cex=0.5,col=J)
  
  lines(data.collapsed.tails$ramet.count[data.collapsed.tails$pop==J],data.collapsed.tails$pred.genets[data.collapsed.tails$pop==J],col=J)
}

legend("topright", inset=.02, 
       legend=c("Pelham Bay (PB)", "Caumsett SP (CS)","Ambro Preserve (JA)","Wading River (WR)",
                "Sammys Beach (SB)","William T Davis (WT)","Alley Creek (AC)","Jones Beach (JB)", 
                "Cedar Beach (GI)","Cupsogue Beach (CB)"), 
       pch = c(1:10), lty= c(1:10), col = c(1:10), cex=0.9, pt.cex = 1
)

