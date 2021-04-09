

# Code to compute the host of meta-analytic results given in the main text and in the supplement.
#	The code then produces Figure S3 and Tables S1-S3.

# The Bayesian section implements the hierarchical model described in Gelman et al Bayesian Data Analysis, Chap 5 
#  		based on code from Appendix C in the same book
# 	Parts of Bayesian code here is then implemented in the script which generates Figures 4 and 5.


# Load some required libraries.  The user will need to install these if not already installed
rm(list=ls())
library(arm)
library(xtable)
library(Hmisc)
library(class)

# First replace the following directory with the directory where replication files were unzipped
setwd("/users/emily/Desktop/replication")

# Read in data
data=read.csv("data/standardized_effects.csv") 

# Capture any text output in a log file
sink("output/metaAnalysis_log.txt")


#-----------------------------------------------------------------------------
#------------ CLASSICAL APPROACH ---------------------------------------------
#-----------------------------------------------------------------------------
# Define a function that will estimate the median, the mean effect, the std error on the mean, and the precision-weighted distribution of a given sample of studies.  The function returns these things in a list, as well as the original effect sizes and their se, and the estimated density for plotting purposes

meta.anal <- function(toplot) {
  bt=as.numeric(as.character(data$effect_1sd))[toplot]*100
  sd=as.numeric(as.character(data$SE_1sd))[toplot]*100
  med <- median(bt) # Get median effect 
  ww=1/(sd)^2/sum(1/(sd)^2)  #inverse variance weights
  wt=matrix(ww,nrow=1)  #put weights in matrix
  beta=matrix(bt,ncol=1) #put estimates in matrix
  V=diag((sd)^2) #make variance-covariance matrix, setting covariances to zero
  mn<-wt%*%beta  #variance-weighted mean effect
  se<-sqrt(wt%*%V%*%t(wt))  #variance-weighted standard error
  
  # now generate precision-weighted distribution of effects, and get quantiles of this distribution
  x= seq(-50,100,0.1)
  mass <- matrix(nrow=length(bt),ncol=length(x))
  for (i in 1:length(bt)) {
    mass[i,]<-dnorm(x,bt[i],sd[i])
  }
  precision <- matrix(rep(ww,length(bt)),nrow=length(bt),byrow=T)
  dens <- apply(precision%*%mass,2,sum)
  cdf <- as.matrix(cumsum(dens)/sum(dens))  #get cdf
  pr <- as.matrix(c(0.025,0.050,0.250,0.500,0.750,0.950,0.975))  #quantiles we want
  qts = as.numeric(as.character(knn(cdf,pr,x,k=1)))
  zeroval <- cdf[x==0]  #how much probability mass below zero
  return(list(med,mn,se,qts,zeroval,bt,sd,x,dens,ww))  #return the stuff we want in a list
}


#-----------------------------------------------------------------------------
#------------ BAYESIAN APPROACH ---------------------------------------------
#-----------------------------------------------------------------------------
# Our approach here borrows code directly from the Gelman et al textbook, Bayesian Data Analysis, Appendix C. The reader is referred to Chapter 5 and Appendix C in that book for more information.

# First we define a function to compute the Bayesian statistics we want.  The function takes three arguments:  what sample we want, the max of the uniform prior on tau, and the number of simulations we want to run. 
# The function outputs a 5 element list, where the first element is the n x J matrix of estimated effect sizes (J is number of studies and n is number of simulation replications), the second element is the values of tau from each of the n runs, 3rd and 4th are the vector of original treatment effects and their se, and the 5th is the vector of estimated mu's

bayespost <- function(toplot, gridmax, n.sims) {
  
  y <- as.numeric(as.character(data$effect_1sd))[toplot]*100  
  sigma.y <- as.numeric(as.character(data$SE_1sd))[toplot]*100
  J=length(y)
  
  #define some other functions as in Gelman
  mu.hat <- function(tau, y, sigma.y) {
    sum(y/(sigma.y^2 + tau^2))/sum(1/(sigma.y^2 + tau^2))
  }
  V.mu <- function(tau, y, sigma.y) {
    1/sum(1/(tau^2 + sigma.y^2)) 
  }
  n.grid <- 2000
  tau.grid <- seq(0, gridmax, length=n.grid) #discrete distribution for tau
  log.p.tau <- rep(NA, n.grid)
  
  for (i in 1:n.grid){
    mu <- mu.hat(tau.grid[i], y, sigma.y) 
    V <- V.mu(tau.grid[i], y, sigma.y) 
    log.p.tau[i] <- 0.5*log(V) - 0.5*sum(log(sigma.y^2 + tau.grid[i]^2)) - 0.5*sum((y-mu)^2/(sigma.y^2 + tau.grid[i]^2))  # Equation 5.21 in Gelman et al. in log form, giving log density of tau
  }
  log.p.tau <- log.p.tau - max(log.p.tau) 
  p.tau <- exp(log.p.tau)
  p.tau <- p.tau/sum(p.tau)
  tau <- sample (tau.grid, n.sims, replace=T, prob=p.tau)  #sample n times from the discrete distribution, based on the probabilities just calculated
  
  mu <- rep (NA, n.sims)
  theta <- array (NA, c(n.sims,J)) 
  for (i in 1:n.sims){
    mu[i] <- rnorm (1, mu.hat(tau[i],y,sigma.y), sqrt(V.mu(tau[i],y,sigma.y)))
    if (tau[i]>0) {
      theta.mean <- (mu[i]/tau[i]^2 + y/sigma.y^2)/ (1/tau[i]^2 + 1/sigma.y^2)
      theta.sd <- sqrt(1/(1/tau[i]^2 + 1/sigma.y^2)) 
    } else {
      theta.mean <- mu[i]
      theta.sd <- 0	
    }
    theta[i,] <- rnorm (J, theta.mean, theta.sd) 
  }
  return(list(theta,tau,y,sigma.y, mu))
}


###############################################################
# CREATE FIGURE 4
###############################################################
set.seed(8675309)
toplot <- data$indiv_grou==2
numg=sum(toplot)  #number of group-level estimates
out <- meta.anal(toplot)

eff <- out[[6]] # Get point estimate
effsd <- out[[7]] # Get standard deviation 
hi=(eff+1.96*effsd)  #95% CI
lo=(eff-1.96*effsd)
med <- out[[1]] # Get median
m1 <- out[[2]] # Get mean
sd1 <- out[[3]] # Get standard deviation
ww <- out[[10]]  #inverse variance weights

#make color scheme:  
iv=data$Independent.variable[toplot]
colz=rep("red",length(iv)) # Temperature-related studies are in red
colz[iv=="rainfall"]="blue" # Rainfall-related studies are in blue
colz[iv=="rainfall dev."]="green" # Drought-related studies are in green
colz[iv=="ENSO"]="brown"
colz[iv=="storms"]="grey" # Storm-related studies are in grey
colz[iv%in%c("SPEI","PSDI")]="orange" # Other studies are in orange
obs=sum(data$obs[toplot])


pdf(file="output/Figure4.pdf",width=7.5,height=5,useDingbats=F)
par(mar=c(8,4,3,3))
ptsize <- 1.2
xl=c(0.2,numg+4.7)  #define xlim of full plot
yl=c(-8,13)  #define ylim of plot
plot(1,xlim=xl,ylim=yl,type="n",ylab=bquote(paste("% change per 1",sigma," change in climate",sep="")),xlab="",axes=F,xaxs="i",yaxs="i") # Plot the effect points
abline(v=numg+0.5,lwd=1,lty=1) # Line at mean 
axis(2,at=seq(-8,12,4),las=1) #  Add extra horizontal line
axis(4,at=seq(-8,12,4),las=1) # Add extra horizontal line
segments(0,seq(-20,20,4),numg+0.5,seq(-20,20,4),col="grey85",lty=1,lwd=0.5) # Add vertical bars
segments(0,0,numg+0.5,0,lty=1,lwd=1.5) # Add vertical bars
rect(0,m1-1.96*sd1,numg+0.5,m1+1.96*sd1,col="grey90",border=NA)  #CI on mean
segments(0,m1,numg+0.5,m1,lty=1,lwd=1.2)  #mean for group
segments(0,med,numg+0.5,med,lty=2)  #median for group
ll=1:sum(toplot)
segments(ll,lo,ll,hi,lwd=3.5,lty=1,col="black")
segments(ll,lo,ll,hi,lwd=3,lty=1,col=colz)
points(ll,eff,cex=ptsize,pch=21,bg=colz) # Point at effect size
mtext(paste(data$Study[toplot],data$study_year[toplot],sep=" "),side=1,at=ll-0.22,cex=0.5,las=2,line=0.5) # Add x labels with names of studies
mtext(data$subgroup[toplot],side=1,at=ll,cex=0.5,las=2,line=0.5)
mtext(data$Region[toplot],side=1,at=ll+0.22,cex=0.5,las=2,line=0.5)
arx=c(9.7,10.7)
ary=c(12,12)
arrows(arx,ary,arx,ary+0.7,length=0.05) # Add arrows for last two bars
text(arx,ary-0.5,c("16%","20%"),cex=0.5) # Add text for last two bars

#set up density plots
rect(numg+1,yl[1],numg+2.6,yl[2],col="grey93")
rect(numg+3,yl[1],numg+5,yl[2],col="grey93")

#density for all studies
rect(numg+1,-30,numg+2.6,60,col="grey93",border=NA) # Add sublot for density
x <- out[[8]]
dens <- out[[9]]
dens = dens*1.5/max(dens) #rescale to be height = 1.5
lines(dens+numg+1,x,lwd=2) # Add vertical lines
points(rep(numg+1.1,length(eff)),eff,pch="-",col="grey",cex=1.2) # Add point at mean
abline(v=numg+1,col="grey")

#bayes density
td<- bayespost(toplot,30,10000)
th <- td[[1]]
wd <- rep(ww,each=10000)/sum(rep(ww,10000)) # precision weights for each point
yy <- density(th,weights=wd)
ht=(yy$y - min(yy$y))/(max(yy$y)-min(yy$y))*1.5+numg+1
lines(ht,yy$x,lwd=2,lty=2)

#density for temperature only
toplot <- data$indiv_grou==2 & data$Independent.variable%in%c("temperature","SPEI","ENSO","PSDI")
out <- meta.anal(toplot)
eff <- out[[6]]
m2 <- out[[2]]
sd2 <- out[[3]]
ww <- out[[10]]  #inverse variance weights
x <- out[[8]]
dens <- out[[9]]
dens = dens*1.5/max(dens) #rescale to be height = 1.5
lines(dens+numg+3,x,col="red")  

td<- bayespost(toplot,30,10000)
th <- td[[1]]
wd <- rep(ww,each=10000)/sum(rep(ww,10000)) # precision weights for each point
yy <- density(th,weights=wd)
ht=(yy$y - min(yy$y))/(max(yy$y)-min(yy$y))*1.5+numg+3
lines(ht,yy$x,lwd=2,lty=2,col="red")
points(rep(numg+3.1,sum(toplot)),eff,pch="-",col="grey",cex=1.2)

# pooled point estimates
mns=c(m1,m2)
sds=c(sd1,sd2)
cilo=mns-1.96*sds
cihi=mns+1.96*sds
dd=c(1.3,3.3)
segments(numg+dd,cilo,numg+dd,cihi,col=c("black","red"),lwd=3)
points(numg+dd,mns,pch=21,bg=c("white","red"),cex=ptsize)
rect(numg+1,yl[1],numg+2.6,yl[2],lwd=1.5)
rect(numg+3,yl[1],numg+5,yl[2],lwd=1.5)
rect(xl[1],yl[1],numg+0.5,yl[2],lwd=1.5)
dev.off()



##################################################
# 	Create Figure 5
##################################################
set.seed(8675309)
toplot <- data$indiv_grou==1  #sample to plot
numg=sum(toplot)  #number of group-level estimates
out <- meta.anal(toplot)

eff <- out[[6]]
effsd <- out[[7]]
hi=(eff+1.96*effsd)  #95% CI
lo=(eff-1.96*effsd)
med <- out[[1]]
m1 <- out[[2]]
sd1 <- out[[3]]
ww <- out[[10]]  #inverse variance weights

#make color scheme:  
colz=rep("red",sum(toplot))
iv=data$Independent.variable[toplot]
colz[iv=="rainfall"]="blue"
colz[iv=="rainfall dev."]="green"
colz[iv=="ENSO"]="brown"
colz[iv=="storms"]="grey"
colz[iv%in%c("SPEI","PSDI")]="orange"
obs=sum(data$obs[toplot])

pdf(file="output/Figure5.pdf",width=12,height=5,useDingbats=F)
par(mar=c(8,4,3,3))
ptsize <- 1.3
xl=c(0.2,numg+4.7)  #define xlim of full plot
yl=c(-20,55)  #define ylim of plot
plot(1,xlim=xl,ylim=yl,type="n",ylab=bquote(paste("% change per 1",sigma," change in climate",sep="")),xlab="",axes=F,xaxs="i",yaxs="i")
abline(v=numg+0.5,lwd=1,lty=1)
axis(2,at=seq(-100,100,10),las=1)
axis(4,at=seq(-100,100,10),las=1)
segments(0,seq(-100,100,10),numg+0.5,seq(-100,100,10),col="grey85",lty=1,lwd=0.5)
segments(0,0,numg+0.5,0,lty=1,lwd=1.5)
rect(0,m1-1.96*sd1,numg+0.5,m1+1.96*sd1,col="grey90",border=NA)  #CI on mean
segments(0,m1,numg+0.5,m1,lty=1,lwd=1.2)  #mean for group
segments(0,med,numg+0.5,med,lty=2)  #median for group
ll=1:sum(toplot)
segments(ll,lo,ll,hi,lwd=3.5,lty=1,col="black")
segments(ll,lo,ll,hi,lwd=3,lty=1,col=colz)
points(ll,eff,cex=ptsize,pch=21,bg=colz[toplot])
mtext(paste(data$Study[toplot],data$study_year[toplot],sep=" "),side=1,at=ll-0.25,cex=0.5,las=2,line=0.5)
mtext(data$subgroup[toplot],side=1,at=ll,cex=0.5,las=2,line=0.5)
mtext(data$Region[toplot],side=1,at=ll+0.25,cex=0.5,las=2,line=0.5)
arx=c(19.7,20.7)
ary=c(yl[2]-6,yl[2]-2)
arrows(arx[1],ary[1],arx[1],ary[2],length=0.05)
arrows(arx[2],ary[1],arx[2],ary[2],length=0.05)
text(arx,ary[1]-1,c("71%","93%"),cex=0.5)

#set up density plots
rect(numg+1,yl[1],numg+2.6,yl[2],col="grey93")
rect(numg+3,yl[1],numg+5,yl[2],col="grey93")

#density plots for all studies
rect(numg+1,-30,numg+2.6,60,col="grey93",border=NA)
x <- out[[8]]
dens <- out[[9]]
dens = dens*1.5/max(dens) #rescale to be height = 1.5
lines(dens+numg+1,x,lwd=2)
points(rep(numg+1.1,length(eff)),eff,pch="-",col="grey",cex=1.2)
abline(v=numg+1,col="grey")

#bayes density
td<- bayespost(toplot,30,10000)
th <- td[[1]]
wd <- rep(ww,each=10000)/sum(rep(ww,10000)) # precision weights for each point
yy <- density(th,weights=wd)
ht=(yy$y - min(yy$y))/(max(yy$y)-min(yy$y))*1.5+numg+1
lines(ht,yy$x,lwd=2,lty=2)

#density for temperature only
toplot <- data$indiv_grou==1 & data$Independent.variable%in%c("temperature","SPEI","ENSO","PSDI")
out <- meta.anal(toplot)
eff <- out[[6]]
m2 <- out[[2]]
sd2 <- out[[3]]
ww <- out[[10]]  #inverse variance weights
x <- out[[8]]
dens <- out[[9]]
dens = dens*1.5/max(dens) #rescale to be height = 1.5
lines(dens+numg+3,x,col="red")  

td<- bayespost(toplot,30,10000)
th <- td[[1]]
wd <- rep(ww,each=10000)/sum(rep(ww,10000)) # precision weights for each point
yy <- density(th,weights=wd)
ht=(yy$y - min(yy$y))/(max(yy$y)-min(yy$y))*1.5+numg+3
lines(ht,yy$x,lwd=2,lty=2,col="red")
points(rep(numg+3.1,sum(toplot)),eff,pch="-",col="grey",cex=1.2)

# pooled point estimates
mns=c(m1,m2)
sds=c(sd1,sd2)
cilo=mns-1.96*sds
cihi=mns+1.96*sds
dd=c(1.3,3.3)
segments(numg+dd,cilo,numg+dd,cihi,col=c("black","red"),lwd=3)
points(numg+dd,mns,pch=21,bg=c("white","red"),cex=ptsize)
rect(numg+1,yl[1],numg+2.6,yl[2],lwd=1.5)
rect(numg+3,yl[1],numg+5,yl[2],lwd=1.5)
rect(xl[1],yl[1],numg+0.5,yl[2],lwd=1.5)
dev.off()



##################################################
# 	Create Table S1
##################################################

smpls <- list(data$indiv_grou==1,data$indiv_grou==1 & data$Independent.variable%in%c("temperature","SPEI","ENSO","PSDI"),data$indiv_grou==2,data$indiv_grou==2 & data$Independent.variable%in%c("temperature","SPEI","ENSO","PSDI"))
meds = mns = ses = qnt = c()
set.seed(8675309)
for (j in 1:4) {
  out <- meta.anal(smpls[[j]])
  meds = c(meds,out[[1]])
  mns = c(mns,out[[2]])
  ses = c(ses,out[[3]])
  qnt = rbind(qnt,out[[4]])
}
tbl<-data.frame(meds,mns,ses,qnt)
names(tbl)=c("Median","Mean","SE","2.5%","5%","25%","50%","75%","95%","97.5%")
rownames(tbl)=c("Intergroup","Intergroup (Temp.)","Interpersonal","Interpersonal (Temp)")
xtable(tbl,align=c("l",rep("c",10)))  #write out in latex form

# check probability mass below zero for both sets of studies
out <- meta.anal(smpls[[1]])
out[[5]]  # 10% of probability mass below zero for intergroup conflict
out <- meta.anal(smpls[[3]])
out[[5]]  # 0.04% of probability mass below zero



# ESTIMATES OF VARIANCE OF MEAN EFFECT UNDER DIFFERENT ASSUMPTIONS ABOUT COVARIANCES
# These results are reported in Section B of the SOM
toplot <- data$indiv_grou==1  #which estimates to focus on
numg=sum(toplot)  #number of group-level estimates
bt=as.numeric(as.character(data$effect_1sd))[toplot]*100
sd=as.numeric(as.character(data$SE_1sd))[toplot]*100
ww=1/(sd)^2/sum(1/(sd)^2)  #inverse variance weights
wt=matrix(ww,nrow=1)
beta=matrix(bt,ncol=1)
V=diag((sd)^2)
mn=wt%*%beta  #variance-weighted mean effect
sqrt(wt%*%V%*%t(wt))  #variance-weighted standard error
# what if all pairs of studies had betas correlated at r = 0.5. multiplying sigma_i*sigma_j*r = cov(i,j) because r = cov(i,j)/sigma_i*sigma_j
rr=c(0.1,0.3,0.5,0.7)
for (r in rr) {
  sd1=matrix(sd,nrow=1)
  V=t(sd1)%*%sd1
  for (i in 1:numg) {
    for (j in 1:numg)	{
      if (i!=j) {V[i,j]=V[i,j]*r}
    }
  }
  print(sqrt(wt%*%V%*%t(wt)))  #variance-weighted standard error
}


##################################################
# 	FIGURE S3, left panel
##################################################
set.seed(8675309)  #initialize seed since we are taking random draws and want to be able to replicate
toplot <- data$indiv_grou==1
mm <- bayespost(toplot,30,10000)
theta <- mm[[1]]
tau <- mm[[2]]
J = dim(theta)[2]
mu <- mm[[5]]

pdf(file="output/FigureS3a.pdf",width=5,height=5)
plot(1,xlim=range(tau),ylim=c(-20,50),type="n",ylab="Estimated treatment effect",xlab="tau",las=1,xaxs="i",yaxs="i",main="Intergroup conflict")
for (i in 1:J) {
  ll=lowess(tau[tau!=0],theta[tau!=0,i],f=0.4) 
  lines(ll$x,ll$y)
}
bk=seq(min(tau),max(tau),length=100)
th <-hist(tau,breaks=bk,plot=F)
yy=-20
for (i in 1:(length(bk)-1)) {  #manually draw in the density of tau
  rect(bk[i],yy,bk[i+1],yy+th$counts[i]/50,col="grey")
}
dev.off()


##################################################
#  TABLE S3. Uses parameters from what was just run to generate Figure S3a.
##################################################
qnt=c()
for (i in 1:J) {
  q=quantile(theta[,i],probs=c(0.025,0.25,0.5,0.75,0.975))
  qnt=rbind(qnt,q)
}
qnt=rbind(qnt,quantile(mu,probs=c(0.025,0.25,0.5,0.75,0.975)))
qnt=rbind(qnt,quantile(tau,probs=c(0.025,0.25,0.5,0.75,0.975)))
qnt=rbind(qnt,quantile(theta,probs=c(0.025,0.25,0.5,0.75,0.975)))
add=c("","","")
addnm=cbind(c("mu","tau","theta"),c("","",""))  #stuff we're adding
colnames(addnm)=names(data[toplot,1:2])  #so rbind will work in next line
out<- data.frame(rbind(data[toplot,1:2],addnm),as.numeric(c(mm[[3]],add)),as.numeric(c(mm[[4]],add)),qnt)
xtable(out) 


##################################################
# FIGURE S3, right panel
##################################################
toplot <- data$indiv_grou==2
mm <- bayespost(toplot,30,10000)
theta <- mm[[1]]
tau <- mm[[2]]
J = dim(theta)[2]
mu <- mm[[5]]

pdf(file="output/FigureS3b.pdf",width=5,height=5)
plot(1,xlim=c(0,max(tau)),ylim=c(-2,12),type="n",ylab="Estimated treatment effect",xlab="tau",las=1,yaxs="i",xaxs="i",main="Interpersonal violence")
for (i in 1:J) {
  ll=lowess(tau,theta[,i],f=0.5) 
  lines(ll$x,ll$y)
}
bk=seq(min(tau),max(tau),length=100)
th <-hist(tau,breaks=bk,plot=F)
yy=-2
for (i in 1:(length(bk)-1)) {  #manually draw in the density of tau
  rect(bk[i],yy,bk[i+1],yy+th$counts[i]/350,col="grey")
}
dev.off()


##################################################
# TABLE S2. uses parameters from what was just run to generate Figure S3b
##################################################
qnt=c()
for (i in 1:J) {
  q=quantile(theta[,i],probs=c(0.025,0.25,0.5,0.75,0.975))
  qnt=rbind(qnt,q)
}
qnt=rbind(qnt,quantile(mu,probs=c(0.025,0.25,0.5,0.75,0.975)))
qnt=rbind(qnt,quantile(tau,probs=c(0.025,0.25,0.5,0.75,0.975)))
qnt=rbind(qnt,quantile(theta,probs=c(0.025,0.25,0.5,0.75,0.975)))
add=c("","","")
addnm=cbind(c("mu","tau","theta"),c("","",""))  #stuff we're adding
colnames(addnm)=names(data[toplot,1:2])  #so rbind will work in next line
out<- data.frame(rbind(data[toplot,1:2],addnm),as.numeric(c(mm[[3]],add)),as.numeric(c(mm[[4]],add)),qnt)
rownames(out) <- NULL  #resets row numbers to start with 1
xtable(out) 

# end the log file
sink()


















