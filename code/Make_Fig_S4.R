
# R script to make Figure S4 in Hsiang, Burke, Miguel 2013

# First replace the following directory with the directory where replication files were unzipped
setwd("/Documents/Dropbox/Marshall-Sol/drafts/Science_review_v4/replication/")

# Read in data
data=read.csv("data/standardized_effects.csv") 

# Make figure
iv=data$Independent.variable
eff=as.numeric(as.character(data$their_effect))*100
effsd=as.numeric(as.character(data$their_se))*100
dof <- as.numeric(as.character(data$their_dof))
tstat=log(abs(eff)/effsd)
rdf=log(sqrt(dof))
toplot=1:dim(data)[1]  
zz=lm(tstat[toplot] ~ rdf[toplot])
summary(zz)
pchs=rep(19,dim(data)[1])
pchs[iv%in%c("temperature","SPEI","ENSO","PSDI")]=24
pdf(file="output/FigureS4.pdf",width=6,height=6)
plot(rdf,tstat,xlim=c(0,8),ylim=c(0,8),las=1,xlab="log sqrt degrees of freedom",ylab="log t-stat",pch=pchs)
abline(a=zz$coefficients[1],b=1)  #set intercept of 45 degree line to intercept of OLS all
abline(a=zz$coefficients[1],b=zz$coefficients[2],lty=2)
abline(h=log(2))
toplot=iv%in%c("temperature","SPEI","ENSO","PSDI")  
tt=lm(tstat[toplot] ~ rdf[toplot])
abline(a=tt$coefficients[1],b=tt$coefficients[2],lty=2)
dev.off()

