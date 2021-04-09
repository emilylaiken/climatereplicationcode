
# R code to replicate Figures S1 and S2 in Hsiang, Burke, and Miguel 2013

# First replace the following directory with the directory where you unzipped
setwd("/Documents/Dropbox/Marshall-Sol/drafts/Science_review_v4/replication/")

# get median effect size
data=read.csv("data/standardized_effects.csv") 
eff=as.numeric(as.character(dta$effect_1sd[dta$indiv_grou==1]))
med <- median(eff)*100

###############################################
# 		Figure S1
###############################################

data=read.csv("data/Supplement_Data/FigS1_data.csv")  
toplot=1:dim(data)[1]  ##this is the plot order which we can change to whatever.
						# right now it's just in the order of the spreadsheet, which we sorted by hand
eff=as.numeric(as.character(data$effect_1sd))
effsd=as.numeric(as.character(data$SE_1sd))
hi=(eff+1.96*effsd)*100
lo=(eff-1.96*effsd)*100

pdf(file="output/FigureS1.pdf",width=8,height=5,useDingbats=F)
par(mar=c(8,4,3,1))
plot(1,xlim=c(1,dim(data)[1]),ylim=c(-15,46),type="n",ylab=bquote(paste("% change per 1",sigma," change in climate",sep="")),xlab="",xaxt="n",yaxt="n")
#rect(numg+0.5,-30,dim(data)[1]+2,60,col="grey93",border=NA)
axis(2,at=seq(-100,100,10),las=1)
abline(h=seq(-100,100,10),col="grey85",lty=1)
abline(h=0,lty=1,lwd=1.5)
ll=1:length(toplot)
wd=0.08  #width of end whiskers
segments(ll,lo[toplot],ll,hi[toplot],lwd=3,lty=1)
segments(ll,lo[toplot],ll,hi[toplot],lwd=1.5,lty=1,col="red")
segments(ll-wd,lo[toplot],ll+wd,lo[toplot],lwd=3)
segments(ll-wd,lo[toplot],ll+wd,lo[toplot],lwd=0.5,col="red")
segments(ll-wd,hi[toplot],ll+wd,hi[toplot],lwd=3)
segments(ll-wd,hi[toplot],ll+wd,hi[toplot],lwd=0.5,col="red")
points(ll,eff[toplot]*100,cex=1,pch=21,bg="red",lwd=2)
mtext(paste(data$Study[toplot],data$study_year[toplot],sep=" "),side=1,at=ll-0.15,cex=0.7,las=2,line=1)
mtext(data$notes[toplot],side=1,at=ll+0.15,cex=0.7,las=2,line=1)
abline(h=med,lty=2,lwd=1.5)
box()
dev.off()


###############################################
# 		Figure S2
###############################################

data=read.csv("data/Supplement_Data/FigS2_data.csv")  
toplot=1:dim(data)[1]  ##this is the plot order which we can change to whatever.
eff=as.numeric(as.character(data$temp_stdeff))
effsd=as.numeric(as.character(data$temp_stdeff_se))
hi=(eff+1.96*effsd)
lo=(eff-1.96*effsd)

# make shape of point
type=as.character(unique(data$climvar))
pchs=c(1:2)
ptype=pchs[match(data$climvar,type)]
lags=paste("lags = ",data$lags[toplot],sep="")
lags[5]=paste(lags[5],"**",sep="")

pdf(file="output/FigureS2.pdf",width=7,height=4,useDingbats=F)
par(mar=c(5,4,3,1))
plot(1,xlim=c(1,dim(data)[1]),ylim=c(-90,170),type="n",ylab=bquote(paste("% change per 1",sigma," change in climate",sep="")),xlab="",xaxt="n",yaxt="n")
rect(6.5,-100,17,200,col="grey93",border=NA)
abline(v=6.5)
abline(h=med,lty=2,lwd=1.5)  #median
axis(2,at=seq(-100,150,25),las=1)
abline(h=seq(-100,150,25),col="grey85",lty=1)
abline(h=0,lty=1,lwd=1.5)
ll=1:length(toplot)
wd=0.08  #width of end whiskers
segments(ll,lo[toplot],ll,hi[toplot],lwd=3,lty=1)
segments(ll,lo[toplot],ll,hi[toplot],lwd=1.5,lty=1,col="red")
segments(ll-wd,lo[toplot],ll+wd,lo[toplot],lwd=3)
segments(ll-wd,lo[toplot],ll+wd,lo[toplot],lwd=0.5,col="red")
segments(ll-wd,hi[toplot],ll+wd,hi[toplot],lwd=3)
segments(ll-wd,hi[toplot],ll+wd,hi[toplot],lwd=0.5,col="red")
points(ll,eff[toplot],cex=1,pch=ptype,bg="red",lwd=2)
mtext(lags,side=1,at=ll,cex=0.7,las=2,line=1)
mtext("Onset",side=3,at=3.5,cex=1.4,las=1,line=1)
mtext("Incidence",side=3,at=9.5,cex=1.4,las=1,line=1)
box()
dev.off()

