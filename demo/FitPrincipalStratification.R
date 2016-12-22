#install.packages('/Users/coryzigler/Dropbox/Research/HEI Accountability Grant/PM10 Accountability/R Programs/HEIfunctions', repos=NULL, type="source")
#remove.packages('HEIfunctions')
# wd <- "/Users/coryzigler/Dropbox/Research/HEI Accountability Grant/PM10 Accountability/R Programs/"
# setwd(wd)
load(file="demo//pm10_withps.Rda")
library(splines)
library(HEIfunctions)
library(xtable)
library(coda)
library(corpcor)

###########################
##### FINALIZE DATA #######
###########################

### Prune observations that are out of range
dat=subset(dat, outofrange==0)

dat$y = log(dat$pm10fu)    #dat$fu - dat$base
dat$h = dat$Total_death.2001 #Respiratory.2001   ALL_CVD.2001
dat$denom = dat$Tot_Beneficiary.2001 #pyear.2001
dat$ytemp=dat$y
dat$ytemp[is.na(dat$ytemp)]= mean(dat$y, na.rm=TRUE)
mh0=glm(h~as.factor(pscat)  + pm10base1990 + ytemp + offset(log(denom)), family=poisson(link="log"), data=subset(dat, a==0))
mh1=glm(h~as.factor(pscat) + pm10base1990 + ytemp + offset(log(denom)), family=poisson(link="log"), data=subset(dat, a==1))


q=2
n=dim(dat)[[1]]
nburn= 100 #5000
nsamp= 105 #25000 ## nburnin + niter
thin=1
alpha0prop=summary(mh0)$cov.unscaled
alpha1prop=summary(mh1)$cov.unscaled
ypropsd=rep(.75, n*q)
ypropsd[seq(1,n*q,q)] = log(dat$pm10base1990)*.5  # for individual sampling
ypropsd[seq(2,n*q,q)] = log(dat$pm10base1990)*.5  # for individual sampling

coords   	<- cbind(dat$Longitude, dat$Latitude)
phistart 	<- makephi(coords, 10)
philower 	<- makephi(coords, 4) 
phiupper 	<- makephi(coords, 30)

tuning = list(A = 0.1, psi = 0.2, theta = 1, alpha0=.09*alpha0prop, alpha1=.09*alpha1prop, Y=ypropsd)
prior = list(KIG = rep(.5,2), psi = rep(.5,2), theta1 = rep(philower,2), theta2 = rep(phiupper,2))
starting = list(B=NULL, A = NULL, psi = NULL, theta = phistart, alpha0=mh0$coef, alpha1=mh1$coef)

starttime <- proc.time()
pstratamod <- principalstratmod(formula = y ~ as.factor(pscat) + pm10base1990, data = dat,
                                trt = "a", coords, 
                                formula_h = h ~ as.factor(pscat)  + pm10base1990, 
                                denom = "denom", 
                                nsamp, nburn, thin, tuning, prior, starting)
rtime=proc.time() - starttime
save(pstratamod, file = "pstratamod_smallmods_20150312.RData")
#load("/Users/coryzigler/Dropbox/Research/HEI Accountability Grant/PM10 Accountability/R Programs/results/20130411PrincipalStratification_medmods/pstratamod_medmods.RData")
#> rtime
#    user    system   elapsed 
#31376.418  1443.212 32598.076 

## Calculating causal effects
getrate=function(vec){
  return(vec*1000/dat$denom)
	}
y0out=pstratamod$y0[seq(1,(nsamp-nburn),10),]
y1out=pstratamod$y1[seq(1,(nsamp-nburn),10),]
h0out=t(apply(pstratamod$h0[seq(1,(nsamp-nburn),10),], 1, getrate)) #pstratamod$h0[seq(5001,20000,10),]#*1000/dat$Total_den
h1out=t(apply(pstratamod$h1[seq(1,(nsamp-nburn),10),], 1, getrate)) #pstratamod$h1[seq(5001,20000,10),]*1000/dat$Total_den

ceymat=y1out-y0out
atey=rowMeans(ceymat)
atty=rowMeans(ceymat[, dat$a==1])
mean(atey)
sd(atey)
mean(atty)
sd(atty)

cehmat=h1out-h0out
ateh=rowMeans(cehmat)
atth=rowMeans(cehmat[, dat$a==1])
mean(ateh)
sd(ateh)
mean(atth)
sd(atth)

Ceae = -5
Cede = 5

ede=rep(NA, dim(ceymat)[[1]])
eae=ede
for (i in 1:length(ede)){
	ede[i] = mean(cehmat[i,dat$a==1][abs(ceymat[i,dat$a==1])<=Cede])
	eae[i] = mean(cehmat[i,dat$a==1][ceymat[i,dat$a==1] < Ceae])
	}	
#time2=proc.time()-s

mean(ede)
sd(ede)
mean(eae)
sd(eae)


#### Check that Predicted = Observed
table(colMeans(h1out)[dat$a==1]==dat$Tot_ALL_cause_death[dat$a==1]*1000/dat$Total_den[dat$a==1])
table(colMeans(h0out)[dat$a==0]==dat$Tot_ALL_cause_death[dat$a==0]*1000/dat$Total_den[dat$a==0])

##### Check Range of Observed vs. Predicted Health Outcomes ######
par(mfrow=c(1,2))
plot(c(colMeans(h0out)[dat$a==1], dat$Tot_ALL_cause_death[dat$a==0]*1000/dat$Total_den[dat$a==0]), col=c(rep(2,sum(dat$a==1)), rep(1,sum(dat$a==0))), main="H0")
plot(c(colMeans(h1out)[dat$a==0], dat$Tot_ALL_cause_death[dat$a==1]*1000/dat$Total_den[dat$a==1]), col=c(rep(1,sum(dat$a==0)), rep(2,sum(dat$a==1))), main="H1")


#### Plot H0 vs. H1 for both groups #####
x11()
par(mfrow=c(1,2))
plot(colMeans(h0out)[dat$a==1], dat$Tot_ALL_cause_death[dat$a==1]*1000/dat$Total_den[dat$a==1], main="A=1 Areas", xlab="H(0) (predicted)", ylab="H(1) (observed)", xlim=range(c(colMeans(h0out), colMeans(h1out), dat$Tot_ALL_cause_death*1000/dat$Total_den)), ylim=range(c(colMeans(h0out), colMeans(h1out), dat$Tot_ALL_cause_death*1000/dat$Total_den)))
lines(c(-100000, 100000), c(-100000,100000))
plot(dat$Tot_ALL_cause_death[dat$a==0]*1000/dat$Total_den[dat$a==0], colMeans(h1out)[dat$a==0], main="A=0 Areas", xlab="H(0) (observed)", ylab="H(1) (predicted)", xlim=range(c(colMeans(h0out), colMeans(h1out), dat$Tot_ALL_cause_death*1000/dat$Total_den)), ylim=range(c(colMeans(h0out), colMeans(h1out), dat$Tot_ALL_cause_death*1000/dat$Total_den)))
lines(c(-100000, 100000), c(-100000,100000))

x11()
plot(colMeans(h0out), colMeans(h1out), col=dat$a+1, xlab="H(0)", ylab="H(1)")

x11()
plot(h0out[1,],h1out[1,], col=dat$a+1, xlab="H(0)", ylab="H(1)")


#Follow up pollution under observed nonattainment designations
dat$fu_obs=dat$y+dat$base
dat$fu_obs[is.na(dat$y) & dat$a==0]=colMeans(pstratamod$y0)[is.na(dat$y) & dat$a==0] + dat$base[is.na(dat$y) & dat$a==0] 
dat$fu_obs[is.na(dat$y) & dat$a==1]=colMeans(pstratamod$y1)[is.na(dat$y) & dat$a==1] + dat$base[is.na(dat$y) & dat$a==1] 

#Follow up pollution if nowhere nonattainment
dat$fu_allattain=rep(NA, n)
dat$fu_allattain[dat$a==0] = dat$fu_obs[dat$a==0]
dat$fu_allattain[dat$a==1] = colMeans(pstratamod$y0)[dat$a==1] + dat$base[dat$a==1]


meanbasea0=mean(dat$base[dat$a==0])
meanbasea1=mean(dat$base[dat$a==1])

meanfua1_obs=mean(dat$fu_obs[dat$a==1])
meanfua1_allattain=mean(dat$fu_allattain[dat$a==1])
meanfua0_obs=mean(dat$fu_obs[dat$a==0])

with(dat, plot(fu_obs, fu_allattain, col=a+1, pch=16, xlim=range(c(fu_obs, fu_allattain)),ylim=range(c(fu_obs, fu_allattain))   ))






x11()
#pdf(paste(picdir, '/pollutiondrops.pdf', sep=""))
plot(c(1,2.1), c(22, 41), type="n", axes=FALSE, xlab="", ylab="PM10")
lines(c(1,2), c(meanbasea0, meanfua0_obs), lwd=2)
lines(c(1,2), c(meanbasea1, meanfua1_obs), lwd=2, col=2)
lines(c(1,2), c(meanbasea1, meanfua1_allattain), lwd=2, col=2, lty=2)
arrows(2.01,24.3,2.01,25.5,code=3,col=2, lwd=2, length=.05)
text(2.07, 24.9, "Effect", col=2, font=2)
text(1.5,34.5,"n=157", col=2, font=2)
text(1.5,27.6,"n=192", font=2)
axis(side=1, at=c(1,2), labels=c("1987-1989", c("1999-2001")))
axis(side=2, at=seq(22,41,by=2), labels=seq(22,41,by=2))

#legend(1.3, 40, lty=c(1,1,2), lwd=c(2,2,2), col=c(1,2,2), c("Attainment", "Nonattainment (observed)", "Nonattanment (counterfactual)"))
#dev.off()

