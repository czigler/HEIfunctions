plotigammaprior <-
function(vec, xlimit=c(0,100), nsamp = 10000, sd=T){
require(pscl)
if(sd==T){
samp <- sqrt(rigamma(nsamp,vec[1], vec[2]))
plot(density(samp), xlim=xlimit, main="Prior for SD")
}
if(sd==F){
samp <- rigamma(nsamp,vec[1], vec[2])
plot(density(samp), xlim=xlimit, main="Prior for Variance")
}
}
