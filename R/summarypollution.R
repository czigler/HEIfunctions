summarypollution <-
function(base, fu, a, data, na.rm=T){
meanfunc <- function(a) mean(a, na.rm=na.rm)
sdfunc<- function(a) sd(a, na.rm=na.rm)

basename <- base
funame <- fu
aname <- a
base <- data[,base]
fu   <- data[,fu]
a    <- data[,a]

out <- rbind( tapply(base, a, mean, na.rm=na.rm),
tapply(base, a, sd, na.rm=na.rm),
tapply(fu, a, mean, na.rm=na.rm),
tapply(fu, a, sd, na.rm=na.rm),
tapply(fu-base, a, mean, na.rm=na.rm),
tapply(fu-base, a, sd, na.rm=na.rm))
rownames(out) <- c(paste(basename, "mean"), paste(basename, "sd"), 
 paste(funame, "mean"),   paste(funame, "sd"), 
"diffmean", "diffsd")
return(out)
}
