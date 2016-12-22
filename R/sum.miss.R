sum.miss <-
function(whichvars, dat) {
missind <- rep(0, nrow(dat))
out <-  data.frame(nmiss=0, ntotal=0, pctmiss=0)
for(i in 1:length(whichvars)){
cov <- dat[,names(dat) %in% whichvars[i]]
out[i,] <- c(a <- length(bb <- which(is.na(cov)==T)), length(cov), 100*a/length(cov))
missind[bb] <- 1
}
row.names(out) <- whichvars
out[(length(whichvars) + 1),] <- c(sum(missind), nrow(dat), sum(missind)/nrow(dat))
row.names(out)[length(whichvars) + 1] <- "TotalMissing"
return(out) 
}
