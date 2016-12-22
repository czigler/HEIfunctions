mergepredictions <-
function(dat, preddat, cov, predcov, id, predid = NULL, name = NULL) {
if(is.null(predid) == T) predid <- id
newdat <- merge(dat, preddat, by.x = id, by.y = predid, all=T)
col <- (ncol(newdat) + 1)
newdat[,col] <- newdat[,cov]
newdat[which(is.na(newdat[,cov]) == T),col] <- newdat[which(is.na(newdat[,cov]) == T),predcov]
if(is.null(name) ==T) colnames(newdat)[col] <- paste(cov, "new", sep="")
if(is.null(name) ==F) colnames(newdat)[col] <- name
return(newdat)
}
