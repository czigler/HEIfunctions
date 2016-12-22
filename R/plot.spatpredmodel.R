plot.spatpredmodel <-
function(outmodel, origdat=NULL, preddat= NULL, covname= NULL, plotnum=1) {


#if(plotnum== 3 & outmodel$validate == FALSE)
#stop("Can only use plotnum 1 and 2 unless you are validating.")

if(plotnum== 2 & (is.null(origdat) == T | is.null(covname) ==T))
stop("Need to include origdat and covname for extrapolation plot (2).")

if(plotnum==1){
plot(outmodel$model$p.theta.samples)
}

if(plotnum==2){
require(MBA)
obs.surf <- mba.surf(cbind(outmodel$coords, origdat[,covname]), no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(obs.surf, xaxs = "r", yaxs = "r", main="Observed outcome, extrapolation")
points(outmodel$coords)
contour(obs.surf, add=T)
points(outmodel$predcoords, pch=16)
}
if(plotnum==3){
plot(preddat[,covname], outmodel$pred$predmean, pch=16, xlab="True pollution", ylab="Mean Predicted Pollution")
}
if(plotnum==4){
sdpred <- sapply(1:dim(preddat)[1], function(a) mean(sqrt((outmodel$pred$predmean[a] - origdat[a,covname])^2)))
plot(outmodel$pred$predsd, sdpred, xlab="SD of predicted pollution samples",
ylab="SD of predicted means from true pollution")
}
if(plotnum==5){
plot(rbind(outmodel$coords, outmodel$predcoords), xlab="Longitude", ylab="Latitude")
findextreme <- which((preddat[,covname]-outmodel$pred$predmean) == max((preddat[,covname]-outmodel$pred$predmean)))
points(outmodel$predcoords[c(findextreme, findextreme),], pch=16, cex=2)
}
}
