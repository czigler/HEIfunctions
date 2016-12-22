plot.ps.spatmodel <-
function(psmod, plot = 1, group = NULL, contour=F){

require(MBA)

coords.1 <- psmod$m.1$coords
coords.0 <- psmod$m.0$coords
y.1 <- psmod$m.1$Y
y.0 <- psmod$m.0$Y
if(is.null(group) == T) group <- -1
if(group == -1) par(mfrow=c(1,2))

if(plot ==1){
if(group != 1){
obs.surf.0 <-mba.surf(cbind(coords.0, y.0), no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(obs.surf.0, xaxs = "r", yaxs = "r", main="Observed, No intervention", cex.main=0.9)
points(coords.0, pch=16, cex=0.5)
points(coords.1, pch=16, cex=0.5, col="blue")
if(contour==T){contour(obs.surf.0, add=T)}
}
if(group != 0){
obs.surf.1 <- mba.surf(cbind(coords.1, y.1), no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(obs.surf.1, xaxs = "r", yaxs = "r", main="Observed, intervention", cex.main=0.9)
points(coords.1, pch=16, cex=0.5)
points(coords.0, pch=16, cex=0.5, col="blue")
if(contour==T){contour(obs.surf.1, add=T)}
}
}

if(plot==2){
if(group != 1){
w.hat <- rowMeans(psmod$m.0$sp.effects)
w.surf <- mba.surf(cbind(coords.0, w.hat), no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(w.surf, xaxs = "r", yaxs = "r", main="Est. random effects, No Intervention", cex.main=0.9)
points(coords.0, pch=16, cex=.5)
points(coords.1, pch=16, cex=0.5, col="blue")
if(contour==T){contour(w.surf, add=T)}
}
if(group != 0){
w.hat <- rowMeans(psmod$m.1$sp.effects)
w.surf <- mba.surf(cbind(coords.1, w.hat), no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(w.surf, xaxs = "r", yaxs = "r", main="Est. random effects, Intervention", cex.main=0.9)
points(coords.1, pch=16, cex=.5)
points(coords.0, pch=16, cex=0.5, col="blue")
if(contour==T){contour(w.surf, add=T)}
}
}

if(plot==3){
if(group != 0){
ypred.0 <- psmod$ypred.0
obs.surf <- mba.surf(cbind(rbind(coords.0, coords.1), c(y.0, ypred.0)), no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(obs.surf, xaxs = "r", yaxs = "r", main="Predicted response, No intervention", cex.main=0.9)
points(coords.0, pch=16, cex=0.5)
points(coords.1, pch=16, cex=0.5, col="blue")
if(contour==T){contour(obs.surf, add=T)}
}
if(group != 1){
ypred.1 <- psmod$ypred.1
obs.surf <- mba.surf(cbind(rbind(coords.1, coords.0), c(y.1, ypred.1)), no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(obs.surf, xaxs = "r", yaxs = "r", main="Predicted response, Intervention", cex.main=0.9)
points(coords.1, pch=16, cex=0.5)
points(coords.0, pch=16, cex=0.5, col="blue")
if(contour==T){contour(obs.surf, add=T)}
}


}
}
