summary.spatpredmodel <-
function(outmodel, covname = NULL, origdat, preddat, full=T){

if(full==T){
cat("Summary of spLM Model \n")
print(summary(outmodel$model$p.theta.samples))
cat("\n")
}

cat("\nSummary of Prediction Model \n\n")
a <- data.frame(N = length(origdat[,covname]), Mean = mean(origdat[,covname]), SD = sd(origdat[,covname]))
a <- cbind(a, t(quantile(origdat[,covname])))
row.names(a) <- "original"

b <- data.frame(N = length(outmodel$pred$predmean), Mean = mean(outmodel$pred$predmean), SD = sd(outmodel$pred$predmean))
b <- cbind(b, t(quantile(outmodel$pred$predmean)))
row.names(b) <- "prediction"

c <- data.frame(N = length(outmodel$pred$predsd), Mean = mean(outmodel$pred$predsd), SD = sd(outmodel$pred$predsd))
c <- cbind(c, t(quantile(outmodel$pred$predsd)))
row.names(c) <- "SDprediction"

if(outmodel$validate == F){
out <- rbind(a, b, c)
}
if(outmodel$validate == T){
orig <- preddat[,covname]
pred <- outmodel$pred$predmean
d <- data.frame(N = length(orig), Mean = mean(orig-pred), SD = sd(orig-pred))
d <- cbind(d, t(quantile(orig-pred)))
row.names(d) <- "orig - pred"
out <- rbind(a, b, c, d)
}

cat("\n\n")

return(out)
}
