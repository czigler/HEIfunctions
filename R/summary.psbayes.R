
summary.psbayes <- function(out, full=T) {
	if(full==T){
		cat("Summary of Model 0 \n")
		print(summary(out$model.0))
		cat("\n \n Summary of Model 1 \n")
		print(summary(out$model.1))
		cat("\n")
	}
	cat("\n Causal Estimates \n \n")
	a <- data.frame(ate = mean(out$ate), sdate = sd(out$ate),
		 y.0 = mean(apply(out$ypred.0, 1, mean)), sdy.0 = sd(apply(out$ypred.0, 1, mean)), 
		 y.1 = mean(apply(out$ypred.1, 1, mean)), sdy.1 = sd(apply(out$ypred.1, 1, mean)))
	row.names(a) <- ""
	return(a)		
}
