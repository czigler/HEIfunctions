
summary.psregression <- function(out, full=T) {
	if(full==T){
		cat("Summary of Model 0 \n")
		print(summary(out$model.0))
		cat("\n \n Summary of Model 1 \n")
		print(summary(out$model.1))
		cat("\n")
	}
	cat("\n Causal Estimates \n \n")
	a <- data.frame(ate = out$ate, sdate = sd(out$boot), y.0 = mean(out$ypred.0), y.1 = mean(out$ypred.1))
	row.names(a) <- ""
	return(a)		
}
