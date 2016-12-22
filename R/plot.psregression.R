
plot.psregression <- function(psmodel, data, a, ps, cov, plot=1,...) {

	fitted.0 	<- as.numeric(psmodel$model.0$fitted.values)
	fitted.1	<- as.numeric(psmodel$model.1$fitted.values)
	names(fitted.0) <- "Fitted"

	rows.0 	<- names(psmodel$model.0$y)
	rows.1 	<- names(psmodel$model.1$y)

	a 		<- data[,a]
	
	ps.0		<- data[rows.0,ps]
	ps.1		<- data[rows.1,ps]
	names(ps.0) <- "Propensity score"
	
	if(plot==1){
		par(mar = c(1,1,1,1))
		range <- range(c(fitted.0, fitted.1))
		plot(ps.0, fitted.0, col=1, pch=16, ylim=range,xlim=c(0,1),...)
		points(ps.1,fitted.1, col=2, pch=16,...)
	}
	if(plot==2){
		y0 	<- psmodel$ypred.0
		y1 	<- psmodel$ypred.1
		range <- range(c(y0, y1))
		pspred.0 <- data[c(rows.0, rows.1),ps]
		pspred.1 <- data[c(rows.1, rows.0),ps]
		names(pspred.0) <- "Propensity score"
		names(y0) <- "Outcome"
		plot(pspred.0, y0, col=1, pch=16, ylim=range,...)
		points(pspred.1, y1, col=2, pch=16,...)
	}
	if(plot==3){
		covname <- cov
		cov.0 <- data[rows.0,cov]
		cov.1 <- data[rows.1,cov]
		plot(cov.0, fitted.0, col=1, pch=16, xlab=covname, ylab="Fitted",...)
		points(cov.1, fitted.1, col=2, pch=16,...)
	}
	if(plot==4){
		plot(psmodel$model.0$y, fitted.0, col=a+1, pch=16, xlab="Observed outcome", ylab="Fitted",...)
		points(psmodel$model.1$y, fitted.1, col=a+1, pch=16,...)
	}	
}


