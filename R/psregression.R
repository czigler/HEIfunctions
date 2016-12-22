


psregression <- function(formula, data, trt,nboot = NULL,...){ 
	
	require(splines)
	boot 	<- T
	if(is.null(nboot)==T){
		boot <- F
	}
	#remove missing y
	if(length(which(is.na(data[,all.vars(formula)[1]])==T)) > 0) {
		data <- data[-which(is.na(data[,all.vars(formula)[1]])==T),]
	}

	psinternal <- function(data) {
		a 	<- data[,trt]

		if(all.equal(unique(a), c(0,1)) == F)
			stop("a must be coded as values 0, 1.")
	
			model.1 	<- glm(formula, data = data[a==1,],...)
			mm.1 		<- t(model.matrix(formula, data=data[a==0,]))

			model.0 	<- glm(formula, data = data[a==0,],...)
			mm.0 		<- t(model.matrix(formula, data=data[a==1,]))

			if(is.null(model.1$offset) ==F){
				ypred.1 	<- c(model.1$y/model.1$family$linkinv(model.1$offset), 
						     model.1$family$linkinv(coef(model.1)%*%mm.1))
				ypred.0 	<- c(model.0$y/model.0$family$linkinv(model.0$offset), 
						     model.0$family$linkinv(coef(model.0)%*%mm.0))
			}
			if(is.null(model.1$offset) ==T){
				ypred.1 	<- c(model.1$y, 
						     model.1$family$linkinv(coef(model.1)%*%mm.1))
				ypred.0 	<- c(model.0$y, 
						     model.0$family$linkinv(coef(model.0)%*%mm.0))
			}
			
			ate 			<- mean(ypred.1 - ypred.0)
			internal 		<- list()
			internal$ypred.0 	<- ypred.0
			internal$ypred.1 	<- ypred.1
			internal$ate 	<- ate
			internal$model.0 <- model.0
			internal$model.1 <- model.1
			return(internal)
	}
	internal <- psinternal(data)

	if(boot == T){
		bootfunc <- function(nboot){
			sapply(1:nboot, function(i){
			n 	 <- dim(data)[1]
			msg <- "try-error"
			while(msg == "try-error"){
				newdat <- data[sample(1:n, n, replace=T),]
				out  	 <- try(psinternal(newdat), silent=T)
				msg 	 <- class(out)
			}
			return(out$ate)
			})
		}
		bootout <- bootfunc(nboot)
	}	

	out 			<- list()
	out$model.0 	<- internal$model.0
	out$model.1		<- internal$model.1
	out$ypred.0 	<- as.vector(internal$ypred.0)
	out$ypred.1 	<- as.vector(internal$ypred.1)
	out$ate 		<- internal$ate
	if(boot==T){
		out$boot <- bootout
	}
	class(out) 	<- "psregression"
	return(out)
}
