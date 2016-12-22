

psbayes <- function(formula, data, trt, family = "gaussian", approx=F, ...){
	
	#remove missing y
	if(length(which(is.na(data[,all.vars(formula)[1]])==T)) > 0) {
		data <- data[-which(is.na(data[,all.vars(formula)[1]])==T),]
	}

	a  	<- data[,trt]
	n.0 	<- length(which(a == 0))
	n.1 	<- length(which(a == 1))
	n 	<- n.0 + n.1

	if(approx==F){
		require(MCMCpack)
		if(family == "gaussian"){
			model.0 	<- MCMCregress(formula, data=data[a==0,],...)
			mm.0 		<- t(model.matrix(formula, data=data[a==1,]))
			ypred.0 	<- model.0[,1:(dim(mm.0)[1])]%*%mm.0
			y.0 		<- as.numeric(attr(model.0, "y"))
	
			model.1 	<- MCMCregress(formula, data=data[a==1,],...)
			mm.1 		<- t(model.matrix(formula, data=data[a==0,]))
			ypred.1	<- model.1[,1:(dim(mm.1)[1])]%*%mm.1
			y.1 		<- as.numeric(attr(model.1, "y"))
		}
		if(family == "poisson"){
			stop("Not currently supported.")
			#model.0 	<- MCMCpoisson(formula, data=data[a==0,],...) # doesn't allow offsets!
			#mm.0 		<- t(model.matrix(formula, data=data[a==1,]))
			#frame.0 	<- model.frame(formula, data=data[a==1,])
			#offset.0 	<- model.offset(frame.0)
			#ypred.0 	<- exp(model.0[,1:(dim(mm.0)[1])]%*%mm.0)
			#model.1 	<- MCMCpoisson(formula, data=data[a==1,],...)
			#mm.1 		<- t(model.matrix(formula, data=data[a==0,]))
			#frame.1 	<- model.frame(formula, data=data[a==0,])
			#offset.1 	<- model.offset(frame.1)
			#ypred.1	<- exp(model.1[,1:(dim(mm.1)[1])]%*%mm.1)
		}
	}
	if(approx==T){
		require(arm)
		model.0 	<- bayesglm(formula, data=data[a==0,], family = family, ...)
		mm.0 		<- t(model.matrix(formula, data=data[a==1,]))
		if(family == "gaussian") ypred.0 	<- coef(sim(model.0,...))%*%mm.0
		if(family == "poisson") {	
			ypred.0 	<- exp(coef(sim(model.0,...))%*%mm.0)
			off.0 	<- as.numeric(model.0$offset)
		}
		y.0 		<- as.numeric(model.0$y)

		model.1 	<- bayesglm(formula, data=data[a==1,], family = family, ...)
		mm.1 		<- t(model.matrix(formula, data=data[a==0,]))
		ypred.1 	<- coef(sim(model.1,...))%*%mm.1
		if(family == "gaussian") ypred.1 	<- coef(sim(model.1,...))%*%mm.1
		if(family == "poisson") {
			ypred.1	<- exp(coef(sim(model.1,...))%*%mm.1)
			off.1 	<- as.numeric(model.1$offset)
		}
		y.1 		<- as.numeric(model.1$y)

	}


	out 		<- list()
	out$model.0 <- model.0
	out$model.1 <- model.1
	out$ypred.0 <- ypred.0
	out$ypred.1 <- ypred.1
	out$y.0 	<- y.0
	out$y.1 	<- y.1
	if(family=="gaussian")
		out$ate 	<- (n.1*mean(y.1) + n.0*apply(ypred.1,1, mean))/n - 
			   	   (n.0*mean(y.0) + n.1*apply(ypred.0, 1, mean))/n
	if(family=="poisson") {
		out$off.1 	<- off.1
		out$off.0	<- off.0
		out$ate 	<- (n.1*mean(y.1/exp(off.1)) + n.0*apply(ypred.1,1, mean))/n - 
			   	   (n.0*mean(y.0/exp(off.0)) + n.1*apply(ypred.0, 1, mean))/n
	}
	out$family 	<- family
	class(out)  <- "psbayes"
	return(out)
}
