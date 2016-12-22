
###############################################################
######### SPATIAL PROPENSITY SCORE REGRESSION MODEL  ##########
###############################################################

psspatmodelglm <- function(formula, data, trt, coords, na.rm=T, ...){

	require(spBayes)

	row.names(coords) <- row.names(data) 
	############ MAKE DATA FOR MODEL FITTING ##############
	frame 	<- model.frame(formula, data)
	y 	  	<- model.extract(frame, "response")	
	dats 	  	<- data[names(y),]
	coordss  	<- coords[names(y),]

	dat.0 	<- subset(dats, dats[,trt]==0)
	dat.1 	<- subset(dats, dats[,trt]==1)

	frame.0 	<- model.frame(formula, dat.0)
	frame.1 	<- model.frame(formula, dat.1)	

	y.0 	  	<- model.extract(frame.0, "response")
	y.1 	  	<- model.extract(frame.1, "response")
	#y.0 	  	<- model.extract(frame.0, "response")/exp(model.offset(frame.0))
	#y.1 	  	<- model.extract(frame.1, "response")/exp(model.offset(frame.1))

	coords.1 	<- subset(coordss,dats[,trt]==1)
	coords.0 	<- subset(coordss,dats[,trt]==0)

	#phiest.1 	<- -log(0.05)/(max(dist(coords.1))/5)
	#phiest.0 	<- -log(0.05)/(max(dist(coords.0))/5)
	#if(is.null(phiprior) == T){
	#	phiprior.0 <- c(phiest.0 - phiest.0/1000, phiest.0 + phiest.0/1000)
	#	phiprior.1 <- c(phiest.1 - phiest.1/1000, phiest.1 + phiest.1/1000)
	#}
	#if(is.null(phiprior)==F){
	#	phiprior.0 <- phiprior
	#	phiprior.1 <- phiprior
	#}

	if(na.rm==F){
		ap		<- data[,trt]
		yname 	<- names(frame)[1]	
		p.1 		<- (ap==1 | (is.na(data[,yname]) == T))
		p.0 		<- (ap==0 | (is.na(data[,yname]) == T))
		coords.1p 	<- coords[p.1,]
		coords.0p 	<- coords[p.0,]
		dat.0p	<- data[p.0,]
		dat.1p	<- data[p.1,]	
	}

	
	############ NONATTAINMENT AREAS ##############

	m.1 	<- spGLM(formula, data=dat.1, coords=coords.1,...)

	if(na.rm == T){
		out.1 <- spPredict(m.1, pred.coords=coords.0, pred.covars= model.matrix(formula, dat.0))
		#off.1 <- model.offset(model.frame(formula, dat.0))

	}
	if(na.rm == F){
		out.1 <- spPredict(m.1, pred.coords=coords.0p, pred.covars= model.matrix(formula, dat.0p))
		#off.1 <- model.offset(model.frame(formula, dat.0p))
	}
	ypred.1 	<- apply(out.1$y.pred, 1, mean)

	############ ATTAINMENT AREAS ##############


	m.0 	<- spGLM(formula, data=dat.0, coords=coords.0,...)

	if(na.rm == T){
		out.0 	<- spPredict(m.0, pred.coords=coords.1, pred.covars=model.matrix(formula, dat.1))
		#off.0 	<- model.offset(model.frame(formula, dat.1))
	}
	if(na.rm == F){
		out.0 	<- spPredict(m.0, pred.coords=coords.1p, pred.covars=model.matrix(formula, dat.1p))
		#off.0 	<- model.offset(model.frame(formula, dat.1p))
	}
	
	ypred.0 	<- apply(out.0$y.pred, 1, mean)

	############ ATE ##############

	nsamp	 	<- dim(out.0$y.pred)[2]
	ysamp.0 	<- sapply(1:nsamp, function(i) mean(c(y.0, out.0$y.pred[,i])))
	ysamp.1 	<- sapply(1:nsamp, function(i) mean(c(y.1, out.1$y.pred[,i])))

	atesamp 	<- ysamp.1 - ysamp.0

	out 		<- list()
	out$m.1 	<- m.1
	out$m.0 	<- m.0
	out$ate 	<- atesamp
	out$ypred.1 <- ypred.1
	out$ypred.0 <- ypred.0
	out$y.0	<- y.0
	out$y.1 	<- y.1
	#out$off.0	<- off.0	
	#out$off.1	<- off.1
	class(out) 	<- "ps.spatmodel"
	return(out)
}


