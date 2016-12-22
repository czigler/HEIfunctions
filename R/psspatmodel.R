


###############################################################
######### SPATIAL PROPENSITY SCORE REGRESSION MODEL  ##########
###############################################################


psspatmodel <- function(formula, data, trt, coords, na.rm=T, ...){

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

	coords.1 	<- subset(coordss,dats[,trt]==1)
	coords.0 	<- subset(coordss,dats[,trt]==0)

	if(na.rm==F){
		ap		<- data[,trt]
		yname 	<- names(frame)[1]	
		p.1 		<- (ap==1 | (is.na(data[,yname]) == T))
		p.0 		<- (ap==0 | (is.na(data[,yname]) == T))
		coords.1p 	<- coords[p.1,]
		coords.0p 	<- coords[p.0,]
		dat.0p	<- data[p.0,]
		dat.1p	<- data[p.1,]	
		dat.0p[all.vars(formula)[1]]	<- -40
		dat.1p[all.vars(formula)[1]]	<- -40
	}

	############ NONATTAINMENT AREAS ##############

	m.1 	<- spLM(formula, dat.1, coords=coords.1,...)

	if(na.rm == T)
		out.1 <- spPredict(m.1, pred.coords=coords.0, pred.covars= model.matrix(formula, dat.0))
	if(na.rm == F)
		out.1 <- spPredict(m.1, pred.coords=coords.0p, pred.covars= model.matrix(formula, dat.0p))
	
	ypred.1 	<- apply(out.1$p.predictive.samples, 1, mean)


	############ ATTAINMENT AREAS ##############

	m.0	<- spLM(formula, dat.0, coords=coords.0,...)

	if(na.rm == T)
		out.0 	<- spPredict(m.0, pred.coords=coords.1, pred.covars=model.matrix(formula, dat.1))
	if(na.rm == F){
		out.0 	<- spPredict(m.0, pred.coords=coords.1p, pred.covars=model.matrix(formula, dat.1p))
	}
	
	ypred.0 	<- apply(out.0$p.predictive.samples, 1, mean)

	############ ATE ##############
	#ate1 	<- mean(c(y.1, ypred.1)) - mean(c(y.0, ypred.0))

	nsamp	 	<- dim(out.0$p.predictive.samples)[2]
	atesamp 	<- sapply(1:nsamp, function(i) 
				mean(c(y.1, out.1$p.predictive.samples[,i]) - c(y.0, out.0$p.predictive.samples[,i])))

	out 		<- list()
	out$m.1 	<- m.1
	out$m.0 	<- m.0
	out$ate 	<- atesamp
	out$ypred.1 <- ypred.1
	out$ypred.0 <- ypred.0
	class(out) 	<- "ps.spatmodel"
	return(out)
}
