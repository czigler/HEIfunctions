

principalstratmod_monoind <- function(formula, data, trt, formula_h, denom, nsamp, nburn, thin,
		tuning = list(psi = 0.2, alpha0=alpha0prop, alpha1=alpha1prop, B0=B0prop, B1=B1prop, Y=ypropsd), 
		prior = list(psi = rep(.01,2)),
		starting = list(B=NULL, psi = NULL, alpha0=NULL, alpha1=NULL)){

	require(MCMCpack)
	require(corpcor)
	require(MASS)
	require(msm)

	ncov <- dim(model.matrix(formula, data))[2]

	print("Note: function does not accept missing data in covariates.")
	logit	= function(theta, a, b){
		return(log((theta-a)/(b-theta)))
	}


	logitInv	= function(z, a, b){
	  return ( b-(b-a)/(1+exp(z)) );
	}

	###########################
	##### Format the Data #####
	###########################

	index <- seq(from = (nburn+1), to = nsamp, by = thin)
	n 	<- dim(data)[1]
	q	<- 2
	
	#pollution model
	names  <- colnames(model.matrix(formula, data))
	data$a <- data[,trt]
	data$y <- data[,all.vars(formula)[1]]
	data$ytemp 	<- 1
	nterm 	<- length(labels(terms(formula)))
	terms 	<- NULL
	for(ii in 1:nterm){
		if(ii == 1) terms <- labels(terms(formula))[ii]
		if(ii > 1) terms 	<- paste(terms, "+", labels(terms(formula))[ii], sep=" ")
	}
	formulatemp <- as.formula(paste("ytemp ~ ", terms))
	covs 	 	<- model.matrix(formulatemp, data)
	
	#health model
	names_h  <- colnames(model.matrix(formula_h, data))
	data$h <- data[,all.vars(formula_h)[1]]
	data$ytemp 	<- 1
	nterm 	<- length(labels(terms(formula_h)))
	terms 	<- NULL
	for(ii in 1:nterm){
		if(ii == 1) terms <- labels(terms(formula_h))[ii]
		if(ii > 1) terms 	<- paste(terms, "+", labels(terms(formula_h))[ii], sep=" ")
	}
	formulatemp <- as.formula(paste("ytemp ~ ", terms))
	covs_h 	 	<- model.matrix(formulatemp, data)
	
	makeXmat	= function(Xvals){
		#Xvals	= cbind(rep(1,dim(Xvals)[[1]]), Xvals)
		p	= dim(Xvals)[[2]]
		Xout	= matrix(0, n*q, p*q)
		Xout[seq(1,n*q,q),1:p]		= Xvals
		Xout[seq(2,n*q,q),(p+1):(2*p)]= Xvals
		return(Xout)
	}

	X	  = makeXmat(covs)
	p	  = ncol(X)
	
	## Design matrix for the health model WIITHOUT the posttreatment y
	X_h   = covs_h
	p_h = ncol(X_h) + 1 ## the +1 is to account for the posttreatment y to be added later

	y0			= rep(NA,n)
	y1			= rep(NA,n)
 	y0[data$a==0]	= data$y[data$a==0]
	y1[data$a==1]	= data$y[data$a==1]

	
	Y			= rep(-999,n*q)
	Y[seq(1,n*q,q)]	= y0
	Y[seq(2,n*q,q)]	= y1
	ismissy		= (is.na(Y))
	nwithmissing	= sum(ismissy)
	
	H			= rep(-999,n*2)
	H[seq(1,n*2,2)][data$a==0]	= data$h[data$a==0]
	H[seq(2,n*2,2)][data$a==1]	= data$h[data$a==1]
	
	denom=data[,denom]
	
	a			= rep(NA, n*q)
	a[seq(1,n*q,q)]	= data$a
	a[seq(2,n*q,q)]	= data$a

	

	###########################
	##### Select Priors #######
	###########################
	
	####Priors for Psi (Inverse gamma)
	psiig_a	= prior$psi
	psiig_b	= prior$psi
	
	
	##################################
	##### Select Starting Values #####
	##################################

	my0	= summary(lm(formula, data=subset(data, a[seq(1,n*q,q)]==0)))
	my1	= summary(lm(formula, data=subset(data, a[seq(1,n*q,q)]==1)))
	
	missy1	= ismissy[seq(1,n*q,q)]
	missy2	= ismissy[seq(2,n*q,q)]
	
	# simulate one starting value for locations with both Y missing
	missboth = rep(FALSE, n)
	missboth[is.na(Y[seq(1,n*q,q)]) & is.na(Y[seq(2,n*q,q)])] = TRUE
	startm0 = X[seq(1,n*q,q),1:(p/2)]%*%my0$coef
	startm1 = X[seq(2,n*q,q),((p/2)+1):p]%*%my1$coef
	Y[seq(1,n*q,q)][missboth==TRUE] = rnorm(sum(missboth==TRUE) ,startm0[missboth==TRUE], my0$sigma)
	
	# simulate starting values subject to monotonicity
	Y[seq(2,n*q,q)][missy2==1]	= rtnorm(sum(missy2==1), startm1[missy2==1], my1$sigma, lower = -Inf, upper = Y[seq(1,n*q,q)][missy2==1])	
	Y[seq(1,n*q,q)][missy1==1]	= rtnorm(sum(missy1==1), startm0[missy1==1], my0$sigma, lower = Y[seq(2,n*q,q)][missy1==1], upper = Inf)	
					   							   
	H[seq(1,n*2,2)][data$a==1]=rpois(sum(data$a==1), 10)
	H[seq(2,n*2,2)][data$a==0]=rpois(sum(data$a==0), 10)

	#Regression Coefficients
	if(is.null(starting$B) == F){
		B	= starting$B
	}
	if(is.null(starting$B) == T){
		B 	= rep(0, ncov*2)
	}


	#Vector of log of residual variances (Psi)
	if(is.null(starting$psi) == T){
		initpsis	= log(c(my0$sigma^2, my1$sigma^2))
	}else{
		initpsis 	= starting$psi
	}
	
	#Health Model Coefficients
	if(is.null(starting$alpha0) == T){
		initalpha	= rep(0, 2*p_h)
	}else{
		initalpha 	= c(starting$alpha0, starting$alpha1)
	}


	##################################
	##### Select Tuning Parameters ###
	##################################

	####Proposal SDS for log(Psi)
	psipropsds	= tuning$psi
	
	###Propoisal covariance matrices for alpha from poisson models
	alpha0propcov   = tuning$alpha0
	alpha1propcov   = tuning$alpha1
	
	###Propoisal covariance matrices for B
	B0propcov   = tuning$B0
	B1propcov   = tuning$B1
	Bpropcov = list()
	Bpropcov[[1]] = B0propcov
	Bpropcov[[2]] = B1propcov
	
	###Proposal SDS for missing Y
	propysds = tuning$Y
	
	#################################
	##### Define Some Functions #####
	#################################

loglike_norm_monoind	= function(Yvals,Bvals,psivals, lowervals, uppervals){
		#REMEMBER: psivals is on the log scale!!
				
		m = as.vector(X%*%Bvals)
		index = seq(1,n*q,q)
		temp = dtnorm(Yvals[index], mean=m[index], sd=sqrt(exp(psivals[1])), lower=lowervals[index], upper=uppervals[index], log=TRUE)
		bigindex = which(temp==Inf)
		temp[bigindex] = log( dnorm(Yvals[index][bigindex], mean=m[index][bigindex], sd=sqrt(exp(psivals[1]))) / (1 - 0.9999999999999999) )
		llike = sum(temp)
		
		index = seq(2,n*q,q)
		temp = dtnorm(Yvals[index], mean=m[index], sd=sqrt(exp(psivals[2])), lower=lowervals[index], upper=uppervals[index], log=TRUE)
		bigindex = which(temp==Inf)
		temp[bigindex] = log( dnorm(Yvals[index][bigindex], mean=m[index][bigindex], sd=sqrt(exp(psivals[1]))) / (1 - 0.9999999999999999) )
		llike = llike + sum(temp)
		
		## q inverse gamma priors
		###Jacobian for log transformation: \sum log(sigmasq)
		llike = llike + -(psiig_a[1] +1)*psivals[1] - psiig_b[1]/exp(psivals[1]) + psivals[1]
		llike = llike + -(psiig_a[2] +1)*psivals[2] - psiig_b[2]/exp(psivals[2]) + psivals[2]
		
		return(llike)
	}
	
loglike_pois     = function(Yvals,alphavals,h0vals,h1vals){
		##Poisson Part
		ya0  = Yvals[seq(1,n*q,q)]
		ya1  = Yvals[seq(2,n*q,q)]
		Xmat0  = as.matrix(cbind(X_h, ya0))
		alpha0  = alphavals[1:p_h]
		lambda0  = exp(Xmat0%*%alpha0 + log(denom))
		
		Xmat1  = as.matrix(cbind(X_h, ya1))
		alpha1  = alphavals[(p_h+1):length(alphavals)]
		lambda1  = exp(Xmat1%*%alpha1 + log(denom))

		llike  = sum( h0vals*log(lambda0) - lambda0 - lfactorial(h0vals) )
		llike  = llike + sum( h1vals*log(lambda1) - lambda1 - lfactorial(h1vals) )

		return(llike)
	}	


MHstep_pois  = function(whichalphas, Yvals, alphavals, h0vals, h1vals, currentll){
	accept=0
	llret=currentll
	alpharet=alphavals
	alphaprop=alphavals
	index=1:p_h
	if (whichalphas==1){
		index=(p_h+1):length(alphavals)
		}
	
	alphaprop[index]=mvrnorm(1, alpharet[index], alphapropcovs[[whichalphas+1]])
	llprop=loglike_pois(Yvals, alphaprop, h0vals, h1vals)
	ratio=exp(llprop-currentll)
	ratio[ratio>1]=1
	if (runif(1)<=ratio){
		alpharet=alphaprop
		llret=llprop
		accept=1
		}
		
	outlist=list(alpharet, llret, accept)	
	names(outlist)=c("alpha", "ll","accepted")
	return(outlist)
	}



MHstep_B = function(whichb, Yvals, Bvals, psivals, lowervals, uppervals, currentll){
		accept=0
		llret=currentll
		Bret=Bvals
		Bprop=Bvals
		index = 1:(p/2)
		if (whichb==1)
			index = ((p/2)+1):p
		
		Bprop[index] = mvrnorm(1, Bret[index], Bpropcov[[whichb + 1]])
		llprop = loglike_norm_monoind(Yvals, Bprop, psivals, lowervals, uppervals)
		ratio = exp(llprop - currentll)
		if (runif(1)<=ratio){
		Bret=Bprop
		llret=llprop
		accept=1
		}
		
	outlist=list(Bret, llret, accept)	
	names(outlist)=c("B", "ll","accepted")
	return(outlist)
	}

MHstep_psi = function(whichpsi,Yvals, Bvals, psivals, lowervals, uppervals, currentll){
		accept=0
		llret=currentll
		psiret=psivals
		psiprop=psivals
		
		psiprop[whichpsi+1] = rnorm(1, psiret[whichpsi+1], psipropsds[whichpsi+1])
		llprop = loglike_norm_monoind(Yvals, Bvals, psiprop, lowervals, uppervals)
		ratio = exp(llprop - currentll)
		if (runif(1)<=ratio){
		psiret=psiprop
		llret=llprop
		accept=1
		}
		
	outlist=list(psiret, llret, accept)	
	names(outlist)=c("psi", "ll","accepted")
	return(outlist)
	}


sampleH=function(whicha, alphavals, Yvals){
		ya0  = Yvals[seq(1,n*q,q)]
		ya1  = Yvals[seq(2,n*q,q)]
		index  =  1:p_h
		Xmat  =  cbind(X_h, ya0)
		
		if (whicha==1){
			index = (p_h+1) : length(alphavals)
			Xmat  = cbind(X_h, ya1)
		}
	
	
	alphs=alphavals[index]
	lamda=exp(Xmat%*%alphs + log(denom))
		
	h=rpois(sum(data$a==(1-whicha)), lamda[data$a==(1-whicha)])
	return(h)
	}
	
sampleY_monoind=function(ids, Yvals, Bvals, psivals, alphavals, h0vals, h1vals, lowervals, uppervals, currentllpois, currentllnorm){
		accept = rep(0, length(which(ids)))
		currentll = currentllnorm + currentllpois
		llretnorm = currentllnorm
		llretpois = currentllpois
		Yret = Yvals
		Yprop = Yvals
				
		Yprop[ids] = rtnorm( length(which(ids)), Yvals[ids], propysds[ids], lower=lowervals[ids], upper=uppervals[ids])

		llpropnorm = loglike_norm_monoind(Yprop, Bvals, psivals, lowervals, uppervals)
		llproppois = loglike_pois(Yprop, alphavals, h0vals, h1vals)
		
		propdens_prop = dtnorm(Yprop[ids], Yvals[ids], propysds[ids], lower=lowervals[ids], upper=uppervals[ids], log=TRUE)
		propdens_current = dtnorm(Yvals[ids], Yprop[ids], propysds[ids], lower=lowervals[ids], upper=uppervals[ids], log=TRUE)

		llprop=llpropnorm+llproppois
		ratio=exp(llprop-currentll + sum(propdens_current - propdens_prop) )
		ratio[ratio>1]=1
		if (runif(1)<=ratio){
			Yret=Yprop
			llretnorm=llpropnorm
			llretpois=llproppois
			accept= rep(1, length(which(ids)))
		}
		
	outlist=list(Yret, llretnorm,llretpois, accept)	
	names(outlist)=c("Y", "llnorm", "llpois","accepted")
	return(outlist)
	}
	
	########################################
	##### Initialize Sampling Matrices #####
	########################################

	params	= c(initpsis, initalpha)  #Everything here is on the transformed scale for proposals
	nparams	= length(params)
	propsds	= c(psipropsds)
	alphapropcovs = list(alpha0propcov, alpha1propcov)
	
	psiindex	= 1:2
	alpha0index = 3:(2+p_h)
	alpha1index = (max(alpha0index)+1) : (max(alpha0index)+ p_h)
	accepted 	= rep(0, nparams)
	accepted_y  = rep(0, n*q)
	accepted_B  = rep(0,p)

	binsize 	= 10
	rho		= 0
	notpd		= 0

	samples	= matrix(NA, nrow=length(index), ncol=nparams+length(B))
	dimnames(samples)[[2]]	= c(paste("Psi", 1:q, sep=""), paste("alpha0", (0:(p_h-1)), sep=""), paste("alpha1", (0:(p_h-1)), sep=""),
					     paste("B0", 0:((p/2)-1), sep=""), paste("B1", 0:((p/2)-1), sep=""))
	ysims		= matrix(NA, nrow=length(index), ncol=length(Y))
	hsims       = matrix(NA, nrow=length(index), ncol=length(H))

	##Calculate initial log likelihood
	lowertrunc = rep(NA, n*q)
	uppertrunc = lowertrunc
	lowertrunc[seq(1,n*q,q)] = Y[seq(2,n*q,q)]  #a=0 pollution has to be >= a=1 value
	uppertrunc[seq(1,n*q,q)] = Inf
	uppertrunc[seq(2,n*q,q)] = Y[seq(1,n*q,q)]  #a=1 pollution has to be <= a=0 value
	lowertrunc[seq(2,n*q,q)] = -Inf
	llnorm	= loglike_norm_monoind(Yvals=Y,Bvals=B,psivals=params[psiindex], lowervals=lowertrunc, uppervals=uppertrunc) 
	llpois  = loglike_pois(Y, params[c(alpha0index, alpha1index)], H[seq(1,n*2,2)], H[seq(2,n*2,2)])
	ll  = llnorm + llpois

	iterno	= 1
	donesampling= FALSE
	kk 		= 1
	
	while (donesampling==FALSE){
		
		#MH step for sampling B
		mhstep = MHstep_B(0,Y, B, params[psiindex],lowertrunc,uppertrunc,llnorm)
		B = mhstep$B
		llnorm = mhstep$ll
		accepted_B[1:(p/2)] = accepted_B[1:(p/2)] + mhstep$accepted
		if (llnorm == Inf) print(paste("Infinite Likelihood while/after updating B0 at iterno=", iterno))
		
		mhstep = MHstep_B(1,Y, B, params[psiindex], lowertrunc, uppertrunc,llnorm)
		B = mhstep$B
		llnorm = mhstep$ll
		accepted_B[((p/2)+1):p] = accepted_B[((p/2)+1):p] + mhstep$accepted
		if (llnorm == Inf) print(paste("Infinite Likelihood while/after updating B1 at iterno=", iterno))

		ll	= llnorm+llpois
	
		###### Metropolis steps for Psi #########
		mhstep	= MHstep_psi(whichpsi=0, Yvals=Y, Bvals=B, psivals=params[psiindex], lowertrunc, uppertrunc, currentll=llnorm)
		params[psiindex]	= mhstep$psi
		llnorm	= mhstep$ll
		accepted[psiindex[1]] = accepted[psiindex[1]]+mhstep$accepted
		if (llnorm == Inf) print(paste("Infinite Likelihood while/after updating Psi0 at iterno=", iterno))
	
		mhstep	= MHstep_psi(1, Y, B, params[psiindex], lowertrunc, uppertrunc, llnorm)
		params[psiindex]	= mhstep$psi
		llnorm	= mhstep$ll
		accepted[psiindex[2]] = accepted[psiindex[2]]+mhstep$accepted
		if (llnorm == Inf) print(paste("Infinite Likelihood while/after updating Psi1 at iterno=", iterno))
	
		ll 		= llnorm + llpois	
			
		###### Metropolis steps for Poisson parameters in h0/h1 model ##########
		#a0 parameters
		mhstep = MHstep_pois(0, Y, params[c(alpha0index, alpha1index)], H[seq(1,2*n,2)], H[seq(2,2*n,2)], llpois)
		params[c(alpha0index,alpha1index)] = mhstep$alpha
		llpois                             = mhstep$ll
		accepted[alpha0index]              = accepted[alpha0index] + mhstep$accepted
		
		#a1 parameters
		mhstep = MHstep_pois(1, Y, params[c(alpha0index, alpha1index)], H[seq(1,2*n,2)], H[seq(2,2*n,2)], llpois)
		params[c(alpha0index,alpha1index)] = mhstep$alpha
		llpois                             = mhstep$ll
		accepted[alpha1index]              = accepted[alpha1index] + mhstep$accepted
	
		##################################	
		###### Simulate Missing Y ########
		##################################
			
		######### Fully Bayesian sampling (with H) using R ##########
		for ( whichysamp in which(ismissy) ){	
		whichy=rep(FALSE, n*q)
		whichy[whichysamp]=TRUE
		mhstep = sampleY_monoind(whichy, Y, B, params[psiindex], params[c(alpha0index, alpha1index)], H[seq(1,n*2,2)], H[seq(2,n*2,2)], lowertrunc, uppertrunc, llpois, llnorm)
		Y      = mhstep$Y
		llnorm = mhstep$llnorm
		llpois = mhstep$llpois
		accepted_y[whichy] = accepted_y[whichy] + mhstep$accepted
		
		lowertrunc[seq(1,n*q,q)] = Y[seq(2,n*q,q)]  #a=0 pollution has to be >= a=1 value
		uppertrunc[seq(1,n*q,q)] = Inf
		uppertrunc[seq(2,n*q,q)] = Y[seq(1,n*q,q)]  #a=1 pollution has to be <= a=0 value
		lowertrunc[seq(2,n*q,q)] = -Inf
		}

		
		#Simulate Missing H
		H[seq(1,2*n,2)][data$a==1] = sampleH(0, params[c(alpha0index, alpha1index)], Y)
		H[seq(2,2*n,2)][data$a==0] = sampleH(1, params[c(alpha0index, alpha1index)], Y)
		llpois                     = loglike_pois(Y, params[c(alpha0index, alpha1index)], H[seq(1,n*2,2)], H[seq(2,n*q,2)])
		
		ll = llnorm + llpois
		
		if(iterno %in% index){
			samples[kk,]	= c(exp(params[psiindex]),params[c(alpha0index, alpha1index)], B) ##Transform log(Psi) to Psi
			ysims[kk,]		= Y
			hsims[kk,]      = H
			kk <- kk + 1
		}
	
		if (iterno%%binsize==0){
			print(paste("Iteration:", iterno))
			print(c("Accepted:", round(accepted/iterno,3), round(accepted_B/iterno,3), round(mean(accepted_y[ismissy]/iterno),3)))
		}
	
		if (iterno>=nsamp){donesampling=TRUE}	
	
		iterno = iterno+1

	}#while donesampling==FALSE
		
	out		<- list()
	out$samples <- samples
	out$y0	<- ysims[,seq(1,n*q,q)]
	out$y1 	<- ysims[,seq(2,n*q,q)]
	out$h0  <- hsims[,seq(1,n*2,2)]
	out$h1  <- hsims[,seq(2,n*2,2)]
	out$trt	<- data$a
	out$formula <- formula
	out$formula_h <- formula_h
	return(out)
}