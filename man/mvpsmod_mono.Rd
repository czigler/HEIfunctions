\name{mvpsmod_mono}
\alias{mvpsmod_mono}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MV linear model, cory edit.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mvpsmod_mono(formula, data, trt, coords, nsamp, nburn, thin, tuning = list(A = 0.1, psi = 0.2, theta = 1, Y = 1), prior = list(KIG = rep(0.01, 2), psi = rep(0.01, 2), theta1 = rep(0.6, 2), theta2 = rep(10, 2)), starting = list(B = NULL, A = NULL, psi = NULL, theta = NULL))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
%%     ~~Describe \code{formula} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{trt}{
%%     ~~Describe \code{trt} here~~
}
  \item{coords}{
%%     ~~Describe \code{coords} here~~
}
  \item{nsamp}{
%%     ~~Describe \code{nsamp} here~~
}
  \item{nburn}{
%%     ~~Describe \code{nburn} here~~
}
  \item{thin}{
%%     ~~Describe \code{thin} here~~
}
  \item{tuning}{
%%     ~~Describe \code{tuning} here~~
}
  \item{prior}{
%%     ~~Describe \code{prior} here~~
}
  \item{starting}{
%%     ~~Describe \code{starting} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (formula, data, trt, coords, nsamp, nburn, thin, tuning = list(A = 0.1, 
    psi = 0.2, theta = 1, Y = 1), prior = list(KIG = rep(0.01, 
    2), psi = rep(0.01, 2), theta1 = rep(0.6, 2), theta2 = rep(10, 
    2)), starting = list(B = NULL, A = NULL, psi = NULL, theta = NULL)) 
{
    require(spBayes)
    require(corpcor)
    require(MASS)
    require(msm)
    ncov <- dim(model.matrix(formula, data))[2]
    print("Note: function does not accept missing data in covariates.")
    logit = function(theta, a, b) {
        return(log((theta - a)/(b - theta)))
    }
    logitInv = function(z, a, b) {
        return(b - (b - a)/(1 + exp(z)))
    }
    index <- seq(from = (nburn + 1), to = nsamp, by = thin)
    n <- dim(data)[1]
    q <- 2
    names <- colnames(model.matrix(formula, data))
    data$a <- data[, trt]
    data$y <- data[, all.vars(formula)[1]]
    data$ytemp <- 1
    nterm <- length(labels(terms(formula)))
    terms <- NULL
    for (ii in 1:nterm) {
        if (ii == 1) 
            terms <- labels(terms(formula))[ii]
        if (ii > 1) 
            terms <- paste(terms, "+", labels(terms(formula))[ii], 
                sep = " ")
    }
    formulatemp <- as.formula(paste("ytemp ~ ", terms))
    covs <- model.matrix(formulatemp, data)
    makeXmat = function(Xvals) {
        p = dim(Xvals)[[2]]
        Xout = matrix(0, n * q, p * q)
        Xout[seq(1, n * q, q), 1:p] = Xvals
        Xout[seq(2, n * q, q), (p + 1):(2 * p)] = Xvals
        return(Xout)
    }
    X = makeXmat(covs)
    p = ncol(X)
    y0 = rep(NA, n)
    y1 = rep(NA, n)
    y0[data$a == 0] = data$y[data$a == 0]
    y1[data$a == 1] = data$y[data$a == 1]
    Y = rep(-999, n * q)
    Y[seq(1, n * q, q)] = y0
    Y[seq(2, n * q, q)] = y1
    ismissy = (is.na(Y))
    nwithmissing = sum(ismissy)
    a = rep(NA, n * q)
    a[seq(1, n * q, q)] = data$a
    a[seq(2, n * q, q)] = data$a
    KIG_a = prior$KIG
    KIG_b = prior$KIG
    psiig_a = prior$psi
    psiig_b = prior$psi
    thetaunif_a = prior$theta1
    thetaunif_b = prior$theta2
    my0 = summary(lm(y0 ~ data$ps))
    my1 = summary(lm(y1 ~ data$ps))
    missy1 = ismissy[seq(1, n * q, q)]
    missy2 = ismissy[seq(2, n * q, q)]
    missingindy = list(missy1, missy2)
    missboth = rep(FALSE, n)
    missboth[is.na(Y[seq(1, n * q, q)]) & is.na(Y[seq(2, n * 
        q, q)])] = TRUE
    Y[seq(1, n * q, q)][missboth == TRUE] = rnorm(sum(missboth == 
        TRUE), data$base[missboth == TRUE], sd(dat$y, na.rm = TRUE))
    Y[seq(2, n * q, q)][missy2 == 1] = rtnorm(sum(missy2 == 1), 
        mean(Y[seq(2, n * q, q)], na.rm = T), sd(Y[seq(2, n * 
            1, q)], na.rm = T), lower = -Inf, upper = Y[seq(1, 
            n * q, q)][missy2 == 1])
    Y[seq(1, n * q, q)][missy1 == 1] = rtnorm(sum(missy1 == 1), 
        mean(Y[seq(1, n * q, q)], na.rm = T), sd(Y[seq(1, n * 
            q, q)], na.rm = T), lower = Y[seq(2, n * q, q)][missy1 == 
            1], upper = Inf)
    if (is.null(starting$B) == F) {
        B = starting$B
    }
    if (is.null(starting$B) == T) {
        B = rep(0, ncov * 2)
    }
    if (is.null(starting$A) == T) {
        initA1 = log(my0$sigma^2)
        initA2 = log(my1$sigma^2)
    }
    else {
        initA1 = starting$A
        initA2 = starting$A
    }
    if (is.null(starting$psi) == T) {
        initpsis = log(c(my0$sigma^2, my1$sigma^2) * 0.1)
    }
    else {
        initpsis = log(rep(starting$psi, 2))
    }
    if (is.null(starting$theta) == T) {
        inittheta = c(1, 1)
        inittheta = logit(inittheta, thetaunif_a, thetaunif_b)
    }
    else {
        inittheta = logit(starting$theta, thetaunif_a, thetaunif_b)
    }
    A1propsds = tuning$A
    A2propsds = tuning$A
    psipropsds = rep(tuning$psi, 2)
    thetapropsds = rep(tuning$theta, 2)
    propysds = tuning$Y
    loglike_sp = function(paramvals, Yvals) {
        llike = 0
        K1 = exp(paramvals[Aindex[1]])
        K2 = exp(paramvals[Aindex[2]])
        K = createspatialsig(K1, K2, rho)
        if (!is.na(K)[[1]]) {
            thetatemp = logitInv(paramvals[thetaindex], thetaunif_a, 
                thetaunif_b)
            det = mvCovInvLogDet(coords = coords, cov.model = "exponential", 
                V = K, Psi = diag(exp(paramvals[psiindex]), nrow = q), 
                theta = thetatemp, modified.pp = FALSE, SWM = TRUE)
            llike = llike + sum(-(KIG_a + 1) * paramvals[Aindex] - 
                KIG_b/exp(paramvals[Aindex]) + paramvals[Aindex])
            llike = llike + sum(-(psiig_a + 1) * paramvals[psiindex] - 
                psiig_b/exp(paramvals[psiindex]) + paramvals[psiindex])
            llike = llike + sum(log(thetatemp - thetaunif_a) + 
                log(thetaunif_b - thetatemp))
            outlist = list(llike, det$C, det$C.inv, det$log.det)
            names(outlist) = c("ll", "C", "Cinv", "logdet")
        }
        if (is.na(K)[[1]]) {
            outlist = list(-Inf, 1)
            names(outlist) = c("ll", "notpd")
        }
        return(outlist)
    }
    loglike_norm = function(Yvals, Bvals, Cinvval, logdetval) {
        YXB = (Yvals - X \%*\% Bvals)
        llike = -0.5 * logdetval - 0.5 * t(YXB) \%*\% Cinvval \%*\% 
            YXB
        return(llike)
    }
    MHstep_sp = function(index, paramvals, Bvals, Yvals, Cval, 
        Cinvval, logdetval, currentllsp, currentllnorm) {
        accept = 0
        notpd = 1
        llretsp = currentllsp
        llretnorm = currentllnorm
        currentll = currentllsp + currentllnorm
        Cret = Cval
        Cinvret = Cinvval
        logdetret = logdetval
        currentparams = paramvals
        props = currentparams
        props[index] = rnorm(length(index), paramvals[index], 
            propsds[index])
        llanddet = loglike_sp(props, Yvals)
        Cprop = llanddet$C
        Cinvprop = llanddet$Cinv
        logdetprop = llanddet$logdet
        llpropsp = llanddet$ll
        if (llanddet$ll != -Inf) {
            notpd = 0
            llpropnorm = loglike_norm(Yvals, Bvals, Cinvprop, 
                logdetprop)
            llprop = llpropsp + llpropnorm
            ratio = exp(llprop - currentll)
            ratio[ratio > 1] = 1
            if (runif(1) <= ratio) {
                currentparams[index] = props[index]
                llretsp = llpropsp
                llretnorm = llpropnorm
                Cret = Cprop
                Cinvret = Cinvprop
                logdetret = logdetprop
                accept = 1
            }
        }
        outlist = list(currentparams, llretsp, llretnorm, accept, 
            Cret, Cinvret, logdetret, notpd)
        names(outlist) = c("params", "llsp", "llnorm", "accepted", 
            "C", "Cinv", "logdet", "notpd")
        return(outlist)
    }
    updateBeta = function(Cinvval, Yvals) {
        S_beta = solve(t(X) \%*\% Cinvval \%*\% X)
        Mu_beta = S_beta \%*\% t(X) \%*\% Cinvval \%*\% Yvals
        return(mvrnorm(1, Mu_beta, S_beta))
    }
    createspatialsig = function(K1mat, K2mat, rho) {
        sig = matrix(0, q, q)
        sig[1, 1] = K1mat
        sig[2, 2] = K2mat
        sig[1, 2] = rho * sqrt(K1mat * K2mat)
        sig[2, 1] = sig[1, 2]
        if (!is.positive.definite(sig)) {
            sig = NA
        }
        return(sig)
    }
    sampleY_mono = function(ids, Yvals, Bvals, Cinvval, logdetval, 
        currentllnorm) {
        accept = rep(0, length(which(ids)))
        currentll = currentllnorm
        llretnorm = currentllnorm
        Yret = Yvals
        Yprop = Yvals
        lowertrunc = rep(NA, n * q)
        uppertrunc = lowertrunc
        lowertrunc[seq(1, n * q, q)] = Yvals[seq(2, n * q, q)]
        uppertrunc[seq(1, n * q, q)] = Inf
        uppertrunc[seq(2, n * q, q)] = Yvals[seq(1, n * q, q)]
        lowertrunc[seq(2, n * q, q)] = -Inf
        Yprop[ids] = rtnorm(length(which(ids)), Yvals[ids], propysds[ids], 
            lower = lowertrunc[ids], upper = uppertrunc[ids])
        llpropnorm = loglike_norm(Yprop, Bvals, Cinvval, logdetval)
        propdens_prop = dtnorm(Yprop[ids], Yvals[ids], propysds[ids], 
            lower = lowertrunc[ids], upper = uppertrunc[ids], 
            log = TRUE)
        propdens_current = dtnorm(Yvals[ids], Yvals[ids], propysds[ids], 
            lower = lowertrunc[ids], upper = uppertrunc[ids], 
            log = TRUE)
        ratio = exp(llpropnorm - currentllnorm + sum(propdens_current - 
            propdens_prop))
        ratio[ratio > 1] = 1
        if (runif(1) <= ratio) {
            Yret = Yprop
            llretnorm = llpropnorm
            accept = rep(1, length(which(ids)))
        }
        outlist = list(Yret, llretnorm, accept)
        names(outlist) = c("Y", "llnorm", "accepted")
        return(outlist)
    }
    params = c(initA1, initA2, initpsis, inittheta)
    nparams = length(params)
    propsds = c(A1propsds, A2propsds, psipropsds, thetapropsds)
    Aindex = 1:2
    psiindex = 3:4
    thetaindex = 5:6
    accepted = rep(0, length(propsds))
    accepted_y = rep(0, n * q)
    binsize = 10
    rho = 0
    notpd = 0
    samples = matrix(NA, nrow = length(index), ncol = nparams + 
        1 + length(B))
    dimnames(samples)[[2]] = c("K[1,1]", "K[2,1]", "K[2,2]", 
        paste("Psi", 1:q, sep = ""), paste("Theta", 1:2, sep = ""), 
        paste("B0", 0:((p/2) - 1), sep = ""), paste("B1", 0:((p/2) - 
            1), sep = ""))
    ysims = matrix(NA, nrow = length(index), ncol = length(Y))
    initsp = loglike_sp(params, Y)
    Cinv = initsp$Cinv
    C = initsp$C
    logdet = initsp$logdet
    llsp = initsp$ll
    llnorm = loglike_norm(Y, B, Cinv, logdet)
    ll = llsp + llnorm
    iterno = 1
    donesampling = FALSE
    kk = 1
    while (donesampling == FALSE) {
        B = updateBeta(Cinv, Y)
        llnorm = loglike_norm(Y, B, Cinv, logdet)
        ll = llsp + llnorm
        paramindex = c(Aindex[1], psiindex[1], thetaindex[1])
        mhstep = MHstep_sp(index = paramindex, paramvals = params, 
            Bvals = B, Yvals = Y, Cval = C, Cinvval = Cinv, logdetval = logdet, 
            currentllsp = llsp, currentllnorm = llnorm)
        params = mhstep$params
        llsp = mhstep$llsp
        llnorm = mhstep$llnorm
        Cinv = mhstep$Cinv
        C = mhstep$C
        logdet = mhstep$logdet
        accepted[paramindex] = accepted[paramindex] + mhstep$accepted
        notpd = notpd + mhstep$notpd
        paramindex = c(Aindex[2], psiindex[2], thetaindex[2])
        mhstep = MHstep_sp(paramindex, params, B, Y, C, Cinv, 
            logdet, llsp, llnorm)
        params = mhstep$params
        llsp = mhstep$llsp
        llnorm = mhstep$llnorm
        Cinv = mhstep$Cinv
        C = mhstep$C
        logdet = mhstep$logdet
        accepted[paramindex] = accepted[paramindex] + mhstep$accepted
        notpd = notpd + mhstep$notpd
        ll = llsp + llnorm
        for (whichysamp in which(ismissy)) {
            whichy = rep(FALSE, n * q)
            whichy[whichysamp] = TRUE
            mhstep = sampleY_mono(ids = whichy, Yvals = Y, Bvals = B, 
                Cinvval = Cinv, logdetval = logdet, currentllnorm = llnorm)
            Y = mhstep$Y
            llnorm = mhstep$llnorm
            accepted_y[whichy] = accepted_y[whichy] + mhstep$accepted
        }
        ll = llsp + llnorm
        K1out = exp(params[Aindex[1]])
        K2out = exp(params[Aindex[2]])
        Kout = createspatialsig(K1out, K2out, rho)
        thetaout = logitInv(params[thetaindex], thetaunif_a, 
            thetaunif_b)
        if (iterno \%in\% index) {
            samples[kk, ] = c(Kout[lower.tri(diag(1, q), TRUE)], 
                exp(params[psiindex]), thetaout, B)
            ysims[kk, ] = Y
            kk <- kk + 1
        }
        if (iterno\%\%binsize == 0) {
            print(paste("Iteration:", iterno))
            print(c("Accepted:", round(accepted/iterno, 3), round(mean(accepted_y[ismissy]/iterno), 
                3)))
            print(c("Number of not PD:", notpd))
        }
        if (iterno >= nsamp) {
            donesampling = TRUE
        }
        iterno = iterno + 1
    }
    out <- list()
    out$samples <- samples
    out$y0 <- ysims[, seq(1, n * q, q)]
    out$y1 <- ysims[, seq(2, n * q, q)]
    out$coords <- coords
    out$trt <- data$a
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
