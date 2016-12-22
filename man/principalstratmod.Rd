\name{principalstratmod}
\alias{principalstratmod}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit principal stratification model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
principalstratmod(formula, data, trt, coords, formula_h, denom, nsamp, nburn, thin, tuning = list(A = 0.1, psi = 0.2, theta = 1, alpha0 = alpha0prop, alpha1 = alpha1prop, Y = ypropsd), prior = list(KIG = rep(0.01, 2), psi = rep(0.01, 2), theta1 = rep(0.6, 2), theta2 = rep(10, 2)), starting = list(B = NULL, A = NULL, psi = NULL, theta = NULL, alpha0 = NULL, alpha1 = NULL))
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
  \item{formula_h}{
%%     ~~Describe \code{formula_h} here~~
}
  \item{denom}{
%%     ~~Describe \code{denom} here~~
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
function (formula, data, trt, coords, formula_h, denom, nsamp, 
    nburn, thin, tuning = list(A = 0.1, psi = 0.2, theta = 1, 
        alpha0 = alpha0prop, alpha1 = alpha1prop, Y = ypropsd), 
    prior = list(KIG = rep(0.01, 2), psi = rep(0.01, 2), theta1 = rep(0.6, 
        2), theta2 = rep(10, 2)), starting = list(B = NULL, A = NULL, 
        psi = NULL, theta = NULL, alpha0 = NULL, alpha1 = NULL)) 
{
    require(spBayes)
    require(corpcor)
    require(MASS)
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
    names_h <- colnames(model.matrix(formula_h, data))
    data$h <- data[, all.vars(formula_h)[1]]
    data$ytemp <- 1
    nterm <- length(labels(terms(formula_h)))
    terms <- NULL
    for (ii in 1:nterm) {
        if (ii == 1) 
            terms <- labels(terms(formula_h))[ii]
        if (ii > 1) 
            terms <- paste(terms, "+", labels(terms(formula_h))[ii], 
                sep = " ")
    }
    formulatemp <- as.formula(paste("ytemp ~ ", terms))
    covs_h <- model.matrix(formulatemp, data)
    makeXmat = function(Xvals) {
        p = dim(Xvals)[[2]]
        Xout = matrix(0, n * q, p * q)
        Xout[seq(1, n * q, q), 1:p] = Xvals
        Xout[seq(2, n * q, q), (p + 1):(2 * p)] = Xvals
        return(Xout)
    }
    X = makeXmat(covs)
    p = ncol(X)
    X_h = covs_h
    p_h = ncol(X_h) + 1
    y0 = rep(NA, n)
    y1 = rep(NA, n)
    y0[data$a == 0] = data$y[data$a == 0]
    y1[data$a == 1] = data$y[data$a == 1]
    Y = rep(-999, n * q)
    Y[seq(1, n * q, q)] = y0
    Y[seq(2, n * q, q)] = y1
    ismissy = (is.na(Y))
    nwithmissing = sum(ismissy)
    H = rep(-999, n * 2)
    H[seq(1, n * 2, 2)][data$a == 0] = data$h[data$a == 0]
    H[seq(2, n * 2, 2)][data$a == 1] = data$h[data$a == 1]
    denom = data[, denom]
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
    Y[seq(1, n * q, q)][missy1 == 1] = rnorm(sum(missy1 == 1), 
        mean(Y[seq(1, n * q, q)], na.rm = T), sd(Y[seq(1, n * 
            q, q)], na.rm = T))
    Y[seq(2, n * q, q)][missy2 == 1] = rnorm(sum(missy2 == 1), 
        mean(Y[seq(2, n * q, q)], na.rm = T), sd(Y[seq(2, n * 
            1, q)], na.rm = T))
    H[seq(1, n * 2, 2)][data$a == 1] = rpois(sum(data$a == 1), 
        10)
    H[seq(2, n * 2, 2)][data$a == 0] = rpois(sum(data$a == 0), 
        10)
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
    if (is.null(starting$alpha0) == T) {
        initalpha = rep(0, 2 * p_h)
    }
    else {
        initalpha = c(starting$alpha0, starting$alpha1)
    }
    A1propsds = tuning$A
    A2propsds = tuning$A
    psipropsds = rep(tuning$psi, 2)
    thetapropsds = rep(tuning$theta, 2)
    alpha0propcov = tuning$alpha0
    alpha1propcov = tuning$alpha1
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
    loglike_pois = function(Yvals, alphavals, h0vals, h1vals) {
        ya0 = Yvals[seq(1, n * q, q)]
        ya1 = Yvals[seq(2, n * q, q)]
        Xmat0 = as.matrix(cbind(X_h, ya0))
        alpha0 = alphavals[1:p_h]
        lambda0 = exp(Xmat0 \%*\% alpha0 + log(denom))
        Xmat1 = as.matrix(cbind(X_h, ya1))
        alpha1 = alphavals[(p_h + 1):length(alphavals)]
        lambda1 = exp(Xmat1 \%*\% alpha1 + log(denom))
        llike = sum(h0vals * log(lambda0) - lambda0 - lfactorial(h0vals))
        llike = llike + sum(h1vals * log(lambda1) - lambda1 - 
            lfactorial(h1vals))
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
    MHstep_pois = function(whichalphas, Yvals, alphavals, h0vals, 
        h1vals, currentll) {
        accept = 0
        llret = currentll
        alpharet = alphavals
        alphaprop = alphavals
        index = 1:p_h
        if (whichalphas == 1) {
            index = (p_h + 1):length(alphavals)
        }
        alphaprop[index] = mvrnorm(1, alpharet[index], alphapropcovs[[whichalphas + 
            1]])
        llprop = loglike_pois(Yvals, alphaprop, h0vals, h1vals)
        ratio = exp(llprop - currentll)
        ratio[ratio > 1] = 1
        if (runif(1) <= ratio) {
            alpharet = alphaprop
            llret = llprop
            accept = 1
        }
        outlist = list(alpharet, llret, accept)
        names(outlist) = c("alpha", "ll", "accepted")
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
    sampleH = function(whicha, alphavals, Yvals) {
        ya0 = Yvals[seq(1, n * q, q)]
        ya1 = Yvals[seq(2, n * q, q)]
        index = 1:p_h
        Xmat = cbind(X_h, ya0)
        if (whicha == 1) {
            index = (p_h + 1):length(alphavals)
            Xmat = cbind(X_h, ya1)
        }
        alphs = alphavals[index]
        lamda = exp(Xmat \%*\% alphs + log(denom))
        h = rpois(sum(data$a == (1 - whicha)), lamda[data$a == 
            (1 - whicha)])
        return(h)
    }
    sampleY_r = function(ids, Yvals, Bvals, Cinvval, logdetval, 
        alphavals, h0vals, h1vals, currentllpois, currentllnorm) {
        accept = rep(0, length(which(ids)))
        currentll = currentllnorm + currentllpois
        llretnorm = currentllnorm
        llretpois = currentllpois
        Yret = Yvals
        Yprop = Yvals
        Yprop[ids] = rnorm(length(which(ids)), Yvals[ids], propysds[ids])
        llpropnorm = loglike_norm(Yprop, Bvals, Cinvval, logdetval)
        llproppois = loglike_pois(Yprop, alphavals, h0vals, h1vals)
        llprop = llpropnorm + llproppois
        ratio = exp(llprop - currentll)
        ratio[ratio > 1] = 1
        if (runif(1) <= ratio) {
            Yret = Yprop
            llretnorm = llpropnorm
            llretpois = llproppois
            accept = rep(1, length(which(ids)))
        }
        outlist = list(Yret, llretnorm, llretpois, accept)
        names(outlist) = c("Y", "llnorm", "llpois", "accepted")
        return(outlist)
    }
    params = c(initA1, initA2, initpsis, inittheta, initalpha)
    nparams = length(params)
    propsds = c(A1propsds, A2propsds, psipropsds, thetapropsds)
    alphapropcovs = list(alpha0propcov, alpha1propcov)
    Aindex = 1:2
    psiindex = 3:4
    thetaindex = 5:6
    alpha0index = 7:(6 + p_h)
    alpha1index = (max(alpha0index) + 1):(max(alpha0index) + 
        p_h)
    accepted = rep(0, nparams)
    accepted_y = rep(0, n * q)
    binsize = 10
    rho = 0
    notpd = 0
    samples = matrix(NA, nrow = length(index), ncol = nparams + 
        1 + length(B))
    dimnames(samples)[[2]] = c("K[1,1]", "K[2,1]", "K[2,2]", 
        paste("Psi", 1:q, sep = ""), paste("Theta", 1:2, sep = ""), 
        paste("alpha0", (0:(p_h - 1)), sep = ""), paste("alpha1", 
            (0:(p_h - 1)), sep = ""), paste("B0", 0:((p/2) - 
            1), sep = ""), paste("B1", 0:((p/2) - 1), sep = ""))
    ysims = matrix(NA, nrow = length(index), ncol = length(Y))
    hsims = matrix(NA, nrow = length(index), ncol = length(H))
    initsp = loglike_sp(params, Y)
    Cinv = initsp$Cinv
    C = initsp$C
    logdet = initsp$logdet
    llsp = initsp$ll
    llnorm = loglike_norm(Y, B, Cinv, logdet)
    llpois = loglike_pois(Y, params[c(alpha0index, alpha1index)], 
        H[seq(1, n * 2, 2)], H[seq(2, n * 2, 2)])
    ll = llsp + llnorm + llpois
    iterno = 1
    donesampling = FALSE
    kk = 1
    while (donesampling == FALSE) {
        B = updateBeta(Cinv, Y)
        llnorm = loglike_norm(Y, B, Cinv, logdet)
        ll = llsp + llnorm + llpois
        paramindex = c(Aindex[1], psiindex[1], thetaindex[1])
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
        ll = llsp + llnorm + llpois
        mhstep = MHstep_pois(0, Y, params[c(alpha0index, alpha1index)], 
            H[seq(1, 2 * n, 2)], H[seq(2, 2 * n, 2)], llpois)
        params[c(alpha0index, alpha1index)] = mhstep$alpha
        llpois = mhstep$ll
        accepted[alpha0index] = accepted[alpha0index] + mhstep$accepted
        mhstep = MHstep_pois(1, Y, params[c(alpha0index, alpha1index)], 
            H[seq(1, 2 * n, 2)], H[seq(2, 2 * n, 2)], llpois)
        params[c(alpha0index, alpha1index)] = mhstep$alpha
        llpois = mhstep$ll
        accepted[alpha1index] = accepted[alpha1index] + mhstep$accepted
        for (whichysamp in which(ismissy)) {
            whichy = rep(FALSE, n * q)
            whichy[whichysamp] = TRUE
            mhstep = sampleY_r(whichy, Y, B, Cinv, logdet, params[c(alpha0index, 
                alpha1index)], H[seq(1, n * 2, 2)], H[seq(2, 
                n * 2, 2)], llpois, llnorm)
            Y = mhstep$Y
            llnorm = mhstep$llnorm
            llpois = mhstep$llpois
            accepted_y[whichy] = accepted_y[whichy] + mhstep$accepted
        }
        H[seq(1, 2 * n, 2)][data$a == 1] = sampleH(0, params[c(alpha0index, 
            alpha1index)], Y)
        H[seq(2, 2 * n, 2)][data$a == 0] = sampleH(1, params[c(alpha0index, 
            alpha1index)], Y)
        llpois = loglike_pois(Y, params[c(alpha0index, alpha1index)], 
            H[seq(1, n * 2, 2)], H[seq(2, n * q, 2)])
        ll = llsp + llnorm + llpois
        K1out = exp(params[Aindex[1]])
        K2out = exp(params[Aindex[2]])
        Kout = createspatialsig(K1out, K2out, rho)
        thetaout = logitInv(params[thetaindex], thetaunif_a, 
            thetaunif_b)
        if (iterno \%in\% index) {
            samples[kk, ] = c(Kout[lower.tri(diag(1, q), TRUE)], 
                exp(params[psiindex]), thetaout, params[c(alpha0index, 
                  alpha1index)], B)
            ysims[kk, ] = Y
            hsims[kk, ] = H
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
    out$h0 <- hsims[, seq(1, n * 2, 2)]
    out$h1 <- hsims[, seq(2, n * 2, 2)]
    out$coords <- coords
    out$trt <- data$a
    out$formula <- formula
    out$formula_h <- formula_h
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
