\name{principalstratmod_monoind}
\alias{principalstratmod_monoind}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit principal stratification model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
principalstratmod_monoind(formula, data, trt, formula_h, denom, nsamp, nburn, thin, tuning = list(psi = 0.2, alpha0 = alpha0prop, alpha1 = alpha1prop, B0 = B0prop, B1 = B1prop, Y = ypropsd), prior = list(psi = rep(0.01, 2)), starting = list(B = NULL, psi = NULL, alpha0 = NULL, alpha1 = NULL))
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
function (formula, data, trt, formula_h, denom, nsamp, nburn, 
    thin, tuning = list(psi = 0.2, alpha0 = alpha0prop, alpha1 = alpha1prop, 
        B0 = B0prop, B1 = B1prop, Y = ypropsd), prior = list(psi = rep(0.01, 
        2)), starting = list(B = NULL, psi = NULL, alpha0 = NULL, 
        alpha1 = NULL)) 
{
    require(MCMCpack)
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
    psiig_a = prior$psi
    psiig_b = prior$psi
    my0 = summary(lm(y0 ~ data$ps))
    my1 = summary(lm(y1 ~ data$ps))
    missy1 = ismissy[seq(1, n * q, q)]
    missy2 = ismissy[seq(2, n * q, q)]
    missboth = rep(FALSE, n)
    missboth[is.na(Y[seq(1, n * q, q)]) & is.na(Y[seq(2, n * 
        q, q)])] = TRUE
    Y[seq(1, n * q, q)][missboth == TRUE] = rnorm(sum(missboth == 
        TRUE), data$base[missboth == TRUE], 0.2 * sd(dat$y, na.rm = TRUE))
    Y[seq(2, n * q, q)][missy2 == 1] = rtnorm(sum(missy2 == 1), 
        mean(Y[seq(2, n * q, q)], na.rm = T), 0.2 * sd(Y[seq(2, 
            n * 1, q)], na.rm = T), lower = -Inf, upper = Y[seq(1, 
            n * q, q)][missy2 == 1])
    Y[seq(1, n * q, q)][missy1 == 1] = rtnorm(sum(missy1 == 1), 
        mean(Y[seq(1, n * q, q)], na.rm = T), 0.2 * sd(Y[seq(1, 
            n * q, q)], na.rm = T), lower = Y[seq(2, n * q, q)][missy1 == 
            1], upper = Inf)
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
    if (is.null(starting$psi) == T) {
        initpsis = c(my0$sigma^2, my1$sigma^2)
    }
    else {
        initpsis = starting$psi
    }
    if (is.null(starting$alpha0) == T) {
        initalpha = rep(0, 2 * p_h)
    }
    else {
        initalpha = c(starting$alpha0, starting$alpha1)
    }
    psipropsds = tuning$psi
    alpha0propcov = tuning$alpha0
    alpha1propcov = tuning$alpha1
    B0propcov = tuning$B0
    B1propcov = tuning$B1
    propysds = tuning$Y
    loglike_norm_monoind = function(Yvals, Bvals, psivals) {
        lowertrunc = rep(NA, n * q)
        uppertrunc = lowertrunc
        lowertrunc[seq(1, n * q, q)] = Yvals[seq(2, n * q, q)]
        uppertrunc[seq(1, n * q, q)] = Inf
        uppertrunc[seq(2, n * q, q)] = Yvals[seq(1, n * q, q)]
        lowertrunc[seq(2, n * q, q)] = -Inf
        m = as.vector(X \%*\% Bvals)
        index = seq(1, n * q, q)
        llike = sum(dtnorm(Yvals[index], mean = m[index], sd = sqrt(psivals[1]), 
            lower = lowertrunc[index], upper = uppertrunc[index], 
            log = TRUE))
        index = seq(2, n * q, q)
        llike = llike + sum(dtnorm(Yvals[index], mean = m[index], 
            sd = sqrt(psivals[2]), lower = lowertrunc[index], 
            upper = uppertrunc[index], log = TRUE))
        llike = llike + log(dinvgamma(psivals[1], psiig_a[1], 
            psiig_b[1])) + log(dinvgamma(psivals[2], psiig_a[2], 
            psiig_b[2]))
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
    MHstep_B = function(Yvals, Bvals, psivals, currentll) {
        accept = 0
        llret = currentll
        Bret = Bvals
        Bprop = Bvals
        Bprop = mvrnorm(1, Bret, Bpropcov)
        llprop = loglike_norm_monoind(Yvals, Bprop, psivals)
        ratio = exp(llprop - currentll)
        if (runif(1) <= ratio) {
            Bret = Bprop
            llret = llprop
            accept = 1
        }
        outlist = list(Bret, llret, accept)
        names(outlist) = c("B", "ll", "accepted")
        return(outlist)
    }
    MHstep_psi = function(whichpsi, Yvals, Bvals, psivals, currentll) {
        accept = 0
        llret = currentll
        psiret = psivals
        psiprop = psivals
        psiprop[whichpsi + 1] = rnorm(1, psiret[whichpsi + 1], 
            psipropsds[whichpsi + 1])
        llprop = loglike_norm_monoind(Yvals, Bvals, psiprop)
        ratio = exp(llprop - currentll)
        if (runif(1) <= ratio) {
            psiret = psiprop
            llret = llprop
            accept = 1
        }
        outlist = list(psiret, llret, accept)
        names(outlist) = c("psi", "ll", "accepted")
        return(outlist)
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
    sampleY_monoind = function(ids, Yvals, Bvals, psivals, alphavals, 
        h0vals, h1vals, currentllpois, currentllnorm) {
        accept = rep(0, length(which(ids)))
        currentll = currentllnorm + currentllpois
        llretnorm = currentllnorm
        llretpois = currentllpois
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
        llpropnorm = loglike_norm_monoind(Yprop, Bvals, psivals)
        llproppois = loglike_pois(Yprop, alphavals, h0vals, h1vals)
        propdens_prop = dtnorm(Yprop[ids], Yvals[ids], propysds[ids], 
            lower = lowertrunc[ids], upper = uppertrunc[ids], 
            log = TRUE)
        propdens_current = dtnorm(Yvals[ids], Yprop[ids], propysds[ids], 
            lower = lowertrunc[ids], upper = uppertrunc[ids], 
            log = TRUE)
        llprop = llpropnorm + llproppois
        ratio = exp(llprop - currentll + sum(propdens_current - 
            propdens_prop))
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
    params = c(initpsis, initalpha)
    nparams = length(params)
    propsds = c(psipropsds)
    alphapropcovs = list(alpha0propcov, alpha1propcov)
    Bpropcov = matrix(0, p, p)
    Bpropcov[1:(p/2), 1:(p/2)] = B0propcov
    Bpropcov[((p/2) + 1):p, ((p/2) + 1):p] = B1propcov
    psiindex = 1:2
    alpha0index = 3:(2 + p_h)
    alpha1index = (max(alpha0index) + 1):(max(alpha0index) + 
        p_h)
    accepted = rep(0, nparams)
    accepted_y = rep(0, n * q)
    accepted_B = rep(0, p)
    binsize = 10
    rho = 0
    notpd = 0
    samples = matrix(NA, nrow = length(index), ncol = nparams + 
        length(B))
    dimnames(samples)[[2]] = c(paste("Psi", 1:q, sep = ""), paste("alpha0", 
        (0:(p_h - 1)), sep = ""), paste("alpha1", (0:(p_h - 1)), 
        sep = ""), paste("B0", 0:((p/2) - 1), sep = ""), paste("B1", 
        0:((p/2) - 1), sep = ""))
    ysims = matrix(NA, nrow = length(index), ncol = length(Y))
    hsims = matrix(NA, nrow = length(index), ncol = length(H))
    llnorm = loglike_norm_monoind(Yvals = Y, Bvals = B, psivals = params[psiindex])
    llpois = loglike_pois(Y, params[c(alpha0index, alpha1index)], 
        H[seq(1, n * 2, 2)], H[seq(2, n * 2, 2)])
    ll = llnorm + llpois
    iterno = 1
    donesampling = FALSE
    kk = 1
    while (donesampling == FALSE) {
        mhstep = MHstep_B(Y, B, params[psiindex], llnorm)
        B = mhstep$B
        llnorm = mhstep$ll
        accepted_B = accepted_B + mhstep$accepted
        if (llnorm == Inf) 
            print(paste("Infinite Likelihood while/after updating B at iterno=", 
                iterno))
        ll = llnorm + llpois
        mhstep = MHstep_psi(whichpsi = 0, Yvals = Y, Bvals = B, 
            psivals = params[psiindex], currentll = llnorm)
        params[psiindex] = mhstep$psi
        llnorm = mhstep$ll
        accepted[psiindex[1]] = accepted[psiindex[1]] + mhstep$accepted
        if (llnorm == Inf) 
            print(paste("Infinite Likelihood while/after updating Psi0 at iterno=", 
                iterno))
        mhstep = MHstep_psi(1, Y, B, params[psiindex], llnorm)
        params[psiindex] = mhstep$psi
        llnorm = mhstep$ll
        accepted[psiindex[2]] = accepted[psiindex[2]] + mhstep$accepted
        if (llnorm == Inf) 
            print(paste("Infinite Likelihood while/after updating Psi1 at iterno=", 
                iterno))
        ll = llnorm + llpois
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
            mhstep = sampleY_monoind(whichy, Y, B, params[psiindex], 
                params[c(alpha0index, alpha1index)], H[seq(1, 
                  n * 2, 2)], H[seq(2, n * 2, 2)], llpois, llnorm)
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
        ll = llnorm + llpois
        if (iterno \%in\% index) {
            samples[kk, ] = c(params[psiindex], params[c(alpha0index, 
                alpha1index)], B)
            ysims[kk, ] = Y
            hsims[kk, ] = H
            kk <- kk + 1
        }
        if (iterno\%\%binsize == 0) {
            print(paste("Iteration:", iterno))
            print(c("Accepted:", round(accepted/iterno, 3), round(mean(accepted_y[ismissy]/iterno), 
                3)))
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
