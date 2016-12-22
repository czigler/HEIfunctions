\name{mvpsmod_ind}
\alias{mvpsmod_ind}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MV linear model, cory edit.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mvpsmod_ind(formula, data, trt, nsamp, nburn, thin, prior = list(psi = rep(0.01, 2)), starting = list(B = NULL, psi = NULL))
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
  \item{nsamp}{
%%     ~~Describe \code{nsamp} here~~
}
  \item{nburn}{
%%     ~~Describe \code{nburn} here~~
}
  \item{thin}{
%%     ~~Describe \code{thin} here~~
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
function (formula, data, trt, nsamp, nburn, thin, prior = list(psi = rep(0.01, 
    2)), starting = list(B = NULL, psi = NULL)) 
{
    require(MCMCpack)
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
    psiig_a = prior$psi
    psiig_b = prior$psi
    my0 = summary(lm(y0 ~ data$ps))
    my1 = summary(lm(y1 ~ data$ps))
    missy1 = ismissy[seq(1, n * q, q)]
    missy2 = ismissy[seq(2, n * q, q)]
    missingindy = list(missy1, missy2)
    Y[seq(1, n * q, q)][missy1 == 1] = rnorm(sum(missy1 == 1), 
        mean(Y[seq(1, n * q, q)], na.rm = T), sd(Y[seq(1, n * 
            q, q)], na.rm = T))
    Y[seq(2, n * q, q)][missy2 == 1] = rnorm(sum(missy2 == 1), 
        mean(Y[seq(2, n * q, q)], na.rm = T), sd(Y[seq(2, n * 
            1, q)], na.rm = T))
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
    updateBeta = function(whichbetas, psivals, Yvals) {
        index1 = seq(1, n * q, q)
        index2 = 1:(p/2)
        if (whichbetas == 1) {
            index1 = seq(2, n * q, q)
            index2 = ((p/2) + 1):p
        }
        S_beta = solve(t(X[index1, index2]) \%*\% X[index1, index2]) * 
            psivals[whichbetas + 1]
        Mu_beta = S_beta \%*\% t(X[index1, index2]) \%*\% Yvals[index1]/psivals[whichbetas + 
            1]
        return(mvrnorm(1, Mu_beta, S_beta))
    }
    updatepsi = function(whichpsi, Bvals, Yvals) {
        index1 = seq(1, n * q, q)
        if (whichpsi == 1) 
            index1 = seq(2, n * q, q)
        YXB = (Yvals[index1] - X[index1, ] \%*\% Bvals)
        anew = psiig_a[whichpsi + 1] + n/2
        bnew = psiig_b[whichpsi + 1] + 0.5 * (t(YXB) \%*\% YXB)
        return(rinvgamma(1, anew, bnew))
    }
    params = c(initpsis)
    nparams = length(params)
    psiindex = 1:2
    binsize = 500
    rho = 0
    notpd = 0
    samples = matrix(NA, nrow = length(index), ncol = nparams + 
        length(B))
    dimnames(samples)[[2]] = c(paste("Psi", 1:q, sep = ""), paste("B0", 
        0:((p/2) - 1), sep = ""), paste("B1", 0:((p/2) - 1), 
        sep = ""))
    ysims = matrix(NA, nrow = length(index), ncol = length(Y))
    iterno = 1
    donesampling = FALSE
    kk = 1
    while (donesampling == FALSE) {
        B[1:(p/2)] = updateBeta(whichbetas = 0, psivals = params[psiindex], 
            Yvals = Y)
        B[((p/2) + 1):p] = updateBeta(1, params[psiindex], Y)
        params[psiindex[1]] = updatepsi(whichpsi = 0, B, Y)
        params[psiindex[2]] = updatepsi(1, B, Y)
        index1 = seq(1, n * q, q)
        m_pred = X[index1, ] \%*\% B
        Y[index1][ismissy[index1]] = rnorm(sum(ismissy[index1]), 
            m_pred[ismissy[index1]], sqrt(params[psiindex[1]]))
        index1 = seq(2, n * q, q)
        m_pred = X[index1, ] \%*\% B
        Y[index1][ismissy[index1]] = rnorm(sum(ismissy[index1]), 
            m_pred[ismissy[index1]], sqrt(params[psiindex[2]]))
        if (iterno \%in\% index) {
            samples[kk, ] = c(params[psiindex], B)
            ysims[kk, ] = Y
            kk <- kk + 1
        }
        if (iterno\%\%binsize == 0) {
            print(paste("Iteration:", iterno))
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
    out$trt <- data$a
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
