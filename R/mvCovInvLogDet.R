mvCovInvLogDet <- function(coords, knots, cov.model, V, Psi, theta, modified.pp=TRUE, SWM=TRUE, ...){
  
  
  ####################################################
  ##Check for unused args
  ####################################################
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
    if(! i %in% formal.args)
      warning("'",i, "' is not an argument")
  }
    
  ####################################################
  ##Distance matrices
  ####################################################
  
  ####################
  ##Coords
  ####################
  if(missing(coords)){stop("error: coords must be specified")}
  if(!is.matrix(coords) || ncol(coords) != 2){
    stop("error: coords is misspecified")
  }

  
  ####################################################
  ##Covariance model
  ####################################################
  if(missing(cov.model)){stop("error: cov.model must be specified")}
  if(!cov.model%in%c("gaussian","exponential","matern","spherical"))
    {stop("error: specified cov.model '",cov.model,"' is not a valid option; choose, from gaussian, exponential, matern, spherical.")}
  
      
    n <- nrow(coords)
    m <- nrow(V)
    mn <- m*n
    coords.D <- as.matrix(dist(coords))
    C.inv <- matrix(0,mn,mn)

    tmp.mm <- matrix(0,m,m)
    tmp.mm1 <- matrix(0,m,m)
    tmp.mm2 <- matrix(0,m,m)
    
    storage.mode(m) <- "integer"
    storage.mode(n) <- "integer"
    storage.mode(coords.D) <- "double"
    storage.mode(C.inv) <- "double"
    storage.mode(V) <- "double"
    storage.mode(Psi) <- "double"
    storage.mode(theta) <- "double"
    storage.mode(tmp.mm) <- "double"
    storage.mode(tmp.mm1) <- "double"
    storage.mode(tmp.mm2) <- "double"

    log.det <-  .Call("mvCovInvDet_wrap", coords.D, C.inv, n, m, Psi, V, theta, 
                      tmp.mm, tmp.mm1, tmp.mm2, cov.model, PACKAGE = "HEIfunctions")
        
    out <- list()
    out$log.det <- log.det
    out$C.inv <- C.inv
    out$C.inv[upper.tri(C.inv,F)] <- t(C.inv)[upper.tri(C.inv,F)]

    out$C <- .Call("mvCov", coords.D, n, m, Psi, V, theta, cov.model)

    return(out);
}
