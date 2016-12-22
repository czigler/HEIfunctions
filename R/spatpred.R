spatpred = function (formula, data, data.p, coords, coords.p, n.samples,uniqueid, 
                      nburn, ...) 
{
  require(spBayes)
  
  frame <- model.frame(formula, data)
  y <- model.extract(frame, "response")
  val <- FALSE
  if (length(y) < dim(data)[1]) 
    val <- TRUE
 
  m.1 <- spLM(formula, data = data, coords = coords, n.samples = n.samples, 
              , ...)

  ######### PREDICTION ###########
  data.p[names(frame)[1]] <- -99
  newcovs <- model.matrix(formula, data.p)
  out <- spPredict(m.1, pred.coords = coords.p, pred.covars = newcovs)
  
  prediction <- data.frame(Monitor = as.character(data.p$Monitor), 
                           predmean = apply(out$p.y.predictive.samples[, nburn:n.samples], 
                                            1, mean), predsd = apply(out$p.y.predictive.samples[, 
                                                                                                nburn:n.samples], 1, sd))
 
  results <- list()
  results$model <- m.1
  results$predmodel <- out
  results$pred <- prediction
  results$coords <- coords
  results$coords.p <- coords.p
  results$validate <- val
  class(results) <- "spatpredmodel"
  return(results)
}