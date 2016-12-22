\name{plot.spatpredmodel}
\alias{plot.spatpredmodel}
\title{
Plotting function for spatial predicitons.
}
\description{
Makes multiple plots to summarize baseline spatial predictions. 
}
\usage{
plot.spatpredmodel(outmodel, origdat = NULL, preddat = NULL, covname = NULL, plotnum = 1)
}
\arguments{
  \item{outmodel}{model with baseline spatial predictions, of class spatpredmodel}
  \item{origdat}{original dataset}
  \item{preddat}{dataset for prediction}
  \item{covname}{string with covariate name in predicted dataset (used for validation)}
  \item{plotnum}{number specifying which plot to draw, see details}
}
\details{
Plot numbers:
1. Traceplots and posterior density from original data model.
2. Extrapolated surface of original outcome data (not that interesting, just gives an idea of original spatial structure). 
3. Use if you are validating the model only.  Plot of true versus predicted mean pollution from model.  
4. Use if you are validating the model only.  Plot of standard deviation of predicted pollution
 samples versus standard deviation of predicted means from true pollution. 
5. Use if you are validating the model only.  Shows the location of the worst prediction in validation model. 
}
\author{
Lauren
}
