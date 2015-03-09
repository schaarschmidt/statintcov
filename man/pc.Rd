\name{pc}
\alias{pc}
\docType{data}
\title{
Yield of a chemical process depending on the amount of a substance and presence of two additves.
}
\description{
A simulated data set (similar to a data set presented by Milliken and Johnson, 2002), in which an interaction of a treatment effect with the quadratic term of a covariate is plausible.
}
\usage{data("pc")}
\format{
  A data frame with 36 observations on the following 3 variables.
  \describe{
    \item{\code{yield}}{a numeric vector, observed yield of a chemical process, the response variable of interest}
    \item{\code{x}}{a numeric vector, substance on which the chemical porocess depends, a potential covariate }
    \item{\code{additive}}{a factor with levels \code{Control} \code{S1} \code{S2}, reprenting the treatments of interest: the process run without additive (control), and the process run in presence of two different additives}
  }
}
\details{
Milliken and Johnson (2002) discuss the analysis of this data set using ANCOVA. It is obvious, that a quadratic term can be used to model the \code{yield}s dependency on the covariante \code{x}. Milliken and Johnson show that either an interaction of the \code{additive} with \code{x} and a common quadratic term meay be used, or a common linear term and an interaction of the \code{additive} with a quadratic term for \code{x} may be used. The latter model is used here to illustrate the comparison of treatments in presence of treatment- covariate interaction, when interaction is with the quadratic term.
}

\source{
The data set is simulated to resemble a data set presented in
Milliken and Johnson, 2002: Table 4.9, page 108; model: Table 4.11, 4.13, page 110.
}

\examples{

data("pc")

fitpc<-lm(yield ~ x + I(x^2) + additive + additive:I(x^2), data=pc)
anova(fitpc)

if(require("ggplot2")){
ggplot(pc, aes(y=yield, x=x)) + geom_point() + facet_wrap(~additive) 
}

}
\keyword{datasets}
