\name{cmratioiacov}
\alias{cmratioiacov}
\title{
Contrast objects for treatment comparison over a grid of covariate values
}
\description{
 Compute contrast object for pairwise comparisons of treatment levels over grids of covariates given a fitted model object of class 'lm' or 'glm', such that it is suitable for further use in function gsci.ratio, package mratios.
 Restricted to ONE factor variable and pairwise comparisons of treatments ("Dunnett", "Tukey", or user defined subsets of Tukey-type comparisons). More than one covariate, quadratic terms etc. are possible.
}
\usage{
cmratioiacov(modelfit, treatment, treatcon = "Dunnett", covset, ...)
}

\arguments{
  \item{modelfit}{a fitted model of class \code{"lm"} or \code{"glm"}, right hand side of the model formula should contain one factor variable, and one (or several) covariate(s), etc.}
  \item{treatment}{a single character string, identifying the factor variable in the model}
  \item{treatcon}{
either a list of two matrices with integer entries, or a single character string. If a list of two matrices: number of columns must be = number of treatment levels in the factor variable; rows correspond to comparisons between treatments: only coefficients 1,0 are allowed, all contrasts must define pairwise treatment differences; if a character string: only "Dunnett" or "Tukey" are currently allowed, these will be passed to the \code{type}-argument in \code{contrMatRatio} (package \code{mratios})
}
  \item{covset}{a named list with numeric vectors: names of list elements must be the names of the covariates in the model, vector in each element is supposed to contain the covariate values at which contrasts sould be computed}
  \item{\dots}{further arguments to be passed to internal functions, currently only:
\code{base} a single integer indicating the control level, if \code{treatcon} is a character string, passed to \code{contrMat}
\code{digits} a single integer, defining the number of visible non-zero digits when \code{round}ing covariate values in contrast names
}
}
\details{
Based on the call and data set in the object \code{modelfit}, matrices are contructed that allow to compute predictions over the grid of covariate values specified in \code{covset} and the factor variable specified in \code{treatment}. These matrices are then combined to define pairwise differences between the levels in treatment, over the grid of covariate values. The type (or set) of pairwise differences can be defined in the argument \code{treatcon}.
}
\value{
A list with items
\item{numC}{matrix of linear combinations of model parameters that define the model predictions to appear in the numerator of the multiple ratios}
\item{denC}{matrix of linear combinations of model parameters that define the model predictions to appear in the denominator of the multiple ratios}
\item{rnames}{vector of names for the comparisons defined by \code{numC} and \code{denC}}
\item{dataratio}{ a data set with names for the ratios, the treatment levels contributing in numerator and denominator of the ratios, and the covariate(s); rows of the data set are ordered as are the rows in the numerator and denominator matrices, \code{denC} and \code{numC}}
\item{Xpos}{model matrix for those predictions used as minuend to compute \code{linfct}}
\item{Xneg}{model matrix for those predictions used as subtrahend to compute \code{linfct}}
\item{datapos}{}
\item{dataneg}{}
}

\author{Frank Schaarschmidt}

\examples{

data(pc)

fitPC<-lm(yield ~ x + I(x^2) + additive + additive:I(x^2), data=pc)

anova(fitPC)

# estimate the ratio of model predictions over a grid of covariates,
# where the covariate x takes vales 1,2,3,...9
# and the ratios of additive S1 to control and S2 to control
# are of interest (Dunnett-type-comparisons)

cmrpc <- cmratioiacov(fitPC, treatment="additive",
   covset=list("x" = seq(from=1, to=9, by=1)), treatcon="Dunnett")

str(cmrpc, max.level=1)

# numerator and denominator matrices for plug-in in gsci.ratio
str(cmrpc$denC)
str(cmrpc$numC)

# a vector of names for the comparisons
str(cmrpc$rnames)

# a data set containing names of comparisons, treatment levels
# contributing and covariate values(for later use in plotting)
str(cmrpc$rnames)

}

\keyword{misc}
