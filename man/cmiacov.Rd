\name{cmiacov}
\alias{cmiacov}

\title{
Contrast matrices for treatments over a grid of covariate values
}

\description{
 Compute suitable contrast matrix for pairwise comparisons of treatment levels over grids of covariates given a fitted model object of class 'lm' or 'glm'
 Restricted to ONE factor variable and pairwise comparisons of treatments ("Dunnett", "Tukey", or user defined subsets of Tukey-type comparisons). More than one covariate, quadratic terms etc. are possible
}
\usage{
cmiacov(modelfit, treatment, treatcon = "Dunnett", covset, ...)
}

\arguments{
  \item{modelfit}{a fitted model of class \code{"lm"} or \code{"glm"}, right hand side of the model formula should contain one factor variable, and one (or several) covariate(s), etc.}
  \item{treatment}{a single character string, identifying the factor variable in the model}
  \item{treatcon}{either a matrix with integer entries, or a single character string. If a matrix: number of columns must be = number of treatment levels in the factor variable rows specifies contrasts between treatments: only coefficients 1,0,-1 are allowed, all contrasts must define pairwise treatment differences: if a character string: "Dunnett" or "Tukey", will be passed to the \code{type}-argument in \code{contrMat} (package \code{multcomp})}
  \item{covset}{a named list with numeric vectors: names of list elements must be the names of the covariates in the model, vector in each element is supposed to contain the covariate values at which contrasts sould be computed}
  \item{\dots}{further arguments to be passed to internal functions, currently only:
\code{base} a single integer indicating the control level, if \code{treatcon} is a character string, passed to \code{contrMat}
\code{digits} a single integer, defining the number of visible non-zero digits when \code{round}ing covariate values in contrast names

}
}
\details{
Based on the call and data set in the model object, matrices are contructed that allow to compute predictions over the grid of covariate values specified in \code{covset} and the factor variable specified in \code{treatment}. These matrices are then combined to define pairwise differences between the levels in treatment, over the grid of covariate values. The pairwise differences can be defined in the argument \code{treatcon}.
}
\value{
A list with items
\item{linfct}{a matrix of linear combinations of the model parameters, for plugging into the \code{linfct}-argument of \code{glht}, package multcomp}
\item{datadiff }{ a data set with names for the differences between treatment levels, the treatment levels contributing as minuend and subtrahend of the differences, and the covariate(s); rows of the data set are ordered as are the rows in the matrix \code{linfct}}
\item{Xpos}{model matrix for those predictions used as minuend to compute \code{linfct}}
\item{Xneg}{model matrix for those predictions used as subtrahend to compute \code{linfct}}
\item{datapos}{}
\item{dataneg}{}

}


\author{
Frank Schaarschmidt
}

\examples{

# Data set simulated acc. to Milliken and Johnson (2002), see ?pc for details

data(pc)

# A model including treatment interaction with the quadratic term

fitpc<-lm(yield ~ x + I(x^2) + additive + additive:I(x^2), data=pc)

cmpc <- cmiacov(fitpc, treatment="additive", 
 covset=list("x" = seq(from=1, to=9, by=1)), treatcon="Dunnett")
str(cmpc)


# treatment interaction with the linear term

fitpclin<-lm(yield ~ x + I(x^2) + additive + additive:x, data=pc)

cmpclin <- cmiacov(fitpc, treatment="additive",
 covset=list("x" = seq(from=1, to=9, by=1)), treatcon="Dunnett")
str(cmpclin)


# daphnids data from package drc
# no: counted no. of immobile daphnids
# mob: corresponding no. mobile daphnids

if(require("drc")){

data(daphnids, package="drc")

daphnids$mob <- daphnids$total - daphnids$no
daphnids$l10dose <- log10(daphnids$dose)

ggplot(daphnids, aes(y=no/total, x=l10dose)) + geom_point(aes(color=time))

fitda <- glm(cbind(no, mob) ~ time * l10dose, data=daphnids, family=quasibinomial())
anova(fitda, test="F")

cmda <- cmiacov(fitda, treatment="time", 
covset=list("l10dose" = seq(from=2, to=4, by=0.5)), treatcon="Dunnett")
str(cmda, max.level=1)
cmda$linfct

confint(glht(fitda, linfct=cmda$linfct))

}

}

\keyword{htest}
\concept{simultaneous confidence intervals}
