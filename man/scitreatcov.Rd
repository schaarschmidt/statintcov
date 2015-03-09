\name{scitreatcov}
\alias{scitreatcov}
\alias{sciratiotreatcov}

\title{
Simultaneous treatment comparisons for a set of pre-specifed covariate values in a linear model.
}

\description{
For a linear model with one factor variable, one covariate, and a factor-covariate interaction, simultaneous confidence intervals for user-defined treatment contrasts are computed for a set of pre-specified covariate values. \code{scitreatcov} compares treatments in terms of differences, code{sciratiotreatcov} compares treatments in terms of ratios. 
}

\usage{
scitreatcov(response, treatment, covariate, data,
 covset = NULL, nocov = 10, treatcon = "Dunnett",
 conf.level = 0.95, alternative = "two.sided", ...)

sciratiotreatcov(response, treatment, covariate, data,
 covset = NULL, nocov = 10, treatcon = "Dunnett",
 conf.level = 0.95, alternative = "two.sided", ...)
}

\arguments{
  \item{response}{a single character string, naming the variable in the data set, containg the numeric response }
  \item{treatment}{a single character string, naming the variable in the data set, containg the factor variable}
  \item{covariate}{a single character string, naming the variable in the data set, containg the numeric covariate}
  \item{data}{the data set (a \code{data.frame} contaning the three variables above)}
  \item{covset}{a numeric vector, containing the values of the covariate for which treatment comparisons should be performed}
  \item{nocov}{a single integer, defines the number of covariate values starting from \code{min(covariate)} ending at \code{max(covariate)}; ignored, if covset is specified}
  \item{treatcon}{either a matrix with numeric entries, or a single character string; if a matrix: number of columns must be = number of treatment levels in the factor variable, rows should specify contrasts between treatments. If a character string: the character string is passed to the \code{type}-argument in \code{contrMat} (package \code{multcomp}), see \code{?contrMat} for detail, for \code{sciratiotreatcov}, see \code{contrMatRatio}}
  \item{conf.level}{a single numeric value ]0;1[, the simultaneous confidence level }
  \item{alternative}{one of "two.sided", "less" (upper limits only), "greater" (lower limits only)}
  \item{\dots}{arguments to be passed to other functions. Currently, in both functions:
 	\code{base}:a single integer indicating the control level, if treatcon is a character stringf, passed to contrMat,
	\code{digits}:a single integer, defining the number of visible non-zero digits when rounding covariate values in contrast names.
In \code{sciratiotreatcov}:
	\code{adjusted}: logical, to be passed to \code{gsci.ratio}.All other arguments are ignored.

}
}
\details{
The functions fit a linear model (using \code{lm}) of the form \code{response ~ 0 + treatment + treatment:covariate}. It constructs a suitable contrast matrix to compute simultaneous confidence intervals (sci) for multiple treatment comparisons (as defined in \code{treatcon}) for the covariate vales defined by \code{covset} or \code{nocov}. \code{scitreatcov} computes sci for differences using function \code{glht} and \code{confint.glht}, package \code{multcomp}. \code{sciratiotreatcov} computes sci for differences using function \code{gsci.ratio} in package \code{mratios}. 

}
\value{
a list with elements:
\item{sci}{a \code{data.frame} with columns:
 \code{comparison} (name of between-treatment comparisons);
 \code{covariate} (covariate values);
 \code{estimate} (point estimate of the contrast);
 \code{lwr, upr} lower and upper confidence limits}
\item{model.fit}{the linear model fit}
\item{glht}{the object returned by \code{glht}, package\code{multcomp}}

\item{treatmentcontrasts}{the matrix defining the contrasts between treatments}
\item{covset}{numeric vector with the covariate values used}
\item{alternative}{as input}
\item{conf.level}{as input}
}

\author{
Frank Schaarschmidt
}

\examples{

if(require("MASS")){

data(anorexia, package="MASS")

anova(lm(Postwt ~ Treat*Prewt, data=anorexia))

dscian <- scitreatcov(response="Postwt", treatment="Treat",
covariate="Prewt", data=anorexia, covset=seq(from=75, to=95, by=5),
treatcon="Tukey", conf.level=0.95)

str(dscian)

 ggplot(dscian$sci, aes(x=covariate, y=Estimate, ymin=lwr, ymax=upr)) +
 geom_errorbar(width=1) + geom_line() +  geom_point(shape=15) +
 facet_grid(.~comparison) + xlab("Preweight") + 
 ylab("Difference in expected postweight") + geom_hline(yintercept=0)

}

# for ratios, :

if(require("MASS")){

data(anorexia, package="MASS")

dscianr <- sciratiotreatcov(response="Postwt", treatment="Treat",
covariate="Prewt", data=anorexia, covset=seq(from=75, to=95, by=5),
treatcon="Tukey", conf.level=0.95)

str(dscianr, max.level=1)
str(dscianr$sci)

ggplot(dscianr$sci, aes(x=covariate, y=estimate, ymin=lower, ymax=upper)) +
geom_errorbar(width=1) + geom_line() + geom_point(shape=15) + 
facet_grid(.~comparison) + xlab("Preweight") + 
ylab("Difference in expected postweight") + geom_hline(yintercept=1)

}

}

\keyword{htest}
\concept{simultaneous confidence intervals}
