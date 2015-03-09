\name{statintcov-package}
\alias{statintcov-package}
\alias{statintcov}
\docType{package}
\title{Multiple comparisons of treatments in presence of treatment-covariate-interactions}

\description{For the case  of treatment-covariate interaction in linear models, this package provides functions to preform multiple comparisons among several treatments over a grid of (many) pre-specified covariate values. It thus may be used as an alternative to confidence bands for differences of multiple treatments. For the special case of (general) linear models with one treatment factor and one covariate, simultaneous confidence intervals for between-treatment-contrasts over a prespecified grid of covariate values can be computed. Such confidence intervals can be computed for differences of regression lines and (approximately) ratios of regression lines. For use in linear or generalized linear models, matrices of linear combinations for multiple pairwise between-treatment comparisons can be computed that can be used for plug-in into the linfct argument of glht, package multcomp.
}
\details{
\tabular{ll}{
Package: \tab statintcov\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-03-09\cr
License: \tab What license is it under?\cr
}

For the special case of one treatment factor and interaction to one covariate, the function \code{scitreatcov}
has a simple interface allowing the specification of response, tretament factor and covariate, the treatment contrast (for differences) of interest and a vector of covariate values. Internally, a suitably parameterized linear model is fit, a corresponding contrast matrix is contructed. Returned are simultaneous confidence intervals for that grid, internally using the function \code{confint.glht} from package \code{multcomp}. The function \code{sciratiotreatcov} allows to perform the same type of analysis for ratios of regression lines, relying on the approximate (!) methods available in function \code{gsci.ratio} from package \code{mratios}. 

For the slightly more general case of (generalized) linear models, that may involve more than one covariate or linear and quadratic dependency on the covariate, the function \code{cmiacov} is provided. For a given fitted model object, a grid of covariate values, and a matrix defining pairwise differences of interest among the levels of a treatment factor, a matrix of coefficients is computed that can be plugged into the argument \code{linfct} of the \code{glht} function in package \code{multcomp}.

}
\author{
Frank Schaarschmidt <schaarschmidt@biostat.uni-hannover.de>
}

\examples{

###############################################
# All pairwise comparisons in a linear model  #
# with baseline as covariate                  #
###############################################

# Loading the data set

if(require("MASS")){

# load data from package MASS:

data(anorexia, package="MASS")

ggplot(anorexia, aes(y=Postwt, x=Prewt)) + geom_point() +
facet_wrap(~Treat) 

fitan <- lm(Postwt ~ Treat*Prewt, data=anorexia)
anova(fitan)

# Comparison to control (cont, which is the second level of factor Treat)
# does any of the two treatments leads to higher postweights compared to control
# at covariate values 75,76,...,95

dscian <- scitreatcov(response="Postwt", treatment="Treat",
covariate="Prewt", data=anorexia, covset=seq(from=75, to=95, by=5),
treatcon="Dunnett", conf.level=0.95, alternative="greater", base=2)

str(dscian, max.level=1)
str(dscian$sci)

 ggplot(dscian$sci, aes(x=covariate, y=Estimate, ymin=lwr, ymax=upr)) + 
 geom_errorbar(width=1) + geom_line() + geom_point(shape=15) +
 facet_grid(.~comparison) + xlab("Preweight") + 
 ylab("Difference in expected postweight") + geom_hline(yintercept=0)

}


###############################################
# Comparisons to control with interaction to  #
# a quadratic regression term                 #
###############################################

data(pc)

ggplot(pc, aes(y=yield, x=x)) + geom_point() + facet_wrap(~additive) 

# Model selected by Milliken and Johnson:

fitpc<-lm(yield ~ x + I(x^2) + additive + additive:I(x^2), data=pc)

anova(fitpc)

# Since this model contains a quadratic term, it can onyl be handled via cmiacov
# 1) computing the contrast matrix, suitable for parameterization in fitpc
# the only covariate in the data is x, although it appears in several terms
# comparsisons to control

cmpc <- cmiacov(fitpc, treatment="additive",
 covset=list("x" = seq(from=1, to=10, by=1)), treatcon="Dunnett")

str(cmpc, max.level=1)
str(cmpc$linfct)


# 2) pass the contrast matrix cmpc$linfct to the linfct argument of glht

scipc <- confint(glht(fitpc, linfct=cmpc$linfct))

# 3) for plotting: write the comparsions and covariate values in a data.frame
# together with the confidence interavls from glht

dscipc <- cbind(cmpc$datadiff, scipc$confint)

ggplot(dscipc, aes(x=x, y=Estimate, ymin=lwr, ymax=upr)) +
 geom_errorbar(width=1) + geom_line() + geom_point(shape=15) + 
 facet_grid(.~comparison) + xlab("Compound x") + 
 ylab("Difference in expected yield") + geom_hline(yintercept=0)



######################################################
# All pairwise comparisons in a binomial generalized #
# linear model with logit link                       #
######################################################

# to load the data set

if(require("drc")){

data(selenium, package="drc")

se <- subset(selenium, conc!=0)
se$typef <- factor(se$type)
levels(se$typef) <- c("Selenate","Selenite","Selenomethionine","Selenocysteine")
se$l10conc <- log10(se$conc)

se$mortality <- se$dead/se$total

ggplot(se, aes(y=mortality, x=conc)) + geom_point() +
facet_wrap( ~typef) + scale_x_log10()

fitse<-glm(cbind(dead, total-dead) ~ typef*l10conc, data=se, family=binomial())
anova(fitse, test="Chisq")

# the binomial assumption, or logit-log linearity may not be adequat
# using the quasibinomial assumption should be more reliable

# 1) computing the contrast matrix, suitable for parameterization in fitse
# all pairwise comparisons ("Tukey")
# a sequence of 7 points over the observed range of the covariate

range(se$l10conc)

cmse <- cmiacov(fitse, treatment="typef", 
covset=list("l10conc" = seq(from=min(se$l10conc),
to=max(se$l10conc), length.out=7)), treatcon="Tukey")

str(cmse, max.level=1)
str(cmse$linfct)


# 2) pass the contrast matrix cman$linfct to the linfct argument of glht

scise <- confint(glht(fitse, linfct=cmse$linfct))

# 3) for plotting: write the comparisons and covariate values in a data.frame
# together with the confidence intervals from glht

dscise <- cbind(cmse$datadiff, scise$confint)

ggplot(dscise, aes(x=l10conc, y=Estimate, ymin=lwr, ymax=upr)) +
geom_errorbar(width=0.1) + geom_line() + geom_point(shape=15) +
facet_wrap(~comparison) + xlab("Log10 Concentration") +
ylab("Difference (logit-scale)") + geom_hline(yintercept=0)


}


}
