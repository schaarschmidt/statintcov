scitreatcov <-
function(response, treatment, covariate, data, covset=NULL, nocov=10, treatcon="Dunnett", conf.level=0.95, alternative="two.sided", ...)
{

dotargs<-list(...)

# check adequacy of data, response, treatment, covariate

if(!is.data.frame(data)){stop("data must be a data.frame")}

if(length(response)!=1 || !is.character(response)){stop("response must be a single character string")}
if(!response %in% names(data)){stop("response could not be found in data")}
if(!(is.numeric(data[,response]) | is.integer(data[,response]))){stop("response must identify a column of data, that contains numeric (or integer) values")}

if(length(covariate)!=1 || !is.character(covariate)){stop("covariate must be a single character string")}
if(!covariate %in% names(data)){stop("covariate could not be found in data")}
if(!(is.numeric(data[,covariate]) | is.integer(data[,covariate]))){stop("covariate must identify a column of data, that contains numeric (or integer) values")}

if(length(treatment)!=1 || !is.character(treatment)){stop("treatment must be a single character string")}
if(!treatment %in% names(data)){stop("treatment could not be found in data")}
if(!is.factor(data[,treatment])){data[,treatment] <- as.factor(data[,treatment]); warning("treatment has been coerced to be a factor variable")}

ntab <- table(data[,treatment])
ngroup <- length(ntab)
if(ngroup < 2){stop("There is only one level in the treatment factor.")}
if(ngroup == 2){warning("There are only two levels in the treatment factor.")}
if(any(ntab<2)){stop(paste("There are less than two observations in treatment level(s)", paste(names(ntab[which(ntab<2)]), collapse=","), sep=" "))}

# check adequate definition of set of covariate values

if(is.null(covset)){
if(is.null(nocov)){warning("Neither nocov nor covset is provided, nocov is set to 10"); NOCOV<-10}else{
if(length(nocov)!=1 || (!(is.numeric(nocov)|is.integer(nocov)))){
warning("nocov must be a single positive number"); NOCOV <- as.integer(nocov[1])}else{NOCOV<-nocov}
if(sign(NOCOV) != 1){ warning("Negative sign in nocov is ignored"); NOCOV <- abs(NOCOV)}
if(NOCOV<1){ warning("nocov smaller 1 is set to 1"); NOCOV<-1}}

rc<-range(data[, covariate])
COVSET <- seq(from=rc[1], to=rc[2], length.out=NOCOV)

}else{

if(!is.numeric(covset) & !is.integer(covset)){stop("covset must be a vector of numeric (or integer) values")}

rc <- range(data[, covariate])
NOCOV <- length(covset)
 
rcs <- range(covset)
if(rcs[1]<rc[1] | rcs[2]>rc[2]){cat("Note: range of covariate values in covset is not completely included in the range of observed covariate values in the data\n")}
COVSET <- sort(covset, decreasing=FALSE)

}

CMAT <- treatcon2cmat(treatcon=treatcon, base=dotargs$base, ntab=ntab)

# check dimensions of contrast matrix

NCOMP <- nrow(CMAT)

if((NOCOV*NCOMP)>1000){warning(paste("the total number of comparisons is too large, may be reduce to ", floor(1000/NCOMP), " or less covariate values or reduce to", floor(1000/NOCOV), "treatment contrasts.", sep=" "))}

NAMCMAT <- rownames(CMAT)

alternative <- match.arg(alternative, choices=c("two.sided","less","greater"))

if(is.null(conf.level) || (!is.numeric(conf.level) | length(conf.level)!=1)){stop("conf.level must be a single numeric value")}
if(conf.level>1|conf.level<0){stop("conf.level must be a number between 0 and 1")}

XMAT <- cbind(matrix(rep(1,NOCOV),ncol=1),  matrix(COVSET, ncol=1) )
#print(XMAT)

CMATCOV <- kronecker(XMAT, CMAT)
NAMCMATR <- factor(rep(NAMCMAT, times=NOCOV), levels=NAMCMAT)
COVR <- rep(COVSET, each=NCOMP)

if(is.null(dotargs$digits)){digits<-3}else{digits<-dotargs$digits}
NAMCMATCOV <- paste(NAMCMATR, ", ", covariate, "=", signif(COVR,digits=digits), sep="") 
rownames(CMATCOV)<-NAMCMATCOV

#print(NAMCMATCOV)

FORM <- as.formula(paste(response, " ~ 0 + ", treatment, paste("+", treatment,":", covariate, sep=""), sep=""))
FIT <- lm(formula=FORM, data=data)
COMP <- glht(model=FIT, linfct=CMATCOV, alternative=alternative, ...)
SCI <- confint(COMP, level=conf.level)

#print(str(SCI))
#print(NAMCMATR)
#print(COVR)

SCIDAT <- cbind(data.frame("comparison"=NAMCMATR, "covariate"=COVR), SCI$confint)


return(list(
sci=SCIDAT,
model.fit=FIT,
glht=COMP,
treatmentcontrasts=CMAT,
covset=COVSET,
alternative=alternative,
conf.level=conf.level
))


}
