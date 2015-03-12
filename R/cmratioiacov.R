cmratioiacov <-
function(modelfit, treatment, treatcon="Dunnett", covset, ...)
{
 
# check the model object and extract data an formula

if(!(class(modelfit)[1] %in% c("lm","glm"))){stop("modelfit must be an object of class 'lm' or 'glm'")}

RHS <- as.formula(modelfit$call)[c(1,3)]
DAT <- modelfit$model

# check the name of the treatment variable

if(length(treatment)!=1 || !is.character(treatment)){ stop("treatment must be a single character string")}
if(!treatment %in% names(DAT)){ stop("treatment could not be found in data used to fit the model")}
if(!is.factor(DAT[,treatment])){ stop("the variable specified in treatment must be a factor")}

dotargs <- list(...)

# check the list of covariate values: correct format, can all covariate names be found in data?

if(!is.list(covset)){stop("covset must be a named list of numeric vectors")}

COVNAMES<-names(covset)
if(!all(COVNAMES %in% names(DAT))){ stop("Not all variables named in covset have been found in the data used for model fitting")}

# check contrast definition/matrix: as many columsn as factor levels, only 0,1,-1, no pooling contrasts 

NTAB <- table(DAT[, treatment])

CMAT <- treatcon2cmatratio(treatcon=treatcon, base=dotargs$base, ntab=NTAB)

if(is.character(treatcon)){
if(!(treatcon %in% c("Dunnett","Tukey"))){stop("currently, only 'Dunnett' and 'Tukey' are allowed as default treatment contrasts")}
}

# weighted means in denominator and numerator matrices?

apcden <- apply(X=CMAT$denC, MARGIN=1, FUN=function(x){!all(x %in% c(0,1))})
if(any(apcden)){stop(paste("currently, contrast coefficients other than 0 or 1 are not allowed, check denominator matrix (denC) in row ", paste(which(apcden), collapse=","), sep=""))}

apcnum <- apply(X=CMAT$numC, MARGIN=1, FUN=function(x){!all(x %in% c(0,1))})
if(any(apcnum)){stop(paste("currently, contrast coefficients other than 0 or 1 are not allowed, check numerator matrix (numC) in row ", paste(which(apcnum), collapse=","), sep=""))}


# more than 1 non-null entry row of denominator and numerator matrices?

nn0den <- apply(X=CMAT$denC, MARGIN=1, FUN=function(x){length(which(x != 0))})
if(any(nn0den>1)){stop(paste("Some rows of the denominator matrix (denC) contain more than 1 non-null coefficients, check denominator matrix (denC) in row ", paste(which(nn0den>1), collapse=","), sep=""))}

nn0num <- apply(X=CMAT$numC, MARGIN=1, FUN=function(x){length(which(x != 0))})
if(any(nn0num>1)){stop(paste("Some rows of the numerator matrix (numC) contain more than 1 non-null coefficients, check numerator matrix (numC) in row ", paste(which(nn0num>1), collapse=","), sep=""))}

#

LEVWF <- levels(DAT[, treatment])

# construct an index vector pointing to those levels of treatment factor that
# are the numerators of the ratios

iposnum  <- apply(X=CMAT$numC, MARGIN=1, FUN=function(x){which(sign(x)==1)})

if(class(iposnum) != "integer" | length(iposnum)!=nrow(CMAT$numC)){warning("At least one row of the numerator matrix contains more than one positive coefficient")}

# construct all combinations of the covariate values and the factor levels
# that contribute positively and negatively to the numerator or denominator,
# data sets, that could be used as newdata in predict.lm for the positive
# and negative part of the differences to estimate

efflistposnum <- c(list(factor(LEVWF[iposnum], levels=LEVWF)), covset)
names(efflistposnum)[1]<-treatment

ndposnum <- expand.grid(efflistposnum)

# construct design matrices acc. to the postive and negative parts new data

Xpredposnum <- model.matrix(RHS, data=ndposnum)

# the difference of the positive and negative design matrices is the 
# contrast matrix for the differences in predicted lines given the actual
# parametrization in the model


# the same procedure for the denominator

iposden  <- apply(X=CMAT$denC, MARGIN=1, FUN=function(x){which(sign(x)==1)})

if(class(iposden) != "integer" | length(iposden)!=nrow(CMAT$denC)){warning("At least one row of the denominator matrix contains more than one positive coefficient")}

efflistposden <- c(list(factor(LEVWF[iposden], levels=LEVWF)), covset)
names(efflistposden)[1]<-treatment

ndposden <- expand.grid(efflistposden)

# construct design matrices acc. to the postive and negative parts new data

Xpredposden <- model.matrix(RHS, data=ndposden)

# the difference of the positive and negative design matrices is the 
# contrast matrix for the differences in predicted lines given the actual
# parametrization in the model

# construct appropriate names for the contrasts

if(is.null(dotargs$digits)){digits<-3}else{digits<-dotargs$digits}

if(length(COVNAMES)==1)
{
CMATNAMR <- paste(paste(ndposnum[,treatment], " / ", ndposden[,treatment],sep=""), ", ", COVNAMES, " = ", signif(ndposnum[,COVNAMES], digits=digits), sep="")
}

if(length(COVNAMES)>1){
CMATNAMCOV <- apply(X=ndposnum, MARGIN=1, FUN=function(x){paste(signif(x, digits=digits), collapse=",")})
CMATNAMR <- paste(paste(ndposnum[,treatment], " - ", ndposden[,treatment],sep=""), "(", paste(COVNAMES, collapse=","), ") = (", signif(ndposnum[,COVNAMES], digits=digits),")", sep=", ")
}

row.names(Xpredposnum) <- row.names(Xpredposden) <- CMATNAMR

cmatRatio <- list(
numC = Xpredposnum,
denC = Xpredposden,
rnames = CMATNAMR)

# create a data set containing with columns naming the positive, negative parts
# as well as the corresponding difference and the covariables, for later plotting

treatratio<-  paste(ndposnum[,treatment], " / ", ndposden[,treatment], sep="")

dataratio <- cbind(ndposnum[,treatment], ndposden[treatment], treatratio, ndposnum[,COVNAMES])
names(dataratio) <- c(paste(treatment, c(".num", ".den") ,sep=""), "comparison", COVNAMES)


return(c(cmatRatio,
  list(
dataratio=dataratio,
Xnum=Xpredposnum,
Xden=Xpredposden,
datanum=ndposnum,
dataden=ndposden)))

}
