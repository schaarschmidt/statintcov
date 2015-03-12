cmiacov <-
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

if(is.character(treatcon)){
if(!(treatcon %in% c("Dunnett","Tukey"))){stop("Currently, only 'Dunnett' and 'Tukey' are allowed as default treatment contrasts")}
}

CMAT <- treatcon2cmat(treatcon=treatcon, base=dotargs$base, ntab=NTAB)

apc <- apply(X=CMAT, MARGIN=1, FUN=function(x){!all(x %in% c(-1,0,1))})
if(any(apc)){stop(paste("currently, contrast coefficients other than 0, 1, or -1 are not allowed, check cmat in row ", paste(which(apc), collapse=","), sep=""))}

nn0 <- apply(X=CMAT, MARGIN=1, FUN=function(x){length(which(x != 0))})
if(any(nn0>2)){stop(paste("Some rows of cmat contain more than 2 non-null coefficients, check cmat in row ", paste(which(nn0>2), collapse=","), sep=""))}

rsc <- rowSums(CMAT)
if(any(rsc-0 > sqrt(.Machine$double.eps))){cat("Some rows in cmat do not sum to 0.\n")}

if(is.null(rownames(CMAT))){NAMCMAT <- paste("c", 1:nrow(CMAT), sep="") }else{NAMCMAT <- rownames(CMAT)}

LEVWF <- levels(DAT[, treatment])

# construct an index vector pointing to those levels of treatment factor that
# are the positive and negative parts of the differences to estimate

ineg <- apply(X=CMAT, MARGIN=1, FUN=function(x){which(sign(x)==-1)})
ipos <- apply(X=CMAT, MARGIN=1, FUN=function(x){which(sign(x)==1)})

if(class(ineg) != "integer" | length(ineg)!=nrow(CMAT)){warning("At least one row of cmat contains more than oen negative coefficient")}
if(class(ipos) != "integer" | length(ipos)!=nrow(CMAT)){warning("At least one row of cmat contains more than one positive coefficient")}

# construct all combinations of the covariate values and the factor levels
# that contribute positively and negatively to the differences to estimate
# data sets, that could be used as newdata in predict.lm for the positive
# and negative part of the differences to estimate

efflistpos <- c(list(factor(LEVWF[ipos], levels=LEVWF)), covset)
names(efflistpos)[1]<-treatment

efflistneg <- c(list(factor(LEVWF[ineg], levels=LEVWF)), covset)
names(efflistneg)[1]<-treatment

ndpos <- expand.grid(efflistpos)
ndneg <- expand.grid(efflistneg)

# construct design matrices acc. to the postive and negative parts new data

Xpredpos <- model.matrix( RHS, data=ndpos)
Xpredneg <- model.matrix( RHS, data=ndneg)

# the difference of the positive and negative design matrices is the 
# contrast matrix for the differences in predicted lines given the actual
# parametrization in the model

cmatpreddiff <- Xpredpos - Xpredneg

#namfabb <- names(formals(abbreviate))
#wabbr <- namfabb %in% names(dotargs)
#if(any(wabbr)){abbargs<-formals(abbreviate)}

# construct appropriate names for the contrasts

if(is.null(dotargs$digits)){digits<-3}else{digits<-dotargs$digits}

if(length(COVNAMES)==1)
{
CMATNAMR <- paste(paste(ndpos[,treatment], " - ", ndneg[,treatment],sep=""), ", ", COVNAMES, " = ", signif(ndpos[,COVNAMES], digits=digits), sep="")
}

if(length(COVNAMES)>1){
CMATNAMCOV <- apply(X=ndpos, MARGIN=1, FUN=function(x){paste(signif(x, digits=digits), collapse=",")})
CMATNAMR <- paste(paste(ndpos[,treatment], " - ", ndneg[,treatment],sep=""), "(", paste(COVNAMES, collapse=","), ") = (", signif(ndpos[,COVNAMES], digits=digits),")", sep=", ")
}

row.names(cmatpreddiff)<-CMATNAMR

# create a data set containing with columns naming the positive, negative parts
# as well as the corresponding difference and the covariables, for later plotting

treatdiff <-  paste(ndpos[,treatment], " - ", ndneg[,treatment], sep="")

datadiff <- cbind(ndpos[,treatment], ndneg[treatment], treatdiff, ndpos[,COVNAMES])
names(datadiff) <- c(paste(treatment, c(".pos", ".neg") ,sep=""), "comparison", COVNAMES)


return(list(
linfct=cmatpreddiff,
datadiff=datadiff,
Xpos=Xpredpos,
Xneg=Xpredneg,
datapos=ndpos,
dataneg=ndneg))

}
