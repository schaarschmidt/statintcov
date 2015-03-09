treatcon2cmatratio <-
function(treatcon, base, ntab){

if(is.character(treatcon)){
if(length(treatcon)>1){warning("treatcon must be a single character string (or a list of two numeric matrices), only the first character string is used."); TYPE <- treatcon[1]}else{TYPE <- treatcon}
if(is.null(base)){CMAT <- contrMatRatio(n=ntab, type=TYPE)}else{
CMAT <- contrMatRatio(n=ntab, type=TYPE, base=base)}
}else{

if(is.list(treatcon) && (is.matrix(treatcon[[1]]) & is.matrix(treatcon[[2]])) & (is.numeric(treatcon[[1]])|is.integer(treatcon[[1]])) & (is.numeric(treatcon[[2]])|is.integer(treatcon[[2]]))){
CMAT <- treatcon

if(ncol(CMAT[[1]])!=length(ntab) | ncol(CMAT[[2]])!=length(ntab)){stop("matrices specified in treatcon must have as many columns, as there are factor levels in the treatment variabel!")}
if(is.null(names(CMAT))){names(CMAT)<-c("numC", "denC")}
if(nrow(CMAT[[1]]) != nrow(CMAT[[2]])){stop("the two matrices must have equal number of rows")}

if(is.null(CMAT$rnames) & (is.null( rownames(CMAT[[1]])) & is.null(rownames(CMAT[[2]])) )){
NAMCMAT <- paste("c", 1:nrow(CMAT), sep=""); rownames(CMAT[[1]])<-rownames(CMAT[[2]]) <- NAMCMAT; CMAT$rnames <- NAMCMAT}else{

if(is.null(CMAT$rnames)){
if(is.null( rownames(CMAT[[1]])) & !is.null(rownames(CMAT[[2]]))){NAMCMAT <- rownames(CMAT[[2]]); rownames(CMAT[[1]]) <- NAMCMAT; CMAT$rnames <- NAMCMAT}
if(is.null( rownames(CMAT[[2]])) & !is.null(rownames(CMAT[[1]]))){NAMCMAT <- rownames(CMAT[[1]]); rownames(CMAT[[2]]) <- NAMCMAT; CMAT$rnames <- NAMCMAT}
if(!is.null( rownames(CMAT[[2]])) & !is.null(rownames(CMAT[[1]]))){NAMCMAT <- paste(rownames(CMAT[[1]]), rownames(CMAT[[2]]), sep="/"); CMAT$rnames<-NAMCMAT}}}


}else{stop("treatcon must be a list of 2 numeric matrices (or a single character string)")}
}
return(CMAT)
}
