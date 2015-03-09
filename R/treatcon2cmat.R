treatcon2cmat <-
function(treatcon, base, ntab){

if(is.character(treatcon)){
if(length(treatcon)>1){warning("treatcon must be a single character string (or a numeric matrix), only the first character string is used."); TYPE <- treatcon[1]}else{TYPE <- treatcon}
if(is.null(base)){CMAT <- contrMat(n=ntab, type=TYPE)}else{
CMAT <- contrMat(n=ntab, type=TYPE, base=base)}
}else{
if(is.matrix(treatcon) & (is.numeric(treatcon)|is.integer(treatcon))){
CMAT <- treatcon
if(ncol(CMAT)!=length(ntab)){stop("matrix specified in treatcon must have as many columns, as there are factor levels in the treatment variabel!")}

if(is.null(rownames(CMAT))){rownames(CMAT) <- paste("c", 1:nrow(CMAT), sep="") }

}else{stop("treatcon must be a numeric matrix (or a single character string)")}
}
return(CMAT)
}
