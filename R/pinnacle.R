# Gabriel Hoffman
# May 25, 2022


# TODO 
# - add gene properties as covariates
# 	- pass as data to lm_each_eclairs()
# 	- specify formula 
# - Compute residuals from fit excluding coef
# - use EB tstats or other


#' Evaluate Gene Set Analysis with Pinnacle
#'
#' Perform Gene Set Analysis (GSA) by comparing t-statistics from a given gene set to genome-wide t-statistics, while accounting for co-expression structure
#'
#' @param fit regression model fit by \code{limma::lmFit()} or \code{variancePartition::dream()}
#' @param coef indicate coefficient or contrast to be extracted from \code{fit} using \code{topTable}
#' @param geneSets \code{GeneSetCollection} from \link{GSEABase}
#' @param data \code{data.frame} storing properties of each gene, with rownames being gene names
#' @param setSize array of two elements specifying the min and max number of genes allowed in a gene set. Only gene sets satisfying these criteria are retained 
#' @param formula formula specifying covariates in regression using t-statistics as response
#' @param quiet suppresss messages 
#' @param ... other arguments passed to \code{lm_each_eclairs()} and then \link{lm}
#'
#' @return \code{data.frame} with results for each gene set
#' 
##' @details
#' 
##' @examples
#'
#' @import variancePartition limma GSEABase Rdpack
#' @importFrom decorrelate eclairs lm_each_eclairs
#' @importFrom stats as.formula model.matrix p.adjust
#' @importFrom methods is
#' @importFrom Matrix colSums
#' 
#' @export
pinnacle = function( fit, coef, geneSets, data, formula = ~ 1, setSize=c(10, 5000), quiet = FALSE,...){ # ecl=NULL,

	# Checks
	########

	# check that geneset is the correct type
	if( !is(geneSets, "GeneSetCollection") ){
		stop("geneSets must be of type 'GeneSetCollection'")
	}

	# check that coef is a string
	if( ! is(coef, "character") ){
		stop("coef must be a string")
	}

	# check that only 1 coef is given
	if( length(coef) != 1 ){
		stop("coef must have exactly one entry")
	}

	# check fit
	if( ! is(fit, 'MArrayLM') ){
		stop("fit must be result from lmFit or dream")
	}

	# check setSize
	if( is(setSize, "array") | length(setSize) != 2){
		stop("setSize must be an array of length 2")
	}

	# Extract tstats
	################

	# extract t-statistics
	tstat = topTable(fit, coef=coef, number=Inf, sort.by="none")[,'t']
	names(tstat) = rownames(fit)

	# Gene Properties
	#################

	if( missing(data) ){
		data = data.frame(Amean = fit$Amean)
	}else{
		data = merge(data, data.frame(Amean = fit$Amean), by="row.names")
	}

	# keep only genes that are in tstat
	data = data[rownames(data) %in% names(tstat),,drop=FALSE]

	# sort genes the same order as in tstat
	data = data[match(names(tstat), rownames(data)),,drop=FALSE]
	# identical(names(tstat), rownames(data))

	form2 = paste('tstat', paste(as.character(formula), collapse=' ')) 
	form2 = as.formula(form2)

	dsgn = model.matrix(form2, data)

	# Process geneSets
	##################

	if( ! quiet ) cat(" Converting gene sets...\n")
	
	# encode geneSets on sparseMatrix X_genesets
	X_genesets = pinnacle::recodeToSparseMatrix( geneSets, names(tstat) )

	# keep only genes that are in tstat
	# X_genesets = X_genesets[rownames(X_genesets) %in% names(tstat),,drop=FALSE]

	# sort genes in X_genesets the same order as in tstat
	# X_genesets = X_genesets[match(names(tstat), rownames(X_genesets)),,drop=FALSE]
	if( ! identical(names(tstat), rownames(X_genesets)) ){
		stop("order is wrong")
	}

	# keep gene sets with number of entries passing setSize cutoffs
	ngenes = colSums(X_genesets)
	idx = (ngenes >= min(setSize)) & (ngenes <= max(setSize)) 
	X_genesets = X_genesets[,idx,drop=FALSE]

	# Perform eclairs decomposition
	###############################

	if( is.null(ecl) ){
		if( ! quiet ) cat(" Performing eclairs decomposition...\n")

		Y_resid = t(residuals(fit))
		ecl = eclairs( Y_resid, compute="correlation")
	}

	# Fit regressions
	#################

	if( ! quiet ) cat(" Fitting regression models...\n")

	# fit pinnacle model using eclairs
	fit_pin = lm_each_eclairs(
					formula 		= form2, 
					data 			= data, 
					X 				= X_genesets, 
					Sigma.eclairs 	= ecl,
					...)

	# post-processing
	fit_pin$FDR = p.adjust(fit_pin$pvalue)

	# return result
	fit_pin
}