


#' Recode GeneSetCollection to data.table 
#'
#' Recode GeneSetCollection to data.table
#'
#' @param gsc GeneSetCollection object
#'
#' @rdname recodeToDataTable
#' @export
setGeneric("recodeToDataTable", function(gsc)
  standardGeneric("recodeToDataTable"))


#' Recode GeneSetCollection to data.table 
#'
#' Recode GeneSetCollection to data.table
#'
#' @param gsc GeneSetCollection object
#'
#' @rdname recodeToDataTable
#' @importFrom data.table data.table setindex
#' @importFrom GSEABase geneIds setName
#' @export
setMethod("recodeToDataTable", c("GeneSetCollection"),
  function( gsc ){

  # convert to data.table
  gsc.dt = lapply( gsc, function( geneSet ){
      # If gene set is empty, return NULL
      #   This drops the geneset from the data.frame 
      if( length(geneIds(geneSet)) == 0){
          res = NULL
      }else{
        # create data.table with the gene set name, and genes
          res = data.table( setName = setName(geneSet), geneIds = geneIds(geneSet) )
      }
      res
      })
  gsc.dt = do.call(rbind, gsc.dt)

  # Create Index keyed on geneIds for binary search
  setindex( gsc.dt, 'geneIds')
})


#' Recode GeneSetCollection to list used by limma
#'
#' Recode GeneSetCollection to list used by limma
#'
#' @param gsc GeneSetCollection object
#'
#' @rdname recodeToList
#' @export
setGeneric("recodeToList", function(gsc)
  standardGeneric("recodeToList"))

#' Recode GeneSetCollection to list used by limma
#'
#' Recode GeneSetCollection to list used by limma
#'
#' @param gsc GeneSetCollection object
#'
#' @rdname recodeToList-GeneSetCollection
#' @importFrom GSEABase geneIds 
#' @export
setMethod("recodeToList", c("GeneSetCollection"),
  function( gsc ){

    # convert to list  
  gsList = lapply( gsc, geneIds) 
  names(gsList) = names(gsc)

  gsList
})


#' Recode GeneSetCollection to sparseMatrix 
#'
#' Recode GeneSetCollection to sparseMatrix with rows as genes, columns as gene sets, and a value of 1 indicating membership of gene i in gene set j.
#'
#' @param gsc GeneSetCollection object
#' @param geneNames order of gene names in result 
#'
#' @rdname recodeToSparseMatrix
#' @export
setGeneric("recodeToSparseMatrix", function(gsc, geneNames)
  standardGeneric("recodeToSparseMatrix"))


#' Recode GeneSetCollection to sparseMatrix 
#'
#' Recode GeneSetCollection to sparseMatrix
#'
#' @param gsc GeneSetCollection object
#' @param geneNames order of gene names in result 
#'
#' @rdname recodeToSparseMatrix-GeneSetCollection
#' @importFrom Matrix sparseMatrix
#' @export
setMethod("recodeToSparseMatrix", c("GeneSetCollection"),
  function(gsc, geneNames){

  	# convert to data.table
  	gsc.dt = recodeToDataTable( gsc )

    if( missing(geneNames) ){
      # if geneNames is not specified, use unique set from Gene set collection
      geneNames = gsc.dt[,unique(geneIds)]
    }else{
      # if specified only retain genes in geneNames
      gsc.dt = gsc.dt[geneIds %in% geneNames,]
    }

  	# For each gene get the index based on all genes
  	#	and for each geneset, get the index based on all genesets
  	# Converting strings to factors does this automatically
  	# use the gene and geneset indeces to set a value in a 
  	#	as sparseMatrix to 1
  	# Set rows (genes) and columns (genesets)

  	# pass R check
  	geneIds = setName = NULL
    `:=` = function(...) NULL

  	# convert to factor
    # use specified gene and set names as levels
    # this allows  
  	gsc.dt[,geneIds := factor(geneIds, levels=geneNames)]
  	gsc.dt[,setName := factor(setName, levels=names(gsc))]

  	# set indeces as integer and set corresponding value to 1
  	i = as.integer( gsc.dt$geneIds )
  	j = as.integer( gsc.dt$setName )
  	value = rep(1, length(i))

  	sparseMatrix( 	
            i       = i,
  					j        = j,
  					x        = value,
            dims     = c(nlevels(gsc.dt$geneId), nlevels(gsc.dt$setName) ),
  					dimnames = list( levels(gsc.dt$geneId), levels(gsc.dt$setName) ))
})


#' Get pairwise Jaccard similarity between genesets
#'
#' Get pairwise Jaccard similarity between genesets
#'
#' @param x collection of gene sets
#'
#' @rdname getJaccardMatrix
#' @export
setGeneric("getJaccardMatrix", function( x)
  standardGeneric("getJaccardMatrix"))



#' Get pairwise Jaccard similarity between genesets
#'
#' Get pairwise Jaccard similarity between genesets
#'
#' @param x collection of gene sets
#'
#' @rdname getJaccardMatrix-GeneSetCollection
#' @importFrom locStra jaccardMatrix
#' @export
setMethod("getJaccardMatrix", c("GeneSetCollection"),
  function( x ){

  # convert to sparseMatrix
  M = recodeToSparseMatrix( x )

  getJaccardMatrix( M )
})





#' Get pairwise Jaccard similarity between genesets
#'
#' Get pairwise Jaccard similarity between genesets
#'
#' @param x collection of gene sets
#'
#' @rdname getJaccardMatrix-sparseMatrix
#' @importFrom locStra jaccardMatrix
#' @export
setMethod("getJaccardMatrix", c("sparseMatrix"),
  function( x ){

  # compute Jaccard distances from sparseMatrix
  C_sim = jaccardMatrix( x )

  rownames(C_sim) = colnames( x )
  colnames(C_sim) = colnames(x)

  C_sim
})




#' Get pairwise Jaccard similarity between genesets
#'
#' Get pairwise Jaccard similarity between genesets
#'
#' @param x collection of gene sets
#'
#' @rdname getJaccardMatrix
#' @export
setGeneric("getJaccardMatrix", function( x)
  standardGeneric("getJaccardMatrix"))



#' Get pairwise Jaccard similarity between genesets
#'
#' Get pairwise Jaccard similarity between genesets
#'
#' @param x collection of gene sets
#'
#' @rdname getJaccardMatrix-GeneSetCollection
#' @importFrom locStra jaccardMatrix
#' @export
setMethod("getJaccardMatrix", c("GeneSetCollection"),
  function( x ){

  # convert to sparseMatrix
  M = recodeToSparseMatrix( x )

  getJaccardMatrix( M )
})





#' Get pairwise Jaccard similarity between genesets
#'
#' Get pairwise Jaccard similarity between genesets
#'
#' @param x collection of gene sets
#'
#' @rdname getJaccardMatrix-sparseMatrix
#' @importFrom locStra jaccardMatrix
#' @export
setMethod("getJaccardMatrix", c("sparseMatrix"),
  function( x ){

  # compute Jaccard distances from sparseMatrix
  C_sim = jaccardMatrix( x )

  rownames(C_sim) = colnames( x )
  colnames(C_sim) = colnames(x)

  C_sim
})




#' Plot gene sets in two dimensions
#'
#' Plot gene sets in two dimensions based on Jaccard distance
#'
#' @param gsc GeneSetCollection object
#' @param res data.frame with score for each gene set
#' @param value indicate which columns of res to use to color points
#' @param zlim limits of colors.  c(0,NA) indicates a range from 0 to the largest observed value
#'
#' @import ggplot2
#' @importFrom stats cmdscale
#' @importFrom ggrepel geom_text_repel
#' @export
plot_mds = function( gsc, res, value="-log10(p.shift)", zlim=c(0, NA) ){

  # compte Jaccard matrix
  M_s = getJaccardMatrix( gsc )

  # Create distance matrix and perform MDS
  loc = cmdscale( as.dist(1 - M_s) )

  df = data.frame(  Geneset = rownames(loc),
            A1    = loc[,1],
            A2    = loc[,2], 
            stringsAsFactors=FALSE)

  res = res[res$Geneset %in% df$Geneset,]
  df = merge(df, res, by="Geneset")

  fig = ggplot(df, aes_string('A1', 'A2', label='Geneset', color=value)) + geom_text_repel() + theme_bw()  + geom_point(size=5) + xlab("Coordinate 1") + ylab("Coordinate 2") + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())#, axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())

  fig = fig + scale_color_gradient(low="white", high="red", limit=zlim)

  fig 
}

#' Plot gene sets in two dimensions
#'
#' Plot gene sets in two dimensions based on Jaccard distance
#'
#' @param gsc GeneSetCollection object
#' @param df_GeneMap result from getGeneMapping() 
#' @param coef indicate contrast column in df_GeneMap
#' @param bicluster apply clustering to rows and columns before plotting
#' @param setOverlap keep genes present in at least this many genes sets
#' @param zlim scalar. specify max value for color scale
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom Matrix rowSums
#' @importFrom stats dist hclust as.dist
#' @export
heatplot = function( gsc, df_GeneMap, coef, setOverlap=2, bicluster=TRUE, zlim=NULL){

  # convert to sparseMatrix
  M = recodeToSparseMatrix( gsc )

  # gene must be present in more than 1 set
  M = M[rowSums(M) > setOverlap,]

  if( nrow(M) < 1){
    warning( "No genes were retained" )
    return(invisible())
  }
  if( ncol(M) < 1){
    warning( "No genesets were retained" )
    return(invisible())
  }

  # pass R CMD check
  Gene = Geneset = statistic = NA

  # new data.frame for differential expression results
  df_de = data.frame( Gene  = df_GeneMap$targetGene,
              stat  = df_GeneMap[[coef]], 
              stringsAsFactors=FALSE)

  # subset to Genes in M
  df_de = df_de[df_de$Gene %in% rownames(M),]

  # subset M to genes in df_de
  M = M[rownames(M) %in% df_de$Gene,]

  # reorder df_de to match M
  df_de = df_de[match(rownames(M), df_de$Gene),]

  # put DE statistic in each gene set for each gene
  M_DE = as.matrix(M * df_de$stat )

  # biclustering
  if( bicluster ){
    hcl1 = hclust(dist(M_DE), method="ward.D2")
    hcl2 = hclust(dist(t(M_DE)), method="ward.D2")
    M_DE = M_DE[hcl1$order,hcl2$order]
  }

  # melt to get tall dataset
  df_melt = melt( M_DE )
  colnames(df_melt) = c("Gene", "Geneset", "statistic")

  if( is.null(zlim) ){
    zlim = max(abs(df_melt$statistic))
  }

  ggplot(df_melt, aes(Gene, Geneset, fill=statistic)) + geom_tile() + theme_bw(8) + theme(aspect.ratio=ncol(M) / nrow(M), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title = element_blank(), axis.text.x=element_text(angle=45, hjust=1)) + scale_fill_gradient2("Statistic", low="navy", high="darkred", mid="white", midpoint=0, limit=c(-zlim, zlim))
}






