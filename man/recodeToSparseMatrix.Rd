% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/recodeGeneSet.R
\name{recodeToSparseMatrix}
\alias{recodeToSparseMatrix}
\title{Recode GeneSetCollection to sparseMatrix}
\usage{
recodeToSparseMatrix(gsc, geneNames)
}
\arguments{
\item{gsc}{GeneSetCollection object}

\item{geneNames}{order of gene names in result}
}
\description{
Recode GeneSetCollection to sparseMatrix with rows as genes, columns as gene sets, and a value of 1 indicating membership of gene i in gene set j.
}
