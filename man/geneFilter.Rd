\name{geneFilter}
\alias{geneFilter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{filter an ExpressionSet using different methods}
\description{
  Create a filtered 'ExpressionSet' based on background, range, or interquartile range 
}
\usage{
geneFilter(object, pct = 0.1, numChip = ceiling(ncol(exprs(object)) * pct), bg = 4, range = 0, iqrPct = 0, output = FALSE, mydir = getwd())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an \code{ExpressionSet}}
  \item{pct}{percentage}
  \item{numChip}{number of chips. If you would like to filter the \code{ExpressionSet} 
      based on at least 3 chips greater than 1 (bg=1), then set numChip = 3}
  \item{bg}{background value. If you would like to filter the \code{ExpressionSet}
      based on at least 3 chips greater than 1, then set bg=1}
  \item{range}{range = max value - min value of each gene}
  \item{iqrPct}{interquartile percentage}
  \item{output}{if output = TRUE, output filtered data in the sepecified directory }
  \item{mydir}{the directory containing the filtered data}
}
\details{
  There are three filtering methods. The User can use either one, two, or three. 
  1). At least a certain number of chips (numChip) are greater than a given background (bg)
  2). The range of the gene have to be greater than a given value (range)
  3). Calulating the interquatile range (IQR) of each gene to create an IQR vector. 
    Based on the given percentage (e.g. iqrPct=0.2), find the corresponding percentile. 
	If IQR is less than percentile, the gene will be filtered
}
\value{
  a filtered \code{ExpressionSet}
}
\author{Xiwei Wu, Arthur Li}
\examples{
data(eSetExample)
filtered <- geneFilter(eSetExample)
}