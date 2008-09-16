\name{createExpressionSet}
\alias{createExpressionSet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Creating an ExpressionSet}
\description{
   Create an \code{ExpressionSet} based on phenotype data and expression data
}
\usage{
createExpressionSet(pData, exprs, ...)
}

\arguments{
  \item{pData}{a data frame contains the phenotype data}
  \item{exprs}{a data frame contains the expression data}
  \item{\dots}{addtional arguments passed to \code{new("ExpressionSet", exprs, phenoData, ...)} if needed }
}
\value{
  an \code{ExpressionSet}
}
\author{Xiwei Wu, Arthur Li }
\seealso{\link[Biobase]{ExpressionSet}}
\references{
  Bioconductor: Open software development for computational biology and
  bioinformatics R. Gentleman, V. J. Carey, D. M. Bates, B.Bolstad, M.
  Dettling, S. Dudoit, B. Ellis, L. Gautier, Y. Ge, and others 2004,
  Genome Biology, Vol. 5, R80
}
\examples{
  data(pDataExample)
  data(exprsExample)
  eSet <- createExpressionSet (pDataExample, exprsExample, 
      annotation = "hugene10st")
}
