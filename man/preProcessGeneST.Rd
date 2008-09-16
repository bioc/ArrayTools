\name{preProcessGeneST}
\alias{preProcessGeneST}
\title{Proprocess genechip ST array}
\description{
  Proprocess genechip ST array by taking the log2 of the expression value.
}
\usage{
preProcessGeneST(object, offset = 1, rmControl = TRUE, output = FALSE, mydir = getwd())
}
\arguments{
  \item{object}{an \code{ExpressionSet}.}
  \item{offset}{The \code{offset} is added to the expression value to avoid log2(0) = -Inf.}
  \item{rmControl}{Setting \code{rmControl = TRUE} to remove control probes.}
  \item{output}{If \code{output = TRUE}, it will output the preprocessed data in the specified direcotry 
    from the \code{mydir} argument.}
  \item{mydir}{specifed directory to contain the output}
} 

\value{
  an \code{ExpressionSet}
}
\author{Xiwei Wu, Arthur Li}
\examples{
   data(eSetExample)
   processedData <- preProcessGeneST(eSetExample)  
}
