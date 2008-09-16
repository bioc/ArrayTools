\name{preProcess3prime}
\alias{preProcess3prime}
\title{A wrapper function to normalize the the 3 prime array}
\description{
  A wrapper function to normalize the 3 prime array by using either RMA or GCRMA method
}
\usage{
preProcess3prime(object, method = c("rma", "gcrma"), output = FALSE, mydir = getwd())
}
\arguments{
  \item{object}{an \code{AffyBatch}.}
  \item{method}{either \code{rma} or \code{gcrma} }
  \item{output}{If \code{output = TRUE}, it will output the preprocessed data in the specified direcotry 
    from the \code{mydir} argument}
  \item{mydir}{specifed directory to contain the output}
}

\value{
  an \code{ExpressionSet}
}
\author{Xiwei Wu, Arthur Li}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
    data(affybatch.example)
    eset <- preProcess3prime(affybatch.example)
}
