\name{Sort}
\alias{Sort}
\title{Sort a regressionResult or an interactionResult}
\description{
  Sort a regressionResult or an interactionResult based on p-value, fold-change,
  or F statistics
}
\usage{
Sort(x, ...)
}
\arguments{
  \item{x}{a regressResult or an interactionResult class}
  \item{\dots}{any other arguments.  See below... }
}
\section{Sort a regressResult or an interactionResult class}{
    \code{Sort(x, sorted.by = c("pValue", "log2Ratio", "F"), top=20)}
	
	\code{x} is a regressResult class or an interactionResult class.  
	\code{sorted.by} can be specified by using
	"pValue" (p value), "log2Ratio" (log2 of fold-change value) or "F" (F statistics).
	\code{top} is used to specified number of genes being printed
}

\value{
  if sorting a regressResult, returned value is a data frame
  if sorting a interactionResult, returned value is a list of data frames
}
\author{Xiwei Wu, Arthur Li}
\seealso{\code{\linkS4class{regressResult}} \code{\linkS4class{interactionResult}} }
\examples{
data(eSetExample)
design<- new("designMatrix", target=pData(eSetExample), covariates = "Treatment")
contrast<- new("contrastMatrix", design.matrix = design, 
    compare1 = "Treated", compare2 = "Control")
result<- regress(eSetExample, contrast)
Sort(result)
}