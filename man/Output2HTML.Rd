\name{Output2HTML}
\alias{Output2HTML}
\title{Creating HTML file for regressResult or interactionResult class}
\description{
  Creating HTML file for regressResult or interactionResult class
}
\usage{
Output2HTML(object, ...)
}
\arguments{
  \item{object}{an regressResult or interactionResult class}
  \item{\dots}{you can specify the directory to store the result
     by using the \code{mydir} argument.  The default value of 
	 \code{mydir} is the current working directory}
}

\value{
  creating an HTML file
}
\author{Xiwei Wu, Arthur Li}
\examples{
data(eSetExample)
design<- new("designMatrix", target=pData(eSetExample), covariates = "Treatment")
contrast<- new("contrastMatrix", design.matrix = design, 
    compare1 = "Treated", compare2 = "Control")
result<- regress(eSetExample, contrast)
sigResult<- selectSigGene(result, fc.value=log2(2)) 
\dontrun{Output2HTML(sigResult)}
}
