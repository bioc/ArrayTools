\name{createIngenuityFile}
\alias{createIngenuityFile}
\title{A Wrapper Function to Create Files for Ingenuity Analysis}
\description{
  A Wrapper Function to Create Files for Ingenuity Analysis
}
\usage{
createIngenuityFile(..., mydir = getwd(), eSet, filename = "IngenuityFile")
}
\arguments{
  \item{\dots}{a list of regressResult class }
  \item{mydir}{the directory where you would like to store the file}
  \item{eSet}{an ExpressionSet}
  \item{filename}{file name}
}
\details{
  This function enable to create the ingenuity upload file
  based on a list of regressResult
}
\value{
  create an Ingenuity upload file
}
\references{\code{http://www.ingenuity.com/}}
\author{Xiwei Wu, Arthur Li}
\examples{
data(eSetExample)
design<- new("designMatrix", target=pData(eSetExample), covariates = "Treatment")
contrast<- new("contrastMatrix", design.matrix = design, 
    compare1 = "Treated", compare2 = "Control")
result<- regress(eSetExample, contrast)
\dontrun{createIngenuityFile(result, eSet = eSetExample)}
}

