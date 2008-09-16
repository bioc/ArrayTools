\name{selectSigGene}
\alias{selectSigGene}
\title{select differentially expressed genes from the regressResult class}
\description{
  select differentially expressed genes based on p value and/or fold change 
  from the regressResult class
}
\usage{
selectSigGene(object, p.value = 0.05, fc.value = 0)
}
\arguments{
  \item{object}{an \code{regressResult} class}
  \item{p.value}{p value}
  \item{fc.value}{fold change cut-off value}
}

\value{
  an regressResult
}
\author{Xiwei Wu, Arthur Li}
\examples{
data(eSetExample)
design<- new("designMatrix", target=pData(eSetExample), covariates = "Treatment")
contrast<- new("contrastMatrix", design.matrix = design, 
    compare1 = "Treated", compare2 = "Control")
result<- regress(eSetExample, contrast)
sigResult<- selectSigGene(result, fc.value=log2(2)) 
}
