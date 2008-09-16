\name{selectSigGeneInt}
\alias{selectSigGeneInt}
\title{select differentially expressed genes from the interactionResult class}
\description{
  select differentially expressed genes based on p value and/or fold change 
  from the interactionResult class
}
\usage{
selectSigGeneInt(object, pGroup = 0.05, fcGroup = 0, pMain = 0.05, fcMain = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an \code{interactionResult} class}
  \item{pGroup}{the p value that used  to select significant genes at each level of the covariate}
  \item{fcGroup}{the fold change value that used  to select significant genes at each level of the covariate}
  \item{pMain}{the p values  that used to select significant genes among genes without any interaction effect}
  \item{fcMain}{the fold change values that used to select significant genes among genes without any interaction effect}
}
\value{
  an interactionResult
}
\author{Xiwei Wu, Arthur Li}
\examples{
data(eSetExample)
design.int<- new("designMatrix", target=pData(eSetExample), covariates = c("Treatment", "Group"),
    intIndex = c(1, 2))
contrast.int<- new("contrastMatrix", design.matrix = design.int, interaction=TRUE)
result.int<- regress(eSetExample, contrast.int)
sigResult.int <- selectSigGene(result.int)
intResult <- postInteraction(eSetExample, sigResult.int, mainVar ="Treatment",
   compare1 = "Treated", compare2 = "Control")
sigResultInt <- selectSigGeneInt(intResult)
}

