\name{createIndex}
\alias{createIndex}
\title{Creating an HTML index file }
\description{
  This HTML index file will link all the ouputted result, including
  Quality Assessment Report, differentially expressed genes, etc...
}
\usage{
createIndex(..., mydir = getwd(), index.file = "index.html", createHeader = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{regressionResults or interactionResult}
  \item{mydir}{the directory to contain the index file }
  \item{index.file}{name of the index file }
  \item{createHeader}{If want to want to create an Header, such as your name,
     company names, etc...}
}
\value{
  creating an HTML index-file in your directory
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

design.int<- new("designMatrix", target=pData(eSetExample), covariates = c("Treatment", "Group"),
    intIndex = c(1, 2))
contrast.int<- new("contrastMatrix", design.matrix = design.int, interaction=TRUE)
result.int<- regress(eSetExample, contrast.int)
sigResult.int <- selectSigGene(result.int)
intResult <- postInteraction(eSetExample, sigResult.int, mainVar ="Treatment",
   compare1 = "Treated", compare2 = "Control")
sigResultInt <- selectSigGeneInt(intResult)
\dontrun{Output2HTML(sigResultInt)}

\dontrun{createIndex(sigResult, sigResultInt, createHeader = c("Arthur Li", "COH"))}
}


