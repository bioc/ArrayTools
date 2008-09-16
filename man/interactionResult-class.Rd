\name{interactionResult-class}
\docType{class}
\alias{interactionResult-class}
\alias{interactionResult}
\alias{class:interactionResult}
\alias{adjustment,interactionResult-method}
\alias{getAdjP,interactionResult-method}
\alias{getAnnotation,interactionResult-method}
\alias{getContrast,interactionResult-method}
\alias{getF,interactionResult-method}
\alias{getFC,interactionResult-method}
\alias{getFCCutoff,interactionResult-method}
\alias{getFileName,interactionResult-method}
\alias{getFilterMethod,interactionResult-method}
\alias{getID,interactionResult-method}
\alias{getIndex,interactionResult-method}
\alias{getLength,interactionResult-method}
\alias{getNormalizationMethod,interactionResult-method}
\alias{getP,interactionResult-method}
\alias{getPCutoff,interactionResult-method}
\alias{Output2HTML,interactionResult-method}
\alias{regressionMethod,interactionResult-method}
\alias{selectSigGene,interactionResult-method}
\alias{show,interactionResult-method}
\alias{Sort,interactionResult-method}
\alias{summary,interactionResult-method}

\title{Class to Contain the Regression Result Based on An Interaction Model}
\description{Class to Contain the Regression Result Based on An Interaction Model.
   Interaction is a statistical term refering to a situation when the
   relationship between the outcome and the variable of the main interest differs at
   different levels of the extraneous variable}
\section{Creating Objects}{
    \code{interactionResult} object is generally created from the \code{postInteraction} function
    See \code{\link{postInteraction}}
}
\section{Object Components}{
   A list of four or more components.  Each component is a reggressResult class.
   The first component contains all the genes.
   The second component contains genes with the interaction effect
   The rest components contains genes with the interaction effect across different
       levels.  Each component contains the result for each level.
}
\section{Extends}{
Class \code{"\linkS4class{list}"}, from data part.
Class \code{"\linkS4class{vector}"}, by class "list", distance 2.
}
\section{Methods}{	 
	 \describe{
    \item{adjustment}{\code{signature(object = "regressResult")}}{: 
	    access the \code{adjustment} slot }
    \item{getAdjP}{\code{signature(object = "regressResult")}}{:  
	    access the \code{adjPVal} slot }
    \item{getAnnotation}{\code{signature(object = "regressResult")}}{:  
	    access the \code{annotation} slot }
    \item{getContrast}{\code{signature(object = "regressResult")}}{:  
	    access the \code{contrast} slot }
    \item{getF}{\code{signature(object = "regressResult")}}{:  
	    access the \code{FValue} slot}
    \item{getFC}{\code{signature(object = "regressResult")}}{:  
	    access the \code{foldChange} slot }
    \item{getFCCutoff}{\code{signature(object = "regressResult")}}{:  
	    access the \code{significantFCCutoff} slot}
    \item{getFileName}{\code{signature(object = "regressResult")}}{:  
	    access the \code{fileName} slot }
    \item{getFilterMethod}{\code{signature(object = "regressResult")}}{:  
	    access the \code{filterMethod} slot }
    \item{getID}{\code{signature(object = "regressResult")}}{:  
	    access the \code{ID} slot  }
    \item{getIndex}{\code{signature(object = "regressResult")}}{:  
	    access the \code{significantIndex} slot}
    \item{getNormalizationMethod}{\code{signature(object = "regressResult")}}{:  
	    access the \code{normalizationMethod} slot}
    \item{getP}{\code{signature(object = "regressResult")}}{:  
	    access the \code{pValue} slot }
    \item{getPCutoff}{\code{signature(object = "regressResult")}}{:  
	    access the \code{significantPvalueCutoff} slot }
    \item{Output2HTML}{\code{signature(object = "regressResult")}}{:  
	    create HTML file for sigificant genes in regressionResult}
    \item{regressionMethod}{\code{signature(object = "regressResult")}}{:  
	    access the \code{regressionMethod} slot}
    \item{selectSigGene}{\code{signature(object = "regressResult")}}{:  
	    select significant genes for \code{regressionResult} class}
    \item{show}{\code{signature(object = "regressResult")}}{:  
	    print \code{regressResult}}
    \item{Sort}{\code{signature(x = "regressResult")}}{:  
	    sort \code{regressResult}}
    \item{summary}{\code{signature(object = "regressResult")}}{:  
	    print the summary for \code{regressResult}}
	\item{getLength}{\code{signature(object = "interactionResult")}}{:  
        calculate the length of the interactionResult class}
	 }
}
\author{Xiwei Wu, Arthur Li}
\seealso{
	\code{\linkS4class{regressResult}}
}
\examples{
## Creating the interactionREsult takes a few steps:
data(eSetExample)
design.int<- new("designMatrix", target=pData(eSetExample), covariates = c("Treatment", "Group"),
    intIndex = c(1, 2))
contrast.int<- new("contrastMatrix", design.matrix = design.int, interaction=TRUE)
result.int<- regress(eSetExample, contrast.int)
sigResult.int <- selectSigGene(result.int)
intResult <- postInteraction(eSetExample, sigResult.int, mainVar ="Treatment",
   compare1 = "Treated", compare2 = "Control")
}
\keyword{classes}
