\name{regressResult-class}
\docType{class}
\alias{regressResult-class}
\alias{regressResult}
\alias{class:regressResult}
\alias{adjustment,regressResult-method}
\alias{getAdjP,regressResult-method}
\alias{getAnnotation,regressResult-method}
\alias{getContrast,regressResult-method}
\alias{getF,regressResult-method}
\alias{getFC,regressResult-method}
\alias{getFCCutoff,regressResult-method}
\alias{getFileName,regressResult-method}
\alias{getFilterMethod,regressResult-method}
\alias{getID,regressResult-method}
\alias{getIndex,regressResult-method}
\alias{getNormalizationMethod,regressResult-method}
\alias{getP,regressResult-method}
\alias{getPCutoff,regressResult-method}
\alias{Output2HTML,regressResult-method}
\alias{regressionMethod,regressResult-method}
\alias{selectSigGene,regressResult-method}
\alias{show,regressResult-method}
\alias{Sort,regressResult-method}
\alias{summary,regressResult-method}

\title{Class to Contain the Regression Result}
\description{Class to Contain the Regression Result}
\section{Creating Objects}{
    \code{regressResult} object is generally created from the \code{regress} function
    See \code{\link{regress}}
}
\section{Slots}{
	 \describe{
    \item{\code{ID}:}{contains probe ID/gene ID}
    \item{\code{foldChange}:}{contains fold change value}
    \item{\code{FValue}:}{contains F statistics}
    \item{\code{pValue}:}{contains p value}
    \item{\code{adjPVal}:}{contains adjusted p value}
    \item{\code{contrast}:}{contains class \code{"contrastMatrix"}}
    \item{\code{regressionMethod}:}{contains regression method: "limma", "regression", 
	    or "permutation"}
    \item{\code{adjustment}:}{contains method for multiple comparison adjustment}
    \item{\code{significantIndex}:}{contains a logical index indicating sigificant genes}
    \item{\code{significantPvalueCutoff}:}{contains a cutoff p-value for choosing significant genes}
    \item{\code{significantFCCutoff}:}{contains a fold change cutoff value for choosing significant genes }
    \item{\code{fileName}:}{contains a file name for output purpose}
    \item{\code{annotation}:}{contains annotation}
    \item{\code{normalizationMethod}:}{contains normalization method - for output purpose}
    \item{\code{filterMethod}:}{contains filtered method - for output purpose}
  }
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
	 }
}
\author{Xiwei Wu, Arthur Li}
\examples{
data(eSetExample)
design<- new("designMatrix", target=pData(eSetExample), covariates = "Treatment")
contrast<- new("contrastMatrix", design.matrix = design, 
    compare1 = "Treated", compare2 = "Control")
result<- regress(eSetExample, contrast)
}
\keyword{classes}
