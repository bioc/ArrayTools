\name{designMatrix}
\docType{class}
\alias{designMatrix-class}
\alias{designMatrix}
\alias{class:designMatrix}

\alias{getCovariates,designMatrix-method}
\alias{getDesign,designMatrix-method}
\alias{getIntIndex,designMatrix-method}
\alias{getTarget,designMatrix-method}
\alias{initialize,designMatrix-method}
\alias{show,designMatrix-method}

\title{Class to Contain the Design Matrix that Used for Linear Regression}
\description{Class to Contain the Design Matrix that Used for Linear Regression}

\section{Creating Objects}{
    \code{new("designMatrix", ..., target, covariates, intIndex=0)}
	
	This create as design matrix class.  \code{target} is a data frame
	that contains chip and covaraite information, or experimental phenotypes 
	recorded in eSet and ExpressionSet-derived classes.  \code{covariates} is
	a list of 1-n covariates.  If \code{intIndex=0}, the interaction effect 
	is not considered; otherwise, use two integers to indicate which 
	covariates are considered for interaction effect. For example, 
	if \code{covariates <- c("estrogen","drug","time")} and you are considering 
	the interaction between "estrogen" and "time", then you would write 
	\code{intIndex=c(1,3)}
	 
}
\section{Slots}{
	 \describe{
    \item{\code{design}:}{contains the design matrix}
    \item{\code{target}:}{contains the \code{target} data}
    \item{\code{covariates}:}{contains the \code{covariates}}
    \item{\code{intIndex}:}{contains the \code{intIndex}}
  }
}
\section{Methods}{
  \describe{
    \item{getCovariates}{\code{signature(object = "designMatrix")}: 
	    access the \code{covariates} slot }
    \item{getDesign}{\code{signature(object = "designMatrix")}:  
	    access the \code{design} slot}
    \item{getIntIndex}{\code{signature(object = "designMatrix")}:  
	    access the \code{intIndex} slot}
    \item{getTarget}{\code{signature(object = "designMatrix")}:  
	    access the \code{target} slot}
    \item{initialize}{\code{signature(.Object = "designMatrix")}:  
        create a new designMatrix class}
    \item{show}{\code{signature(object = "designMatrix")}:  
	    print the designMatrix class}
	 }
}
\author{Xiwei Wu, Arthur Li }
\seealso{
	\code{\linkS4class{contrastMatrix}} 
}
\examples{
data(eSetExample)
## One-way Anova
(design1<- new("designMatrix", target=pData(eSetExample), covariates = "Treatment"))

## Randomized block design
(design2<- new("designMatrix", target=pData(eSetExample), 
   covariates = c("Treatment", "Group")))

## Interaction design
(design3<- new("designMatrix", target=pData(eSetExample), 
   covariates = c("Treatment", "Group"), intIndex=c(1,2)))
   
}
\keyword{classes}
