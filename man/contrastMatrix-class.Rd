\name{contrastMatrix}
\docType{class}
\alias{contrastMatrix-class}
\alias{contrastMatrix}
\alias{class:contrastMatrix}
\alias{getCompare1,contrastMatrix-method}
\alias{getCompare2,contrastMatrix-method}
\alias{getContrast,contrastMatrix-method}
\alias{getInteraction,contrastMatrix-method}
\alias{getLevel,contrastMatrix-method}
\alias{initialize,contrastMatrix-method}
\alias{show,contrastMatrix-method}

\title{Class to Contain the Contrast Matrix that Used for Linear Regression}
\description{Class to Contain the Contrast Matrix that Used for Linear Regression,
    inherited from the designMatrix class}


\section{Creating Objects}{

    \code{new("contrastMatrix", ..., 
	    design.matrix=[designMatrix], 
		compare1=[character], 
		compare2=[character], 
		level=[character],
		interaction==[numeric])}.
    
    This creates a contrast matrix class. \code{design.matrix} is a \code{designMatrix}
	class.  \code{compare1} the first value of the main covariate, and \code{compare2}
	is the second value of the main covariate.  
	For example, suppose that the main covariate is "drug", and there 
	are three unique values: "drug1", "drug2", and "placebo". 
	You would like to compare "drug1" to "drug2". Then you would use "drug1" as 
	\code{compare1} and "drug2" as \code{compare2}.  If \code{interaction==TRUE}, 
	do not specify \code{compare1} and \code{compare2}.  You only specify 
	\code{level} when the design matrix contains an interaction term. 
	Suppose that you would like to compare "drug1" to "drug2" only when estrogen is 
	"present", where "present" is one of the values of the estrogen variable. 
	You will use "present" as level. If \code{interaction==TRUE}, do not specify this 
	value as well.  You only specify \code{interaction=TRUE} when you would like to 
	detect the interaction effect between two covariates. In this case,
	do not provide values for \code{compare1}, \code{compare2}, and \code{level}
}

\section{Slots}{
	 \describe{
    \item{\code{contrast}:}{Object of class \code{"matrix"} contains the contrast matrix}
    \item{\code{compare1}:}{Object of class \code{"character"} contains \code{compare1}}
    \item{\code{compare2}:}{Object of class \code{"character"} contains \code{compare2} }
    \item{\code{level}:}{Object of class \code{"character"} contains \code{level} }
    \item{\code{interaction}:}{Object of class \code{"logical"} contains \code{interaction}}
    \item{\code{design}:}{Object of class \code{"matrix"} contain the design matrix}
    \item{\code{target}:}{Object of class \code{"data.frame"} contains \code{target}}
    \item{\code{covariates}:}{Object of class \code{"character"} contains \code{covariates}}
    \item{\code{intIndex}:}{Object of class \code{"numeric"} contains \code{intIndex} }
  }
}
\section{Extends}{
Class \code{"\linkS4class{designMatrix}"}, directly.
}
\section{Methods}{
  \describe{
    \item{getCompare1}{\code{signature(object = "contrastMatrix")}: 
	    access the \code{compare1} slot}
    \item{getCompare2}{\code{signature(object = "contrastMatrix")}: 
	    access the \code{compare2} slot}
    \item{getContrast}{\code{signature(object = "contrastMatrix")}: 
	    access the \code{contrast} slot }
    \item{getInteraction}{\code{signature(object = "contrastMatrix")}: 
		access the \code{interaction} slot}
    \item{getLevel}{\code{signature(object = "contrastMatrix")}: 
	    access the \code{level} slot }
    \item{initialize}{\code{signature(.Object = "contrastMatrix")}: 
	    create a new contrast matrix class }
    \item{show}{\code{signature(object = "contrastMatrix")}: print the contrast matrix }
	 }
}
\author{ Xiwei Wu, Arthur Li }

\seealso{
	\code{\linkS4class{designMatrix}} 
}
\examples{
data(eSetExample)
## One-way Anova
(design1<- new("designMatrix", target=pData(eSetExample), covariates = "Treatment"))
(contrast1<- new("contrastMatrix", design.matrix = design1, 
    compare1 = "Treated", compare2 = "Control"))
	
## Randomized block design
(design2<- new("designMatrix", target=pData(eSetExample), 
   covariates = c("Treatment", "Group")))
(contrast2<- new("contrastMatrix", design.matrix = design2, 
    compare1 = "Treated", compare2 = "Control"))

## Interaction design
(design3<- new("designMatrix", target=pData(eSetExample), 
   covariates = c("Treatment", "Group"), intIndex=c(1,2)))
# Test for interaction:
(contrast.int<- new("contrastMatrix", design.matrix = design3,
    interaction=TRUE))
# Compare Treated vs Control among group A
(contrast.a<- new("contrastMatrix", design.matrix = design3, 
    compare1 = "Treated", compare2 = "Control", level="A"))
}
\keyword{classes}
