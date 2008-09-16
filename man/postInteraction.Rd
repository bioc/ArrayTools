\name{postInteraction}
\alias{postInteraction}
\title{Create an Object of InteractionResult Class for Testing Interaction}
\description{
   Based on the result from the interaction test by looking at the result
   from the regressResult object, this function partitions tne orignal data, an 
   ExpressionSetinto groups, one contains the genes without the interaction and others
   contains the genes with the interaction across different level of covariates.
}
\usage{
postInteraction(eSet, regressObject, mainVar, compare1, compare2, method = regressionMethod(regressObject), adj = adjustment(regressObject))
}
\arguments{
  \item{eSet}{an ExpressionSet}
  \item{regressObject}{a regressResult }
  \item{mainVar}{variable of main interest}
  \item{compare1}{the first value of the \code{mainVar}. For example, suppose that 
      \code{mainVar} is "drug", and there are three unique values: "drug1", 
	  "drug2", and "placebo". You would like to compare "drug1" to "drug2". 
	  Then you would use "drug1" as \code{compare1}}
  \item{compare2}{Based on the example for \code{compare1}, "drug2" will be
      the \code{compare2}}
  \item{method}{It is used to run regression within each level of the effect modifier.
      choose the follwoing three options: "limma" (LIMMA), 
      "regression" (ordinary linear regression), "permutation" (permutation test)}
  \item{adj}{adjustment method for multiple comparison test, including "holm", 
      "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none". The default 
	  value is "none". Type help(p.adjust) for more detail.}
}
\value{
  an \code{interactionResult} class.  The first component contains all the result
  for all the genes.  The second component contains the genes without the interaction
  effect.  The rest of the components contains genes with the interactions.  
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
}
