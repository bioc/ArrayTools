\name{regress}
\alias{regress}
\title{Run regression to fit genewise linear model}
\description{
  Fit genewise linear model using LIMMA package, ordinary linear regression, or permutation method. 
}
\usage{
regress(object, contrast, method = c("limma", "regression", "permutation"), adj = "none", permute.time = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an ExpressionSet}
  \item{contrast}{a contrastMatrix }
  \item{method}{choose the follwoing three options: "limma" (LIMMA), 
      "regression" (ordinary linear regression), "permutation" (permutation test)}
  \item{adj}{adjustment method for multiple comparison test, including "holm", 
      "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none". The default 
	  value is "none". Type help(p.adjust) for more detail.}
  \item{permute.time}{number of permutation times, only used for the "permutation" method}
}
\value{
  an object of \code{regressResult}
}
\author{Xiwei Wu, Arthur Li}
\examples{
data(eSetExample)
design<- new("designMatrix", target=pData(eSetExample), covariates = "Treatment")
contrast<- new("contrastMatrix", design.matrix = design, 
    compare1 = "Treated", compare2 = "Control")
result<- regress(eSetExample, contrast)
}
