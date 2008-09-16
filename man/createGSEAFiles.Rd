\name{createGSEAFiles}
\alias{createGSEAFiles}
\title{A Wrapper Function to create *.GCT and *.CLS for GSEA analysis}
\description{
  A Wrapper Function to create *.GCT and *.CLS for GSEA analysis
}
\usage{
createGSEAFiles(mydir = getwd(), eSet, catVar)
}
\arguments{
  \item{mydir}{directory where you would like to store the files}
  \item{eSet}{an ExpressionSet}
  \item{catVar}{variable of interest}
}
\value{
  Creating *.GCT and *.CLS for GSEA
}
\author{Xiwei Wu, Arthur Li}
\references{\code{http://www.broad.mit.edu/gsea/}}
\seealso{\code{\link{output.cls}}, \code{\link{output.gct}}}
\examples{
data(eSetExample)
\dontrun{createGSEAFiles (mydir, eSetExample, "Treatment")}
}

