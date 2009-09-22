\name{qa3prime}
\alias{qa3prime}
\title{Creating Quality Assessment Report for 3 Prime Array}
\description{
  Creating Quality Assessment Report for 3 Prime Array in HTML file
}
\usage{
qa3prime(object, parameters, outputFile = "QA.html", mydir = getwd())
}
\arguments{
  \item{object}{an \code{AffyBatch} object}
  \item{parameters}{The names of the variables to be included in the report }
  \item{outputFile}{The name of the outputfile.  Make sure write ".html"}
  \item{mydir}{The name of the directory containing the report}
}
\details{

This function creates quality control report in an HTML file that contains a 
set of 9 assessment figures.

Figure1: The Raw Intensity Plot. The raw intensity should be similar across all chips

Figure2: The Average Background/Percentage Present Plot.  The Average Background 
should be similar across all chips. The Percentage Present should be similar 
across all chips, except that in rare situations transcription is globally shut 
down or turned on under some conditions	
	
Figure3: The Scaling Factor Plot.  The scaling factor should be within 
3-fold across all chips

Figure4: The Hybridization Controls Plot.  BioB, BioC, BioD, CreX should be called 
present, except that it is acceptable if BioB is absent sometimes.

Figure5: The Housekeeping Controls Plot.  The GAPDH ratio should be around 1 
and the actin ratio should be less than 3. Note that if two-cycle amplification 
or NuGen amplification is used, this ratio could be much higher.

Figure6: The RNA Degradation Plot.  On Affymetrix GeneChips, individual probes 
in a probeset are ordered by location relative to the 5` end of the targeted 
RNA molecule. On each chip, probe intensities are averaged by location in the 
probeset, with the average taken over probesets. In an RNA digestion plot, 
these means are plotted side-by-side, making it easy to notice any 5` to 3` trend. 
The trend can be due to RNA degradation or 3`-biased amplification. Since RNA 
degradation typically starts from the 5` end of the molecule and amplification 
starts at the 3` end, we would expect probe intensities to be systematically 
lowered at the 5` end of a probeset when compared to the 3` end. 

Figure7: The Hierarchical Clustering of Samples.  Samples will be grouped using 
hierarchical clustering and principal component analysis (PCA). If the sample 
preparation steps introduced bigger variation than biological variation, 
treatment groups will be mixed up in the plot. This could also happen when the 
samples between groups were mixed up accidentally when the samples were prepared. 
We acknowledge that clinical samples are harder to collect and sometimes impossible 
to control. Therefore, sample QC criteria will be much looser when dealing with 
clinical samples. 


Figure8: The Pseudo-chip Images.  A Pseudo-chip image plots the weights and 
residuals from the model fit. The image plot allows detection of artifacts on the chip.

Figure9: The Normalized Unscaled Standard Error (NUSE) and Relative Log Expression 
(RLE) Plots.  The NUSE is fitted robustly by iteratively reweighted least squares 
(IRLS) so that the standard error of the estimated log2 scale expression can be estimated. 
The boxplots of the NUSE show the differences in hybridization quality most clearly, 
in magnitude as well as variability. A high NUSE likely corresponds to a low signal. 
The RLE plot is a boxplot showing the distribution of Log2 ratio of each chip relative 
to a median chip. A discordant distribution infers a problem with the chip.
}
 
\value{
  no value is returned
}
\references{ \url{http://www.affymetrix.com} }
\author{Xiwei Wu, Arthur Li}
\examples{
  \dontrun{qa3prime(AffyBatchExample, c("var1", "var2"))}
}
