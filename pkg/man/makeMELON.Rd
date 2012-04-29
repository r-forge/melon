\name{makeMELON}
\alias{makeMELON}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Normalize data from DNA methylation experiments according to MELON
}

\description{
Function to create a normalized \code{MELON} object for DNA methylation experimental data
}

\usage{
makeMELON(object, method = c("Max","Clean"), regrwindow = 201, minC = 50, manC = NULL, maxL = 20000, refsample = NULL, doplot = F, savedir = NULL, savefile = NULL, tracefile = NULL)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{raw data, either a dataframe or a matrix}
  \item{method}{indicating the selected criterion to determine the optimal number of used heavily methylated loci. Current implementations foresee \code{Max}, corresponding with the highest CME, and \code{Clean}, corresponding with the highest increase in consistent methylation enrichment (CME) for increasing L.}
  \item{regrwindow}{local regression window, parameter needed for method \code{Clean}. The parameter should be a positive odd integer, default value is 201.}
  \item{minC}{minimum number of stable loci needed in the algorithm. The default is 50.}
  \item{manC}{number of stable loci desired for the algorithm. The default is the number of loci in the full dataset.}
  \item{maxL}{maximal number of considered heavily methylated loci in the algorithm. The default value is 20000.}
  \item{refsample}{index of sample identified as reference sample. If missing, the sample with the largest amount of unique intensities is used as reference.}
  \item{doplot}{logical indicating whether a stability plot should be drawn. This figure plots the consistent methylation enrichment (CME) quantity against the most heavily methylated loci (L) and displays the finally selected L and corresponding CME.}
  \item{savedir}{the directory where to save the plot if \code{doplot == T}}
  \item{savefile}{the filename for the plot if \code{doplot == T}. If \code{doplot == T}, the default name is StabilityPlot (pdf).}
  \item{tracefile}{a writable connection or a character string naming the file to write to, or \code{NULL} to stop sink-ing.}
}

\details{
A technical description of the algorithm and performance is given in the article 'Methylation Enriched LOci Normalization (MELON) for genome-wide DNA methylation experiments' (De Meyer et al. 2012).
}

\value{
An object of class \code{MELON}
}

\references{
De Meyer Tim, Pipelers Peter, Szarc vel Szic Katarzyna, Denil Simon, Trooskens Geert, Renard Jean-Pierre, De Keulenaer Sarah, Vanden Berghe Wim and Van Criekinge Wim (2012) Methylation Enriched LOci Normalization (MELON) for genome-wide DNA methylation experiments.
}

\author{
Peter Pipelers <Peter.Pipelers@ugent.be>%%  ~~who you are~~
}

\note{
## %%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
## %% ~~objects to See Also as \code{\link{help}}, ~~~
}

\examples{
##
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{
##
}
