\name{normScatterPlot}
\alias{normScatterPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create MELON scatter plots
}

\description{
Method to construct MELON scatter plots of normalized samples vs the reference sample 
}

\usage{
normScatterPlot(object, nres = 5, datatype = c("Count", "Continuous"), savedir = NULL, plotabline = T, tracefile = NULL)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{object of class \code{MELON}}
  \item{nres}{percentage, default value is 5\%}
  \item{datatype}{scale of normalized data, either \code{Count} (log) or \code{Continuous}}
  \item{savedir}{the directory where to save the plots}
  \item{plotabline}{logical indicating whether bisector should be shown. Default is TRUE.}
  \item{tracefile}{a writable connection or a character string naming the file to write to, or \code{NULL} to stop sink-ing.}
}

\details{
A technical description of the algorithm and performance is given in the article 'Methylation Enriched LOci Normalization (MELON) for genome-wide DNA methylation experiments' (De Meyer et al. 2012). 
}

\value{
The method provides \code{MELON} scatter plots of normalized samples vs the reference sample. The sample names are embedded in the file names.
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
