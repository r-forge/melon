\name{MELON-class}
\Rdversion{1.1}
\docType{class}
\alias{MELON-class}
\alias{show,MELON-method}
\alias{initialize,MELON-method}
\alias{getNormData,MELON-method}
\alias{getNormFactors,MELON-method}
\alias{getIdStabLoci,MELON-method}
\alias{getDataStabLoci,MELON-method}
\alias{getRefId,MELON-method}
\alias{getStabFeatures,MELON-method}
\alias{normScatterPlot,MELON-method}

\title{
Class "MELON"
}

\description{
Class to represent normalized data from DNA methylation experiments according to the MELON algorithm.
}

\section{Slots}{
  \describe{
    \item{\code{normData}:}{\code{"matrix"} ~~ }
    \item{\code{normFactors}:}{numeric vector ~~ }
    \item{\code{idStabLoci}:}{numeric vector ~~ }
    \item{\code{dataStabLoci}:}{\code{"matrix"} ~~ }
    \item{\code{stabFeatures}:}{\code{"matrix"} ~~ }
    \item{\code{refId}:}{numeric ~~ }
  }
}

\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "MELON")}: ... }
    \item{initialize}{\code{signature(object = "MELON")}: ... }
    \item{getNormData}{\code{signature(object = "MELON")}: ... }
    \item{getNormFactors}{\code{signature(object = "MELON")}: ... }
    \item{getIdStabLoci}{\code{signature(object = "MELON")}: ... }
    \item{getDataStabLoci}{\code{signature(object = "MELON")}: ... }
    \item{getRefId}{\code{signature(object = "MELON")}: ... }
    \item{getStabFeatures}{\code{signature(object = "MELON")}: ... }
    \item{normScatterPlot}{\code{signature(object = "MELON")}: ... }
  }
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
