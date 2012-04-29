\name{getRefId}
\alias{getRefId}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Returns the used reference sample in a normalized \code{MELON} object
}

\description{
Method to return the used reference sample index in a normalized \code{MELON} object
}

\usage{
getRefId(object)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{object of class \code{MELON}}
}

\details{
## %%  ~~ If necessary, more details than the description above ~~
}

\value{
index of used reference sample. Either a user defined sample can be used as a reference sample, but the sample with the largest amount of unique intensities is by default used as reference sample.
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
