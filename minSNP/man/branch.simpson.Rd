% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dMode.R
\name{branch.simpson}
\alias{branch.simpson}
\title{\code{branch.simpson} is used to calculate 1 or more simpson index and list the position with highest index.}
\usage{
branch.simpson(seq, level = 1, included = NULL, excluded = NULL,
  numRes = 1)
}
\arguments{
\item{seq}{is fastaDNA object to analysed.}

\item{level}{is the number of positions.}

\item{included}{is the included position for analysed, these will force the computation
to compute the simpson's index at the position no matter what}

\item{excluded}{is the positions that are excluded from computation}

\item{numRes}{is the number of result to be returned.}
}
\value{
Will returns simpson's index and the position.
}
\description{
\code{branch.simpson} is used to calculate 1 or more simpson index and list the position with highest index.
}
