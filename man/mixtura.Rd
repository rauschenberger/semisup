% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{mixtura}
\alias{mixtura}
\title{Model fitting}
\usage{
mixtura(y, z, dist = "norm",
        phi = NULL, pi = NULL, gamma = NULL,
        test = NULL, iter = 100, kind = 0.05,
        debug = TRUE, ...)
}
\arguments{
\item{y}{\strong{observations:}
numeric vector of length \code{n}}

\item{z}{\strong{class labels:}
integer vector of length \code{n},
with entries \code{0}, \code{1} and \code{NA}}

\item{dist}{distributional assumption\strong{:}
character \code{"norm"} (Gaussian),
\code{"nbinom"} (negative bionomial),
or \code{"zinb"} (zero-inflated negative binomial)}

\item{phi}{dispersion parameters\strong{:}
numeric vector of length \code{q},
or \code{NULL}}

\item{pi}{zero-inflation parameter(s)\strong{:}
numeric vector of length \code{q},
or \code{NULL}}

\item{gamma}{offset\strong{:}
numeric vector of length \code{n},
or \code{NULL}}

\item{test}{resampling procedure\strong{:}
character \code{"perm"} (permutation) or
\code{"boot"} (parametric bootstrap),
or \code{NULL}}

\item{iter}{(maximum) number of resampling iterations \strong{:}
positive integer, or \code{NULL}}

\item{kind}{resampling accuracy\strong{:}
numeric between \code{0} and \code{1}, or \code{NULL}\strong{;}
all \code{p}-values above \code{kind} are approximate}

\item{debug}{verification of arguments\strong{:}
\code{TRUE} or \code{FALSE}}

\item{...}{settings \code{EM} algorithm\strong{:}
\code{starts}, \code{it.em} and \code{epsilon}
(see \code{\link{arguments}})}
}
\value{
This function fits and compares a one-component (\code{H0})
and a two-component (\code{H1}) mixture model.

\item{posterior}{probability of belonging to class 1\strong{:}
numeric vector of length \code{n}}
\item{converge}{path of the log-likelihood\strong{:}
numeric vector with maximum length
\code{it.em}}
\item{estim0}{parameter estimates under \code{H0}\strong{:}
data frame}
\item{estim1}{parameter estimates under \code{H1}\strong{:}
data frame}
\item{loglik0}{log-likelihood under \code{H0}\strong{:}
numeric}
\item{loglik1}{log-likelihood under \code{H1}\strong{:}
numeric}
\item{lrts}{likelihood-ratio test statistic\strong{:}
positive numeric}
\item{p.value}{\code{H0} versus \code{H1}\strong{:}
numeric between \code{0} and \code{1}, or \code{NULL}}
}
\description{
This function fits a semi-supervised mixture model.
It simultaneously estimates two mixture components,
and assigns the unlabelled observations to these.
}
\details{
By default, \code{phi} and \code{pi}
are estimated by the maximum likelihood method,
and \code{gamma} is replaced by a vector of ones.
}
\section{Reference}{

A Rauschenberger, RX Menezes, MA van de Wiel,
NM van Schoor, and MA Jonker (2020).
"Semi-supervised mixture test for
detecting markers associated with a quantitative trait",
\emph{Manuscript in preparation}.
}

\examples{
# data simulation
n <- 100
z <- rep(0:1,each=n/2)
y <- rnorm(n=n,mean=2,sd=1)
z[(n/4):n] <- NA

# model fitting
mixtura(y,z,dist="norm",test="perm")

}
\seealso{
Use \code{\link{scrutor}} for hypothesis testing.
All other functions are \code{\link{internal}}.
}
\keyword{methods}
