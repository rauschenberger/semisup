
#-------------------------------------------------------------------------------
#--- PACKAGE NAME --------------------------------------------------------------
#-------------------------------------------------------------------------------

#' @title
#' Semi-supervised mixture model
#' 
#' @description
#' This R package implements the semi-supervised mixture model.
#' Use \code{\link{mixtura}} for model fitting,
#' and \code{\link{scrutor}} for hypothesis testing.
#' 
#' @name semisup-package
#' @aliases semisup
#' @keywords documentation
#' @docType package
#'  
#' @section Getting started:
#' Please type the following commands\strong{:} \cr
#' \code{utils::vignette("semisup")} \cr
#' \code{?semisup::mixtura} \cr
#' \code{?semisup::scrutor}
#' 
#' @section More information:
#' A Rauschenberger, RX Menezes, MA van de Wiel,
#' NM van Schoor, and MA Jonker (2017).
#' "Detecting SNPs with interactive effects on a quantitative trait",
#' \emph{Manuscript in preparation}.
#' 
#' \email{a.rauschenberger@vumc.nl}
NULL 

#-------------------------------------------------------------------------------
#--- ARGUMENTS -----------------------------------------------------------------
#-------------------------------------------------------------------------------

#' @title
#' Documentation
#'
#' @description
#' This page lists and describes all arguments
#' of the R package \code{\link{semisup}}.
#' 
#' @name arguments
#' @keywords internal
#' 
#' @param y
#' \strong{observations:}
#' numeric vector of length \code{n}
#' 
#' @param Y
#' \strong{observations:}
#' numeric vector of length \code{n},
#' or numeric matrix with \code{n} rows (samples)
#' and \code{q} columns (variables)
#' 
#' @param z
#' \strong{class labels:}
#' integer vector of length \code{n},
#' with entries \code{0}, \code{1} and \code{NA}
#' 
#' @param Z
#' \strong{class labels:}
#' numeric vector of length \code{n},
#' or numeric matrix with \code{n} rows (samples)
#' and \code{p} columns (variables),
#' with entries \code{0} and \code{NA}
#'  
#' @param dist
#' distributional assumption\strong{:}
#' character \code{"norm"} (Gaussian),
#' \code{"nbinom"} (negative bionomial),
#' or \code{"zinb"} (zero-inflated negative binomial)
#' 
#' @param phi
#' dispersion parameter\strong{:}
#' positive numeric,
#' or \code{NULL}
#' 
#' @param phi
#' dispersion parameters\strong{:}
#' numeric vector of length \code{q},
#' or \code{NULL}
#' 
#' @param pi
#' zero-inflation parameter\strong{:}
#' numeric between 0 and 1,
#' or \code{NULL}
#' 
#' @param pi
#' zero-inflation parameter(s)\strong{:}
#' numeric vector of length \code{q},
#' or \code{NULL}
#'
#' @param gamma
#' offset\strong{:}
#' numeric vector of length \code{n},
#' or \code{NULL}
#' 
#' @param test
#' resampling procedure\strong{:}
#' character \code{"perm"} (permutation) or
#' \code{"boot"} (parametric bootstrap),
#' or \code{NULL}
#' 
#' @param iter
#' (maximum) number of resampling iterations \strong{:}
#' positive integer, or \code{NULL}
#' 
#' @param kind
#' resampling accuracy\strong{:}
#' numeric between \code{0} and \code{1}, or \code{NULL}\strong{;}
#' all \code{p}-values above \code{kind} are approximate
#' 
#' @param starts
#' restarts of the \code{EM} algorithm\strong{:}
#' positive integer (defaults to \code{1})
#' 
#' @param it.em
#' (maximum) number of iterations in the \code{EM} algorithm\strong{:}
#' positive integer (defaults to \code{100})
#' 
#' @param epsilon
#' convergence criterion for the \code{EM} algorithm\strong{:}
#' non-negative numeric (defaults to \code{1e-04})
#'
#' @param debug
#' verification of arguments\strong{:}
#' \code{TRUE} or \code{FALSE}
#' 
#' @param pass
#' parameters for parametric bootstrap algorithm
#' 
#' @param ...
#' settings \code{EM} algorithm\strong{:}
#' \code{starts}, \code{it.em} and \code{epsilon}
#' (see \code{\link{arguments}})
#' 
#' @seealso
#' Use \code{\link{mixtura}} for model fitting,
#' and \code{\link{scrutor}} for hypothesis testing.
#' All other functions of the R package \code{\link{semisup}}
#' are \code{\link{internal}}.
NULL

#-------------------------------------------------------------------------------
#--- INTERNAL FUNCTIONS --------------------------------------------------------
#-------------------------------------------------------------------------------

#' @title
#' Documentation
#' 
#' @name internal
#' @keywords internal
#'
#' @description
#' This page lists and describes some internal functions
#' of the R package \code{\link{semisup}}.
#' These functions should not be used for analysing data.
#' 
#' \code{\link{fit.wrap}} multiple restarts
#' \cr \code{\link{fit.norm}} Gaussian mixture model
#' \cr \code{\link{fit.nbinom}} negative binomial mixture model
#' \cr \code{\link{fit.zinb}} zero-inflated negative binomial mixture model
#' \cr \code{\link{estim.nbinom}} dispersion estimation
#' \cr \code{\link{estim.zinb}} dispersion and zero-inflation estimation
#' \cr \code{\link{resam.lrts}} resampling (bootstrap, permutation)
#' 
#' @seealso 
#' Use \code{\link{mixtura}} for model fitting,
#' and \code{\link{scrutor}} for hypothesis testing.
NULL

#-------------------------------------------------------------------------------
#--- DATASETS ------------------------------------------------------------------
#-------------------------------------------------------------------------------

#' @title
#' Table
#' 
#' @name table
#' @keywords internal
#' @docType data
#' 
#' @description
#' This dataset includes tables for the approximate mixture test
#' (\strong{not yet available}).
#'
#' @usage data(table)
#' @format A list of numeric vectors.
#' @return All entries are numeric.
NULL

#' @title
#' Toydata
#' 
#' @name toydata
#' @keywords internal
#' @docType data
#' 
#' @description
#' This dataset allows to reproduce the examples shown in the vignette.
#' 
#' @usage data(toydata)
#' @format A list of numeric vectors and matrices.
#' @return All entries are numeric.
NULL
