% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/femeta.R
\name{femeta}
\alias{femeta}
\title{Generic Fixed-effects Meta-analysis with Several 2x2 Tables}
\usage{
femeta(x, ...)
}
\arguments{
\item{x}{an object used to select a method. Currently x can be either a 2x2xK array, which links to \link{femeta.array} method,
or a integer vector, which leads to the use of \link{femeta.default}.}

\item{...}{further arguments passed to or from other methods.}
}
\description{
Generic function for obtaining the confidence interval of logarithm of Cochran-Mantel-Haenszel estimator,
Woolf's estimator of MLE of an assumed-common log odds ratio without assuming
homogeneity
}