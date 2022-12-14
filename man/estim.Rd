% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Estimation.R
\name{estim}
\alias{estim}
\title{Structural change model estimation}
\usage{
estim(m, q, z, y, b, robust, prewhit, hetomega, hetq, x, p, hetdat, hetvar)
}
\arguments{
\item{m}{number of break}

\item{q}{number of regressors z}

\item{z}{matrix of regressor z with coefficients are allowed to change across
regimes}

\item{b}{break dates}

\item{prewhit}{option to use prewhitening process based on AR(1) approximation}

\item{hetomega, hetq, hetdat, hetvar}{options for assumptions on the error terms}

\item{x}{matrix of regressor x with coefficients are constant across regimes}

\item{p}{number of regressors x}
}
\value{
A list containing the following components:
\itemize{
\item{date} {List of estimated breaks}
\item{CI} {List of Confidence Intervals for each corresponding break}
\item{beta} {Estimated coefficients of the regression. The first
(\code{m}+1)*\code{q} are coefficients of \code{q} variables \code{z} that change across regimes.
The last \code{p} are coefficients of \code{p} variables \code{x}
that are constant across regimes}
\item{SE} {Corrected standard errors of the corresponding coefficients}}
}
\description{
Function to estimate the model by OLS given the obtained break dates
It also computes and reports confidence intervals for the break dates and
corrected standard errors of coefficients estimates given specifications of
errors covariance matrix via robust, hetomega, hetq, hetdat and hetvar
}
