% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/alm_cons.R
\name{cons_fun}
\alias{cons_fun}
\title{bi-factor constraint value function}
\usage{
cons_fun(L, J, G, Pair)
}
\arguments{
\item{L}{the input loading matrix}

\item{J}{the number of items}

\item{G}{the number of group factors}

\item{Pair}{the list indicating which pair of columns have constraints}
}
\value{
a matrix consisting of the values of equality constraints
}
\description{
bi-factor constraint value function
}
