% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{v.ks.test}
\alias{v.ks.test}
\title{Vectorized Kolmogorov-Smirnov Test}
\usage{
v.ks.test(data_vecs, ref_vecs, use_weights = TRUE, weights.pwr = 1)
}
\arguments{
\item{data_vecs}{A list of data vectors, each representing a gene.}

\item{ref_vecs}{A list of reference vectors, each representing a gene.}

\item{use_weights}{Logical, whether to use weights in the KS test. Default is TRUE.}

\item{weights.pwr}{Numeric, the power to which weights are raised. Default is 1.}
}
\value{
A list containing two elements:
\item{D_stats}{A numeric vector of D statistics for each gene}
\item{p_vals}{A numeric vector of p-values for each gene}
}
\description{
This function performs a vectorized Kolmogorov-Smirnov test on multiple gene vectors.
}
\details{
This function applies the Kolmogorov-Smirnov test to multiple gene vectors simultaneously.
It compares each data vector to its corresponding reference vector and calculates
the D statistic and p-value for each gene.

If \code{use_weights} is TRUE, the function applies weights to the KS test, with weights
ranging from 1 to -1 across the ranks.
}
