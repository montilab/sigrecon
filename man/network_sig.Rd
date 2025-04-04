% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/baseline_methods.R
\name{network_sig}
\alias{network_sig}
\title{Finds a simulated network signature}
\usage{
network_sig(net, seeds, sig = c("corr", "rwr"), p = 0.1, limit = 30)
}
\arguments{
\item{net}{network given as an igraph}

\item{seeds}{single gene string or vector of gene strings or named list of sets of gene symbols}

\item{sig}{String, network signature type: random walk, correlation, nearest neighbor}

\item{p}{Numeric, restart value for random walk, default=0.1}

\item{limit}{number of genes to be included in the network signature, default=30}
}
\value{
vector of gene strings
}
\description{
Finds a simulated network signature
}
