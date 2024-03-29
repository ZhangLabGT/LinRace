% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LikelihoodCal_Dif.R
\name{LikelihoodCal_Dif}
\alias{LikelihoodCal_Dif}
\title{Calculate expression likelihood with diffusion loss}
\usage{
LikelihoodCal_Dif(
  tree,
  muts,
  prob_t,
  cell_state_labels,
  state_lineages,
  newick_lookup,
  lambda1 = 10,
  lambda2 = 1
)
}
\arguments{
\item{tree}{input tree}

\item{muts}{lineage barcode matrix}

\item{prob_t}{transition probabilities between cells}

\item{cell_state_labels}{cell states of single cells}

\item{state_lineages}{the lineages that makes the state network from the root state to leaf states}

\item{newick_lookup}{lookup table for internal nodes}

\item{lambda1}{hyperparameter for asymmetric division likelihood}

\item{lambda2}{hyperparameter for neighbor distance likelihood}
}
\description{
This function computes the expression likelihood of a tree based on cell gene expression data
}
