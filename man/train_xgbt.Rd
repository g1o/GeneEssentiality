% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/train_xgbt.R
\name{train_xgbt}
\alias{train_xgbt}
\title{Generate 1 XGBT model}
\usage{
train_xgbt(
  features = data,
  seeds = seed,
  CPU = 4,
  CV = 10,
  nrepeats = 3,
  saveprediction = "final"
)
}
\arguments{
\item{CPU}{Number of threads to use when training the model}

\item{CV}{Number of sets for the cross-validation}

\item{nrepeats}{}

\item{data}{A feature data frame including the label of two Classes}
}
\value{
A list of models
}
\description{
Generate 1 XGBT model
}
