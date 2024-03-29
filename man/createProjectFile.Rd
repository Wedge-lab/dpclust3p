% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocessing.R
\name{createProjectFile}
\alias{createProjectFile}
\title{Creates a DPClust project master file}
\usage{
createProjectFile(
  outputfile,
  donornames,
  samplenames,
  sex,
  purities = NULL,
  rho_and_psi_files = NULL,
  datafiles = paste(samplenames, "_allDirichletProcessInfo.txt", sep = ""),
  cndatafiles = NULL,
  indeldatafiles = NULL
)
}
\arguments{
\item{outputfile}{Full path with filename where the output will be written}

\item{donornames}{A vector with donor identifiers, use the same donor identifier to match multiple samplenames for a multi-sample DPClust run}

\item{samplenames}{A vector with sample identifiers}

\item{sex}{A vector with the sex of each donor}

\item{purities}{A vector with a purity value per sample}

\item{datafiles}{Vector with filenames in which the DPClust input is contained (Default: [samplename]_allDirichletProcessInfo.txt)}

\item{cndatafiles}{A vector with CNA DPClust input files (Default: NULL)}

\item{indeldatafiles}{A vector with indel DPClust input files (Default: NULL)}
}
\description{
Creates a DPClust project master file
}
\author{
sd11
}
