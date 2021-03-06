% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hrm.R
\name{alleles}
\alias{alleles}
\title{\code{alleles} is used to get the sequences fragments of the alleles to be analysed}
\usage{
alleles(loci_info, acceptedChar = c("A", "a", "C", "c", "T", "t", "G",
  "g", "-"))
}
\arguments{
\item{loci_info}{is the csv file with the following header filename,allele_length,fragment_name,fragment_start,fragment_stop
and each rows contain relevant information with ',' seperating the field}

\item{acceptedChar}{is the characters which are considered valid}
}
\value{
Will returns dataframe of different fragment with matching length, their GC content (curve number) and their fragment sequences
}
\description{
\code{alleles} is used to get the sequences fragments of the alleles to be analysed
}
