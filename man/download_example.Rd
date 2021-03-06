% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genie_functions.R
\name{download_example}
\alias{download_example}
\title{Downloads example data for rgenie.}
\usage{
download_example(dir, name = "MUL1", overwrite = FALSE, quiet = FALSE)
}
\arguments{
\item{dir}{Directory where example data should be put.}

\item{name}{The name of the example to download.}

\item{overwrite}{If FALSE, then data are not downloaded if directory 'dir/name' already
exists. Otherwise, all data are downloaded (possibly overwriting files).}

\item{quiet}{If TRUE, then no messages are printing during the analysis.}
}
\value{
Returns a list containing the same tables as in an individual result,
but concatenated across regions.
}
\description{
The example data is a set of BAM files for GenIE replicates.
}
\examples{
\donttest{
download_example(dir = "~/genie_example", name = "MUL1")
# Data are downloaded and we can run an rgenie analysis
setwd("~/genie_example/MUL1/")
regions = readr::read_tsv("mul1.genie_regions.tsv")
replicates = readr::read_tsv("mul1.genie_replicates.tsv")
grep_results = grep_analysis(regions, replicates)
}
}
\seealso{
\code{\link{grep_analysis}}

\code{\link{deletion_analysis}}
}
