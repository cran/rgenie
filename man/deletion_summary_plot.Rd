% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{deletion_summary_plot}
\alias{deletion_summary_plot}
\title{Plots a summary of deletion analysis results for a single region.}
\usage{
deletion_summary_plot(del_result)
}
\arguments{
\item{del_result}{Result from a call to deletion_analysis.}
}
\value{
Returns a ggplot object with a summary of deletion analysis results for a single region.
}
\description{
Plots a summary of deletion analysis results for a single region.
}
\examples{
# Note: First run deletion_analysis()
# mul1_del_results is a pre-loaded result

deletion_summary_plot(mul1_del_results[[1]])
}
\seealso{
\code{\link{deletion_analysis}}
}
