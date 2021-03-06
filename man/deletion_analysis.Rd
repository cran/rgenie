% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genie_functions.R
\name{deletion_analysis}
\alias{deletion_analysis}
\title{Alignment-based GenIE analysis}
\usage{
deletion_analysis(
  regions,
  replicates,
  required_match_left = 10,
  required_match_right = 10,
  crispr_del_window = 100,
  min_mapq = 0,
  max_mismatch_frac = 0.05,
  min_aligned_bases = 50,
  exclude_multiple_deletions = FALSE,
  exclude_nonspanning_reads = TRUE,
  allele_profile = FALSE,
  del_span_start = -20,
  del_span_end = 20,
  quiet = FALSE
)
}
\arguments{
\item{regions}{A data.frame defining GenIE regions.}

\item{replicates}{A data.frame defining GenIE replicates.}

\item{required_match_left}{The length of sequence to the left of the HDR site that must exactly match to identify HDR or WT reads.}

\item{required_match_right}{The length of sequence to the right of the HDR site that must exactly match to identify HDR or WT reads.}

\item{crispr_del_window}{The window around the cut site within which any deletion is considered a CRISPR deletion.
deletions that do not span the region [cut_site - crispr_del_window, cut_site + crispr_del_window]
will be ignored; i.e. such reads can be considered HDR or WT.}

\item{min_mapq}{The minimum mapping quality for reads to be included in the analysis.}

\item{max_mismatch_frac}{The maximum fraction of mismatches a read can have and be included in the analysis.}

\item{min_aligned_bases}{The minimum number of aligned bases within the region of interest for a read to be included in the analysis.}

\item{exclude_multiple_deletions}{If TRUE, then reads with multiple deletions will be excluded from the analysis.}

\item{exclude_nonspanning_reads}{If TRUE, then reads are excluded if their alignment does not overlap the region's
highlight_site (or cut site if no highlight_site is specified)}

\item{allele_profile}{If TRUE, then the result object will contain data.frames named site_profiles and mismatch_profiles, as detailed in the description below.}

\item{del_span_start}{An integer that specifies the start of a window, relative to the region's highlight site, within which deletions are counted.}

\item{del_span_end}{An integer that specifies the end of a window, relative to the region's highlight site, within which deletions are counted.}

\item{quiet}{If TRUE, then no messages are printing during the analysis.}
}
\description{
For each replicate associated with an input region, \code{deletion_analysis}
identifies HDR or WT sequences, as well as deletion alleles, and returns
statistics that indicate for every allele whether the allele:WT ratio
differs in cDNA and gDNA. Statistics are also computed for deletion alleles
aggregated together.
}
\section{Details}{
 For a deletion analysis, the regions parameter is a data.frame
  with a format as follows. All of the column names below must be specified.
  \tabular{lccccclll}{ name \tab sequence_name \tab start \tab end \tab
  highlight_site \tab cut_site \tab hdr_allele_profile \tab wt_allele_profile
  \tab ref_sequence \cr MUL1_rs6700034  \tab MUL1  \tab 1 \tab 21 \tab 11
  \tab 9 \tab ----------A---------- \tab ----------C---------- \tab
  ACCGCACCCCCCCGGCCTAAC \cr }

\tabular{rl}{
  \strong{name} \tab A unique identifier for the region. \cr
  \strong{sequence_name} \tab The chromosome or amplicon sequence name. \cr
  \strong{start} \tab The start coordinate of the amplicon relative to the chromosome or amplicon reference. \cr
  \strong{end} \tab The end coordinate of the amplicon relative to the chromosome or amplicon reference (the end coordinate is included in the region). \cr
  \strong{highlight_site} \tab The relative position of the site of interest, usually the HDR SNP site. \cr
  \strong{cut_site} \tab The relative position of the cut site. \cr
  \strong{ref_sequence} \tab The reference sequence for the amplicon region, which must have length (end - start + 1). \cr
  \strong{hdr_allele_profile} \tab An allele profile describing the HDR allele. See details below. \cr
  \strong{wt_allele_profile} \tab An allele profile describing the WT allele. See details below. \cr
}

If multiple rows are defined for `regions`, then a separate analysis is run
for each region, using matched replicates from `replicates`.

The allele_profile columns indicate the positions in the amplicon sequence
that must match a given nucleotide for a read to be considered either HDR
or WT. This sequence must be the same length as the reference sequence,
and all other positions should be "-". The total sequence region that must
match is determined by both the positions of specified nucleotides and by
the \strong{required_match_left} and \strong{required_match_right}
parameters. These parameters give the length of sequence which must match
the provided reference sequence to the left of the leftmost specified
nucleotide, or to the right of the rightmost specified nucleotide.


The replicates parameter is a data.frame with a format as below.
\tabular{lccl}{ name \tab replicate \tab type \tab bam \cr MUL1_rs6700034
\tab c1.2 \tab cDNA \tab
bam_amplicon/MUL1_rs6700034_cDNA_rep1_pcr2.sortedByCoord.bam \cr
MUL1_rs6700034 \tab c1.3 \tab cDNA \tab
bam_amplicon/MUL1_rs6700034_cDNA_rep1_pcr3.sortedByCoord.bam \cr }

\tabular{rl}{
  \strong{name} \tab Indicates the region that a given replicate
  corresponds with. All replicates matching the name in the regions table
  will be used. \cr
  \strong{replicate} \tab an ID for the replicate, which must
  be unique among replicates for the region. \cr
  \strong{type} \tab Must have the value "cDNA" or "gDNA", indicating
  whether a given replicate contains data for cDNA or gDNA. \cr
  \strong{bam} \tab the path (relative to the working directory) to
  a BAM file with sequencing reads for the replicate. \cr
}

Statistics can only be computed if there are at least 2 replicates of each
type (cDNA and gDNA). Replicates are matched to the region based on the
\strong{name} column.
}

\section{Results}{
 The returned object is a list, where each item is the
result for one region. The result for a region (e.g. results[[1]]) is itself
a list, with the following items:
\tabular{rl}{
  \strong{region_stats} \tab Main analysis output, with statistics indicating whether the HDR/WT levels differ in cDNA relative to gDNA. \cr
  \strong{replicate_stats} \tab A data.frame with a row for each replicate, which has counts of reads in different categories and some summary values. \cr
  \strong{region} \tab 	Details of the input region the result corresponds to. \cr
  \strong{replicates} \tab Details of the input replicates the result corresponds to. \cr
  \strong{opts} \tab A list containing the options that were given for the analysis. \cr
  \strong{type} \tab Has the value “deletion_analysis”, and is used by plotting functions that take a full del_result list as input. \cr
}

The main output of interest is the `region_stats` field, which is a one-row data.frame with the following values:

\tabular{rl}{
  \strong{name} \tab Name of the region. \cr
  \strong{hdr_rate_gDNA} \tab Fraction of reads in gDNA identified as HDR. \cr
  \strong{hdr_rate_cDNA} \tab Fraction of reads in cDNA identified as HDR. \cr
  \strong{wt_rate_gDNA} \tab Fraction of reads in gDNA identified as WT \cr
  \strong{wt_rate_cDNA} \tab Fraction of reads in cDNA identified as WT. \cr
  \strong{del_rate_gDNA} \tab Fraction of reads in gDNA identified as having a CRISPR deletion. \cr
  \strong{del_rate_cDNA} \tab Fraction of reads in cDNA identified as having a CRISPR deletion. \cr
  \strong{hdr_effect} \tab HDR allele: Estimated effect size - the amount by which the HDR:WT ratio differs in cDNA relative to gDNA. \cr
  \strong{hdr_effect_sd} \tab HDR allele: Standard deviation of the estimated effect size. \cr
  \strong{hdr_effect_confint_lo} \tab HDR allele: Lower bound of the 95\% confidence interval for the effect size. \cr
  \strong{hdr_effect_confint_hi} \tab HDR allele: Upper bound of the 95\% confidence interval for the effect size. \cr
  \strong{hdr_df_estimated} \tab THDR allele: he degrees of freedom used in the unequal variances t-test, which is estimated from the data. \cr
  \strong{hdr_pval} \tab HDR allele: A p value from the unequal variants t-test. \cr
}

The fields above beginning with `hdr_` give statistics relating to the HDR allele.
There are 6 equivalent fields that begin with `del_` which relate to all deletions.
There are also 6 equivalent fields that begin with `del_window_`, which relate to
deletions that are contained within a window around the `highlight_site` (defined in
the region input), the extent of which is determined by `crispr_del_window` parameter.
}

\section{Additional fields}{
 The deletion_analysis result object additionally has the following fields:
\tabular{rl}{
  \strong{replicate_qc} \tab A data.frame of summary information for each replicate,
  with counts of reads in different categories, editing rates per replicate, and
  quality control summary information. \cr
  \strong{replicate_alleles} \tab A data.frame of summary information for each unique
  allele in each replicate. \cr
  \strong{region_alleles} \tab A data.frame of summary information for each unique
  allele averaged across all replicates. \cr
  \strong{replicate_allele_fractions} \tab A data.frame of summary information for
  the top 20 alleles across all replicates, which is used for replicate quality control. \cr
  \strong{allele_effect} \tab A data.frame of GenIE effect size estimates (difference
  in allele:WT ratio in cDNA vs. gDNA) for each unique allele with sufficient reads
  across replicates. \cr
  \strong{site_profiles} \tab A data.frame which indicates read counts for each combination of observed
  nucleotides at each specified position in the wt_allele_profile, separately for each replicate. For
  example, if 1 position is defined, then there will be up to 5 unique combinations in each replicate,
  accounting for A, C, G, T, or * (deletion). If 2 positions are defined, there are up to 25 combinations.
  This can be useful to check whether the fraction of WT or edited reads is as expected. In particular,
  if a heterozygous site is targeted, and only one of the haplotype alleles actually gets edited, then
  you expect depletion of only one of the nucleotides at the SNP site. \cr
  \strong{mismatch_profiles} \tab A data.frame which contains every unique read profile,
  including mismatches, for reads considered as WT or HDR. This can help to identify if
  there are an abundance of "HDR" or "WT" reads with mismatches at a particular position,
  or an excess of mismatches in general. \cr
  \strong{replicate_list} \tab Internal data for each replicate. Not likely to be of interest to the user. \cr
  \strong{type} \tab  \cr
}
}

\examples{
# Note: to run an analysis you need BAM files from a GenIE experiment.
# An example is available for download using download_example().

\donttest{
download_example(dir = "~/genie_example", name = "MUL1")
# Data are downloaded and we can run an rgenie analysis
setwd("~/genie_example/MUL1/")
regions = readr::read_tsv("mul1.genie_regions.tsv")
replicates = readr::read_tsv("mul1.genie_replicates.tsv")
delresults = deletion_analysis(regions, replicates)
deletion_plots(delresults[[1]])
}
}
\seealso{
\code{\link{grep_analysis}}

\code{\link{deletion_plots}}

\code{\link{deletion_summary_plot}}

\code{\link{experiment_summary_plot}}

\code{\link{deletion_alleles_plot}}

\code{\link{deletion_profile_plot}}

\code{\link{replicate_summary_plot}}

\code{\link{replicate_qc_plot}}

\code{\link{allele_effect_plot}}

\code{\link{get_variance_components}}

\code{\link{variance_components_plot}}

\code{\link{power_plots}}

\code{\link{bind_results}}
}
