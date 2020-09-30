## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----install, eval=FALSE, echo=TRUE-------------------------------------------
#  devtools::install_github(repo = "jeremy37/rgenie")

## ---- message=FALSE-----------------------------------------------------------
library(rgenie)
library(readr, quietly = TRUE) # Not required, but we use it below
library(magrittr, quietly = TRUE)

## ----download_example_x, echo=FALSE, eval=FALSE-------------------------------
#  # Note that the get the vignette to pass CRAN's checks, we can't download data
#  # as part of the vignette. What we have to do is to provide example data that
#  # has the results objects already. But we would like to display to the user
#  # how to do everything from the downloaded example. To do this, I echo the
#  # code chunks for downloading the example, but don't run them. The results
#  # should already be loaded as part of the package data.

## ----download_example, eval=FALSE, echo=TRUE----------------------------------
#  rgenie::download_example(dir = "~/genie_example", name = "MUL1")

## ----regions, eval=FALSE, echo=TRUE-------------------------------------------
#  mul1_regions = readr::read_tsv("~/genie_example/MUL1/mul1.genie_regions.tsv")
#  print(mul1_regions)

## ---- echo=FALSE--------------------------------------------------------------
print(mul1_regions)

## -----------------------------------------------------------------------------
print(mul1_regions$wt_allele_profile)
print(mul1_regions$hdr_allele_profile)
print(mul1_regions$ref_sequence)

## ----replicates, eval=FALSE---------------------------------------------------
#  mul1_replicates = readr::read_tsv("~/genie_example/MUL1/mul1.genie_replicates.tsv")
#  print(mul1_replicates)

## ---- echo=FALSE--------------------------------------------------------------
print(mul1_replicates)

## ----grep_example, eval=FALSE, echo=TRUE--------------------------------------
#  setwd("~/genie_example/MUL1/")
#  mul1_grep_results = grep_analysis(mul1_regions,
#                                    mul1_replicates,
#                                    required_match_left = 10,
#                                    required_match_right = 10,
#                                    min_mapq = 0,
#                                    quiet = TRUE)

## -----------------------------------------------------------------------------
length(mul1_grep_results)
grep_result = mul1_grep_results[[1]]
names(grep_result)

## -----------------------------------------------------------------------------
grep_result$region_stats

## -----------------------------------------------------------------------------
grep_result$region_stats$effect
grep_result$region_stats$pval

## ----fig.width=7, fig.height=6------------------------------------------------
rgenie::grep_summary_plot(grep_result)

## -----------------------------------------------------------------------------
mul1_replicates = mul1_replicates %>% dplyr::filter(replicate != "c1.1")

## ---- eval=FALSE, echo=TRUE---------------------------------------------------
#  setwd("~/genie_example/MUL1/")
#  mul1_del_results = rgenie::deletion_analysis(mul1_regions,
#                                               mul1_replicates,
#                                               required_match_left = 10,
#                                               required_match_right = 10,
#                                               crispr_del_window = 100,
#                                               min_mapq = 0,
#                                               max_mismatch_frac = 0.05,
#                                               min_aligned_bases = 50,
#                                               exclude_multiple_deletions = FALSE,
#                                               exclude_nonspanning_reads = TRUE,
#                                               allele_profile = FALSE,
#                                               del_span_start = -20,
#                                               del_span_end = 20,
#                                               quiet = TRUE)

## -----------------------------------------------------------------------------
del_result = mul1_del_results[[1]]
names(del_result)

## -----------------------------------------------------------------------------
tbl = tibble::tibble(Name = names(del_result$region_stats), Value = as.character(del_result$region_stats))
knitr::kable(tbl)

## ----fig.width=7, fig.height=6------------------------------------------------
rgenie::deletion_summary_plot(del_result)
del_result$replicate_stats

## -----------------------------------------------------------------------------
del_result$replicate_stats

## ----fig.width=7, fig.height=6------------------------------------------------
rgenie::deletion_alleles_plot(del_result,
                              viewing_window = 40,
                              color_by="window")

## ----fig.width=7, fig.height=6------------------------------------------------
rgenie::deletion_profile_plot(del_result,
                              viewing_window = 40)

## ----fig.width=7, fig.height=6------------------------------------------------
rgenie::replicate_summary_plot(del_result,
                               outlier_threshold = NA)

## ----fig.width=7, fig.height=6------------------------------------------------
rgenie::replicate_qc_plot(del_result,
                          outlier_threshold = NA)

## ----fig.width=7, fig.height=6------------------------------------------------
rgenie::allele_effect_plot(del_result,
                           viewing_window = 40,
                           max_alleles = 40)

## ----fig.width=7, fig.height=6------------------------------------------------
all_plots = rgenie::deletion_plots(del_result,
                                   opts = genie_plot_options(),
                                   variance_components_plot = FALSE,
                                   power_plots = FALSE)

names(all_plots)

## ----eval=FALSE, echo=TRUE----------------------------------------------------
#  pdf("genie_plots.pdf", width=7, height=6)
#  print(all_plots)
#  dev.off()

## -----------------------------------------------------------------------------
plot_opts = genie_plot_options()
print(plot_opts)

## ----fig.width=7, fig.height=6------------------------------------------------
rgenie::experiment_summary_plot(mul1_grep_results, mul1_del_results)

## ----fig.width=7, fig.height=6------------------------------------------------
grep_tables = rgenie::bind_results(mul1_grep_results)
names(grep_tables)

del_tables = rgenie::bind_results(mul1_del_results)
names(del_tables)

