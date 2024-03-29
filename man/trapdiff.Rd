% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools.R
\name{trapdiff}
\alias{trapdiff}
\title{Run trapdiff and save report file}
\usage{
trapdiff(
  counts,
  path_config_json,
  out_path,
  ensembl_version = 102,
  biotypes_filter = "",
  tpms_min = 0,
  padj_cutoff = 0.05,
  save_rds = TRUE,
  save_excel = TRUE,
  save_figures = TRUE,
  deseq_split_size_factors = FALSE,
  splits = "treatment",
  filter_regex = "^mt-"
)
}
\arguments{
\item{counts}{Counts per sample as dataframe}

\item{path_config_json}{Configuration file telling trapdiff group associations}

\item{out_path}{Output folder for report and intermediate files}

\item{ensembl_version}{Ensembl version used for attaching biomart variables}

\item{biotypes_filter}{Filter genes for these biotypes. Put empty string for no filtering.}

\item{tpms_min}{A gene is included when at least one sample contains tpms_min reads}

\item{padj_cutoff}{Minimum padj value for a significantly differentially expressed gene}

\item{save_rds}{Save rds file of results}

\item{save_excel}{Save Excel sheets of results}

\item{save_figures}{Save output figures as png and pdf in out_path/figures}

\item{deseq_split_size_factors}{Create sizeFactor estimates separately for each split condition}

\item{splits}{If deseq_split_size_factors is set to true, sizeFactors are calculated for these groups separately}

\item{filter_regex}{Manually regex out genes. Defaults to all mitochondrial genes starting with "mt-"}
}
\description{
Run trapdiff and save report file
}
