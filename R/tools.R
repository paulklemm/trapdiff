#' Run trapdiff and save report file
#' 
#' @export
#' @import magrittr rmarkdown
#' 
#' @param counts Counts per sample as dataframe
#' @param path_config_json Configuration file telling trapdiff group associations
#' @param out_path Output folder for report and intermediate files
#' @param biotypes_filter Filter genes for these biotypes
#' @param tpms_min A gene is included when at least one sample contains tpms_min reads
#' @param padj_cutoff Minimum padj value for a significantly differentially expressed gene
#' @param save_rds Save rds file of results
#' @param save_excel Save Excel sheets of results
#' @param save_figures Save output figures as png and pdf in out_path/figures
#' @param ensembl_version Ensembl version used for attaching biomart variables
#' @param deseq_split_size_factors Create sizeFactor estimates separately for each split condition
#' @param splits If deseq_split_size_factors is set to true, sizeFactors are calculated for these groups separately
#' @param filter_regex Manually regex out genes. Defaults to all mitochondrial genes starting with "mt-"
trapdiff <- function(
  counts,
  path_config_json,
  out_path,
  ensembl_version = 102,
  biotypes_filter = "protein_coding",
  tpms_min = 1,
  padj_cutoff = 0.05,
  save_rds = TRUE,
  save_excel = TRUE,
  save_figures = TRUE,
  deseq_split_size_factors = FALSE,
  splits = "treatment",
  filter_regex = "^mt-"
) {
  # Be sure the output path exists
  if (!dir.exists(out_path)) {
    paste0(
      "The output path ",
      out_path,
      " does not exist. I will create it now."
    ) %>%
    message()
    dir.create(
      path = out_path,
      showWarnings = TRUE,
      recursive = TRUE
    )
  }

  # Render command with all parameters
  rmarkdown::render(
    system.file("rmd/trapdiff/trapdiff.Rmd", package = "trapdiff"),
    params = list(
      path_out = out_path,
      path_config_json = path_config_json,
      biotypes_filter = biotypes_filter,
      tpms_min = tpms_min,
      padj_cutoff = padj_cutoff,
      save_rds = save_rds,
      save_excel = save_excel,
      save_figures = save_figures,
      ensembl_version = ensembl_version,
      counts = counts,
      deseq_split_size_factors = deseq_split_size_factors,
      splits = splits,
      filter_regex = filter_regex
    ),
    # Change the intermediate path to the output to avoid write access errors
    intermediates_dir = out_path,
    knit_root_dir = out_path,
    # Clean intermediate files created during rendering.
    clean = TRUE,
    output_dir = out_path,
    output_options = list(
      self_contained = TRUE
    )
  )
}