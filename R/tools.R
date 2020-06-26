#' Get test counts. This is for debugging purposes only.
#' @export
#' @import readr dplyr
get_test_counts <- function() {
  dplyr::left_join(
    # HFD
    system.file(
      "extdata/salmon_merged_gene_counts.csv.gz",
      package = "bruening.2020.04.marta.bactrap"
    ) %>%
      readr::read_csv(),
    # CD
    system.file(
      "extdata/salmon_merged_gene_counts.csv.gz",
      package = "bruening.2019.12.marta.bactrap"
    ) %>%
      readr::read_csv(),
    by = "gene_id"
  ) %>%
    # Only include wt samples
    dplyr::select("gene_id", tidyr::contains("wt")) %>%
    dplyr::rename(
      wt_input_ncd_1 = INPUT_WT_1,
      wt_input_ncd_2 = INPUT_WT_2,
      wt_input_ncd_3 = INPUT_WT_3,
      wt_ip_ncd_1 = IP_WT_1,
      wt_ip_ncd_2 = IP_WT_2,
      wt_ip_ncd_3 = IP_WT_3,
      wt_input_hfd_1 = wt_input_879_1,
      wt_input_hfd_2 = wt_input_879_2,
      wt_input_hfd_3 = wt_input_879_3,
      wt_input_hfd_4 = wt_input_879_4,
      wt_ip_hfd_1 = wt_ip_879_1,
      wt_ip_hfd_2 = wt_ip_879_2,
      wt_ip_hfd_3 = wt_ip_879_3,
      wt_ip_hfd_4 = wt_ip_879_4
    ) %>%
    # Round values to comply with DESeq2 requirements
    dplyr::mutate_if(is.numeric, round)
}