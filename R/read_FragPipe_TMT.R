#' Reading MSFragger-generated tmt-report files from a file path as MSnSet object
#'
#' @description Function has only been tested with TMT intensity-based
#'   quantification data. Desired tmt-report output files (e.g., "ratio_multi-site_MD.tsv")
#'   must be properly selected in the FP settings.
#'
#' @param path character; File path to the desired FragPipe-generated tmt-report file.
#'   Any tmt-report file may be used.
#'
#' @return (MSnSet) MSnSet object of MSFragger TMT results
#'
#' @importFrom MSnbase MSnSet
#' @importFrom data.table fread
#' @importFrom tibble column_to_rownames
#' @importFrom dplyr %>% select mutate
#'
#' @examples
#' file_path <- "C:/Users/fakeusr222/Desktop/MSF_TMT_job/ratio_multi-site_MD.tsv"
#'   msnset <- read_FragPipe_TMT(file_path)
#'   show(msnset)
#'
#' @export read_FragPipe_TMT

read_FragPipe_TMT <- function(path = NULL)
{

  path_to_file <- path

  if (!file.exists(path_to_file)) {
    stop(sprintf("file not found in folder: %s", dirname(path_to_file)))
  }

  df <- fread(file = path_to_file, showProgress = FALSE, data.table = FALSE)


  # make featureNames
  if (grepl("multi-site|single-site|peptide", basename(path_to_file))) {
    df <- df %>%
      mutate(rowname = paste(Gene, ProteinID, Peptide, sep = "|"))
  } else if (grepl("gene", basename(path_to_file))) {
    df <- df %>%
      mutate(rowname = paste(Index, ProteinID, sep = "|"))
  } else if (grepl("protein", basename(path_to_file))) {
    df <- df %>%
      mutate(rowname = paste(Gene, Index, sep = "|"))
  }

  df <- df %>%
    mutate(featureName = rowname, .before = colnames(.)[[1]]) %>%
    column_to_rownames(var = "rowname")

  x_data <- df %>%
    select(-c(colnames(.)[[1]]:ReferenceIntensity)) %>%
    as.matrix()

  f_data <- df %>%
    select(c(colnames(.)[[1]]:ReferenceIntensity))

  m <- MSnSet(exprs = x_data, fData = f_data)

  return(m)
}




utils::globalVariables(
  c(".", "featureName")
)

