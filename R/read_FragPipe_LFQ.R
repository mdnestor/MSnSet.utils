#' Reading MSFragger-generated LFQ-based MSstats from a file path as MSnSet object
#'
#' @description Function has only been tested with label-free intensity-based
#'   quantification data. MSstats.csv is
#'   an optional output file which needs to be specified in FP settings.
#'
#' @param path character; File path to the FragPipe-generated MSstats.csv file
#'
#' @return (MSnSet) MSnSet object of MSFragger LFQ results
#'
#' @importFrom MSnbase MSnSet
#' @importFrom data.table fread
#' @importFrom tidyr pivot_wider
#' @importFrom tibble column_to_rownames
#' @importFrom dplyr %>% select filter distinct relocate everything mutate
#'
#' @examples
#' file_path <- "C:/Users/fakeusr222/Desktop/MSF_LFQ_job/MSstats.csv"
#'   msnset <- read_FragPipe_LFQ(file_path)
#'   show(msnset)
#'
#' @export read_FragPipe_LFQ


read_FragPipe_LFQ <- function(path = NULL)
{
  path_to_file <- path

  if (!file.exists(path_to_file)) {
    stop(sprintf("MSstats.csv file not found in folder: %s", dirname(path_to_file)))
  }

  df <- fread(file = path_to_file, showProgress = FALSE, data.table = FALSE) %>%
    filter(!is.na(Intensity)) %>%
    # May add charge col later
    select(ProteinName, PeptideSequence, Run, Intensity) %>%
    mutate(featureName = paste0(ProteinName, "@", PeptideSequence)) %>%
    relocate(featureName, .before = everything())

  # Will sum intensity of unique features.
  x_data <- df %>%
    pivot_wider(id_cols = "featureName",
                names_from = "Run",
                values_from = "Intensity",
                values_fn = sum) %>%
    as.data.frame() %>%
    column_to_rownames(var = "featureName") %>%
    as.matrix()

  f_data <- df %>%
    distinct(featureName, ProteinName, PeptideSequence) %>%
    `rownames<-`(.[["featureName"]])

  p_data <- df %>%
    distinct(Run) %>%
    `rownames<-`(.[["Run"]])

  x_data <- x_data[rownames(f_data), rownames(p_data)]

  m <- MSnSet(exprs = x_data, fData = f_data, pData = p_data)

  return(m)
}




utils::globalVariables(
  c("ProteinName", "PeptideSequence", "Run", "Intensity", ".", "featureName")
)

