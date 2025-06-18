#' Remove Rows with Missing Units or Medium Values from a DataFrame
#'
#' This function takes a dataframe and removes rows where either the 'Unit' or 'Med'
#' column has missing (NA) values. It also converts specified columns to numeric type 
#' if they exist in the dataframe.
#'
#' @param df_param A dataframe with columns 'Unit', 'Med', and others. 
#'                It is expected to contain model parameters.
#' @return A cleaned dataframe with rows having missing 'Unit' or 'Med' values removed 
#'         and specified columns converted to numeric type.
#' @import dplyr
#' @export
remove_param_without_unit <- function(df_param){
  
  # Checking for the presence of required columns
  required_cols <- c("Unit", "Med", "Param")
  if (!all(required_cols %in% names(df_param))) {
    stop("Dataframe does not have the required columns")
  }
  
  # Identifying rows with missing 'Unit' or 'Med' values
  v_row_missing_values <- unique(c(which(is.na(df_param$Unit)), which(is.na(df_param$Med))))
  
  # Removing rows with missing values
  df_param_clean <- df_param[-v_row_missing_values, ]
  
  # Columns specified for conversion to numeric type
  specified_cols <- c("Med", "Lo_alpha", "Hi_beta")
  
  # Finding columns present in df_param to convert
  cols_to_convert <- intersect(specified_cols, names(df_param))
  
  
  # Converting identified columns to numeric type
  df_param_clean[cols_to_convert] <- lapply(df_param_clean[cols_to_convert], as.numeric)
  
  return(df_param_clean) # Returning the cleaned dataframe
}



