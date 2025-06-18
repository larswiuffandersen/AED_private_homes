# Function to validate parameters for each distribution
validate_params <- function(param) {
  errors <- c()  # Initialize an empty vector to collect errors
  
  # Check if required columns are present
  required_columns <- c("Distribution", "Treatment", "Param")
  missing_columns <- setdiff(required_columns, colnames(param))
  if (length(missing_columns) > 0) {
    stop(paste("Missing required columns:", paste(missing_columns, collapse = ", ")))
  }
  
  # Loop through each row in the parameters dataframe
  for (i in 1:nrow(param)) {
    dist <- param[i, "Distribution"]
    treatment <- param[i, "Treatment"]
    param_name <- param[i, "Param"]
    
    # Check for missing values in essential columns
    if (is.na(dist) || dist == "") {
      errors <- c(errors, paste("Missing Distribution in row:", i, "param:", param_name))
      next
    }
    if (is.na(treatment) || treatment == "") {
      errors <- c(errors, paste("Missing Treatment in row:", i, "param:", param_name))
    }
    if (is.na(param_name) || param_name == "") {
      errors <- c(errors, paste("Missing Param in row:", i))
    }
    
    # Handle unknown distribution types
    if (!dist %in% c("Triangular", "Normal", "Lognormal", "Beta", "Uniform", "Gamma", "NA")) {
      errors <- c(errors, paste("Unknown distribution type for treatment:", treatment, "param:", param_name))
      next
    }
    
    # Validate parameters based on the distribution type
    if (dist == "Triangular") {
      if (is.na(param[i, "Lo_alpha"]) || is.na(param[i, "Med"]) || is.na(param[i, "Hi_beta"])) {
        errors <- c(errors, paste("Missing parameters for Triangular distribution in treatment:", treatment, "param:", param_name))
      } else if (param[i, "Lo_alpha"] > param[i, "Med"] || param[i, "Med"] > param[i, "Hi_beta"]) {
        errors <- c(errors, paste("Invalid parameters for Triangular distribution in treatment:", treatment, "param:", param_name, "- Ensure Lo_alpha <= Med <= Hi_beta"))
      }
    } else if (dist == "Normal") {
      if (is.na(param[i, "Med"]) || is.na(param[i, "Hi_beta"])) {
        errors <- c(errors, paste("Missing parameters for Normal distribution in treatment:", treatment, "param:", param_name))
      } else if (param[i, "Hi_beta"] < param[i, "Med"]) {
        errors <- c(errors, paste("Invalid parameters for Normal distribution in treatment:", treatment, "param:", param_name, "- Ensure Hi_beta >= Med"))
      }
    } else if (dist == "Lognormal") {
      if (is.na(param[i, "Med"]) || is.na(param[i, "Hi_beta"])) {
        errors <- c(errors, paste("Missing parameters for Lognormal distribution in treatment:", treatment, "param:", param_name))
      }
    } else if (dist == "Beta") {
      if (is.na(param[i, "Lo_alpha"]) || is.na(param[i, "Hi_beta"])) {
        errors <- c(errors, paste("Missing parameters for Beta distribution in treatment:", treatment, "param:", param_name))
      }
    } else if (dist == "Gamma") {
      if (is.na(param[i, "Lo_alpha"]) || is.na(param[i, "Hi_beta"])) {
        errors <- c(errors, paste("Missing parameters for Beta distribution in treatment:", treatment, "param:", param_name))
      }
    } else if (dist == "Uniform") {
      if (is.na(param[i, "Lo_alpha"]) || is.na(param[i, "Hi_beta"])) {
        errors <- c(errors, paste("Missing parameters for Uniform distribution in treatment:", treatment, "param:", param_name))
      } else if (param[i, "Lo_alpha"] > param[i, "Hi_beta"]) {
        errors <- c(errors, paste("Invalid parameters for Uniform distribution in treatment:", treatment, "param:", param_name, "- Ensure Lo_alpha <= Hi_beta"))
      }
    } else if (dist == "NA") {
      if (is.na(param[i, "Med"])) {
        errors <- c(errors, paste("Missing parameter Med for NA distribution in treatment:", treatment, "param:", param_name))
      }
    }
  }
  
  if (length(errors) == 0) {
    return("no errors")
  } else {
    return(errors)  # Return collected errors
  }
}



# The code used where the PSA function for lognormal was calculated using the mean and high


#' Generate data for the iterations of  probabilistic sensitivity analysis 
#' 
#' \code{make_psa_df} is used to sample the values for each iteration of a probabilistic sensitivity analysis (PSA) and stores it in a dataframe
#' @param param
#' @param n_iter
#' @param seed
#' @keywords Probabilistic sensitivity analysis, PSA, dataframe
#' @section Details:
#' \code{make_psa_df} adds the values for each PSA iteration to the dataframe
#' @return  param_psa A dataframe with the input values of all parameters of each PSA run. 
#' 
# Make a function to make the PSA data set 
# Function to create PSA dataset
make_psa_df <- function(param, n_iter, seed = 123, validate = TRUE) {
  set.seed(seed)  # Set seed for reproducibility
  
  # Validate parameters if the validate option is TRUE
  if (validate) {
    validation_result <- validate_params(param)
    if (is.character(validation_result) && validation_result == "no errors") {
      message("Validation successful: no errors")
    } else {
      stop(paste(validation_result, collapse = "\n"))  # Stop execution and print errors if any
    }
  }
  
  # Draw samples for PSA
  param_psa <- as.data.frame(lapply(param, rep, n_iter))  # Replicate each parameter n_iter times
  param_psa$psa_est <- NA  # Initialize psa_est column with NA
  param_psa$iter <- rep(1:n_iter, each = nrow(param))  # Create iteration column
  
  # Loop through each parameter row to draw samples
  for (i in 1:nrow(param)) {    
    distribution  <- param[i, "Distribution"]
    treatment     <- param[i, "Treatment"]
    param_name    <- param[i, "Param"]
    
    if (distribution == "Triangular") {
      param_psa$psa_est[param_psa$Treatment == treatment & param_psa$Param == param_name] <- with(param[i, ], 
                                                                                                  rtriangle(n = n_iter,
                                                                                                            a = Lo_alpha, 
                                                                                                            b = Hi_beta, 
                                                                                                            c = Med))
    } else if (distribution == "Normal") {
      param_psa$psa_est[param_psa$Treatment == treatment & param_psa$Param == param_name] <- with(param[i, ], 
                                                                                                  rnorm(n = n_iter,
                                                                                                        mean = Med,
                                                                                                        sd = (Hi_beta - Med) / 1.96))
    } else if (distribution == "Lognormal") {
      param_psa$psa_est[param_psa$Treatment == treatment & param_psa$Param == param_name] <- with(param[i, ], 
                                                                                                      rlnorm(n = n_iter,
                                                                                                            meanlog = Lo_alpha,
                                                                                                            sdlog = Hi_beta))
    } else if (distribution == "Beta") {
      param_psa$psa_est[param_psa$Treatment == treatment & param_psa$Param == param_name] <- with(param[i, ],
                                                                                                  rbeta(n = n_iter,
                                                                                                        shape1 = Lo_alpha,
                                                                                                        shape2 = Hi_beta)) 
    } else if (distribution == "Gamma") {
      param_psa$psa_est[param_psa$Treatment == treatment & param_psa$Param == param_name] <- with(param[i, ],
                                                                                                  rgamma(n = n_iter,
                                                                                                        shape = Lo_alpha,
                                                                                                        scale = Hi_beta)) 
    } else if (distribution == "Uniform") {
      param_psa$psa_est[param_psa$Treatment == treatment & param_psa$Param == param_name] <- with(param[i, ],
                                                                                                  runif(n = n_iter,
                                                                                                        min = Lo_alpha,
                                                                                                        max = Hi_beta)) 
    } else if (distribution == "NA") {
      param_psa$psa_est[param_psa$Treatment == treatment & param_psa$Param == param_name] <- with(param[i, ],
                                                                                                  rep(x = Med, n_iter)) 
    }
  }
  
  return(param_psa)  # Return the parameter values for the PSA runs
}


#' Generate probabilistic sensitivity analysis dataframe in DARTH style
#' 
#' \code{gen_psa} is used to compute the probabilistic sensitivity analysis (PSA) dataframe in DARTH style.
#' @param df_param_psa
#' @keywords Probabilistic sensitivity analysis, PSA, DARTH
#' @section Details:
#' \code{gen_psa} reorganises the values from the dataframe that is entered in the function to a data frame with values of each paramters of each PSA run.
#' @return  df_psa_input A data frame with the input values of all parameters of each PSA run. Dimension of the data frame are number of PSA simulations * parameters 
#' 

gen_psa <- function(df_param_psa){
  # Argument
  ## df_param_psa: A dataframe with the parameter values for the PSA
  # Return
  ## df_psa_input: The input values of each PSA iteration in a dataframe of size n_sim x n_parameters
  
  param_names <- unique(df_param_psa$Param)
  n_sim <- length(unique(df_param_psa$iter))
  
  # Create a matrix to store the results 
  df_psa_input <- matrix(data = NA, nrow = n_sim, ncol = length(param_names),
                         dimnames = list(c(1:n_sim), param_names))
  
  for(p in 1:length(param_names)){
    df_psa_input[, p] <- df_param_psa$psa_est[df_param_psa$Param == param_names[p]]
  }
  
  df_psa_input <- as.data.frame(df_psa_input) # Make a dataframe from the matrix 
  
  return(df_psa_input) # Return the psa input dataframe
}


