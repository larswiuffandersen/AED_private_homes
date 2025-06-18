#-------------------------------------------------------------------#
#### Function to simulate the cohort model  ####
#-------------------------------------------------------------------#
#' This code runs the simulation of the cohort for the markov model 
#'
#' \code{run_markov} runs the simulation for the markov model 
#' 
#' @param l_list  List with model parameters 
#' @param n_wtp  Willingness to pay, default is 50000
#' @param return_all Specify if all model output should be returned or just the cost-effectiveness (ce) output, is ce output by default
#' @param verbose Specifies 
#' @return 
#' A dataframe with the results of the cost-effectiveness analysis
#' @import 
#' @export
#' 
#'

calculate_cea_output_private_AED <- function(l_list, n_wtp = 50000, return_all = FALSE, verbose = FALSE){
  with(as.list(l_list), {
    # Strategy names
    v_names_str <- c("AED", "noAED") 
    n_str       <- length(v_names_str)     # Number of strategies
    
    # Markov model parameters
    v_names_states <- c("CPC1","CPC2", "CPC3","CPC4", "D")  # the 6 states of the model: Cardiac arrest, functional neurological status 1 through 4, death
    n_states       <- length(v_names_states)  # number of health states 
    v_names_cycles  <- paste("cycle", 0:n_cycles)    # cycle names
    
    #Probability of shockable rhythm and AED used
    # Risk ratios are reported in the data, they can be applied to the probabilities 
    p_shockable_AED        <- rr_AED_group_shockable * p_shockable_noAED  	  # Convert the probability of shockable rhythm without AED to shockable rhythm with AED
    
    # Probability to survive or die after cardiac arrest with or without AED
    # Calculate the rate of death after shockable cardiac arrest with AED use, based on heart rhythm.
    p_survive_shockable_AED      <- rr_survive_AED_group_shockable * p_survive_shockable_noAED   # Convert the probability of death without AED to probability of death with AED
    
    # Calculate the probability of death after nonshockable cardiac arrest with AED use, based on heart rhythm (hr).
    p_survive_nonshockable_AED       <- rr_survive_AED_group_nonshockable * p_survive_nonshockable_noAED # Convert the probability of death without AED to probability of death with AED
    
    # Calculate probabilities for death from probabilities to survive
    p_arrest_D_shockable_AED      <- 1 - p_survive_shockable_AED      # probability to die after shockable arrest with AED
    p_arrest_D_shockable_noAED    <- 1 - p_survive_shockable_noAED    # probability to die after shockable arrest without AED
    p_arrest_D_nonshockable_AED   <- 1 - p_survive_nonshockable_AED   # probability to die after nonshockable arrest with AED
    p_arrest_D_nonshockable_noAED <- 1 - p_survive_nonshockable_noAED # probability to die after nonshockable arrest without AED
    
    # Create a variable to store probabilities for chance nodes with multiple branches
    # Probabilities by outcome (CPC)
    p_CPC_shockable_AED   <- (p_CPC1_shockable_AED   + p_CPC2_shockable_AED   + p_CPC3_shockable_AED   + p_CPC4_shockable_AED)
    p_CPC_shockable_noAED <- (p_CPC1_shockable_noAED + p_CPC2_shockable_noAED + p_CPC3_shockable_noAED + p_CPC4_shockable_noAED)
    p_CPC_nonshockable    <- (p_CPC1_nonshockable    + p_CPC2_nonshockable    + p_CPC3_nonshockable    + p_CPC4_nonshockable) 
    
    # Calculate probability of shockable rhythm by AED use
    p_nonshockable_AED   <- 1 - p_shockable_AED   # the probability of having a nonshockable event in the AED arm
    p_nonshockable_noAED <- 1 - p_shockable_noAED # the probability of having a nonshockable event in the no AED arm
    
  
    #Probability AED not used
    p_AED_not_used <- 1- p_AED_used
    
    #Probability of death after first year by CPC
    v_CPC1_D <- c(p_CPC1_D_Y1 * p_multiplier, rep(p_CPC1_D * p_multiplier, each = n_cycles - 1))
    
    v_CPC2_D <- c(p_CPC2_D_Y1 * p_multiplier, rep(p_CPC2_D * p_multiplier, each = n_cycles - 1))
    
    v_CPC3_D <- c(p_CPC3_D_Y1 * p_multiplier, rep(p_CPC3_D * p_multiplier, each = n_cycles - 1))
    
    v_CPC4_D <- c(p_CPC4_D_Y1 * p_multiplier, rep(p_CPC4_D * p_multiplier, each = n_cycles - 1))
    
    
    v_m_init_AED <- c("CPC1" = p_AED_used     * p_shockable_AED      * p_survive_shockable_AED      * p_CPC1_shockable_AED /   p_CPC_shockable_AED +   # AED used, shockable, survive, CPC1
                               p_AED_used     * p_nonshockable_AED   * p_survive_nonshockable_AED   * p_CPC1_nonshockable /    p_CPC_nonshockable +    # AED used, nonshockable, survive, CPC1
                               p_AED_not_used * p_shockable_noAED    * p_survive_shockable_noAED    * p_CPC1_shockable_noAED / p_CPC_shockable_noAED + # AED not used, shockable, survive, CPC1
                               p_AED_not_used * p_nonshockable_noAED * p_survive_nonshockable_noAED * p_CPC1_nonshockable /    p_CPC_nonshockable,     # AED not used, non shockable, survive, CPC1 
                      
                      "CPC2" = p_AED_used     * p_shockable_AED      * p_survive_shockable_AED      * p_CPC2_shockable_AED/   p_CPC_shockable_AED +    # AED used, shockable, survive, CPC2
                               p_AED_used     * p_nonshockable_AED   * p_survive_nonshockable_AED   * p_CPC2_nonshockable/    p_CPC_nonshockable +     # AED used, nonshockable, survive, CPC2
                               p_AED_not_used * p_shockable_noAED    * p_survive_shockable_noAED    * p_CPC2_shockable_noAED/ p_CPC_shockable_noAED +  # AED not used, shockable, survive, CPC2
                               p_AED_not_used * p_nonshockable_noAED * p_survive_nonshockable_noAED * p_CPC2_nonshockable/    p_CPC_nonshockable,      # AED not used, non shockable, survive, CPC2
                      
                      "CPC3" = p_AED_used     * p_shockable_AED      * p_survive_shockable_AED      * p_CPC3_shockable_AED/   p_CPC_shockable_AED +    # AED used, shockable, survive, CPC3
                               p_AED_used     * p_nonshockable_AED   * p_survive_nonshockable_AED   * p_CPC3_nonshockable/    p_CPC_nonshockable +     # AED used, nonshockable, survive, CPC3
                               p_AED_not_used * p_shockable_noAED    * p_survive_shockable_noAED    * p_CPC3_shockable_noAED/ p_CPC_shockable_noAED +  # AED not used, shockable, survive, CPC3
                               p_AED_not_used * p_nonshockable_noAED * p_survive_nonshockable_noAED * p_CPC3_nonshockable/    p_CPC_nonshockable,      # AED not used, non shockable, survive, CPC3, 
                      
                      "CPC4" = p_AED_used     * p_shockable_AED      * p_survive_shockable_AED      * p_CPC4_shockable_AED/   p_CPC_shockable_AED +    # AED used, shockable, survive, CPC4
                               p_AED_used     * p_nonshockable_AED   * p_survive_nonshockable_AED   * p_CPC4_nonshockable/    p_CPC_nonshockable +     # AED used, nonshockable, survive, CPC4
                               p_AED_not_used * p_shockable_noAED    * p_survive_shockable_noAED    * p_CPC4_shockable_noAED/ p_CPC_shockable_noAED +  # AED not used, shockable, survive, CPC4
                               p_AED_not_used * p_nonshockable_noAED * p_survive_nonshockable_noAED * p_CPC4_nonshockable/    p_CPC_nonshockable,      # AED not used, non shockable, survive, CPC4
                      
                      "D"    = p_AED_used     * p_shockable_AED      * p_arrest_D_shockable_AED +     # AED used, shockable, die
                               p_AED_used     * p_nonshockable_AED   * p_arrest_D_nonshockable_AED +  # AED used, nonshockable, die
                               p_AED_not_used * p_shockable_noAED    * p_arrest_D_shockable_noAED +   # AED not used, shockable, die
                               p_AED_not_used * p_nonshockable_noAED * p_arrest_D_nonshockable_noAED) # AED not used, non shockable, die  
    
    
    #No AED strategy
    v_m_init_noAED <- c("CPC1" = p_shockable_noAED * p_survive_shockable_noAED    * p_CPC1_shockable_noAED/p_CPC_shockable_noAED +  # shockable, survive, CPC1
                              p_nonshockable_noAED * p_survive_nonshockable_noAED * p_CPC1_nonshockable/p_CPC_nonshockable,         # non shockable, survive, CPC1            
                        
                        "CPC2" = p_shockable_noAED * p_survive_shockable_noAED    * p_CPC2_shockable_noAED/p_CPC_shockable_noAED +  # shockable, survive, CPC2
                              p_nonshockable_noAED * p_survive_nonshockable_noAED * p_CPC2_nonshockable/p_CPC_nonshockable,         # non shockable, survive, CPC2
                        
                        "CPC3" = p_shockable_noAED * p_survive_shockable_noAED    * p_CPC3_shockable_noAED/p_CPC_shockable_noAED +  # shockable, survive, CPC3
                              p_nonshockable_noAED * p_survive_nonshockable_noAED * p_CPC3_nonshockable/p_CPC_nonshockable,         # non shockable, survive, CPC3  
                        
                        "CPC4" = p_shockable_noAED * p_survive_shockable_noAED    * p_CPC4_shockable_noAED/p_CPC_shockable_noAED +  # shockable, survive, CPC4
                              p_nonshockable_noAED * p_survive_nonshockable_noAED * p_CPC4_nonshockable/p_CPC_nonshockable,         # non shockable, survive, CPC4
                        
                        "D"   =  p_shockable_noAED * p_arrest_D_shockable_noAED +                                                   # shockable, survive, CPC1
                              p_nonshockable_noAED * p_arrest_D_nonshockable_noAED)                                                 # non shockable, survive, CPC1
    
    
    
    d_e <- d_c # discount rate per cycle equal discount of costs and QALYs by 3%
    
    # Discount weights for costs and effects
    v_dwc <- 1 / (1 + d_c) ^ (c(0:(n_cycles))) 
    v_dwe <- 1 / (1 + d_e) ^ (c(0:(n_cycles)))  # c(0:(n_cycles - 1)) # Here we have a cycle of "1 day"/"our decision tree" in the cycle 0, therefore, we don't discount. The next cycle after, is the "first year" of our model, meaning we don't discount and after we apply discounting as usual. In other words, our first two cycles (cycle 0 + cycle 1) are representing 1 year as they are the sum of the instantaneous event + a subsequent year after the event (together assumed to be the first "year"). 
  
    
    # Within-cycle correction (WCC) using half-cycle correction
    v_wcc <- darthtools::gen_wcc(n_cycles = n_cycles, 
                                 method = "half-cycle") 
   
    #Calculate AED cost per cardiac arrest
    #Number of cardiac arrests per year
    c_AED_Arrest <- (c_AED + c_AED_training*size_household)/(1-(1 - p_cardiac_arrest)^ size_household)
  
    # create the markov trace matrix M capturing the proportion of the cohort in each state 
    # at each cycle
    m_M_noAED <- m_M_AED <- matrix(NA, 
                                   nrow     = n_cycles + 1, ncol = n_states,
                                   dimnames = list(paste("cycle", 0:n_cycles, sep = " "), v_names_states))
    
    head(m_M_noAED) # show first 6 rows of the matrix 
    
    # The cohort starts after the arrest in the first CPC health state
    m_M_AED[1, ]  <- v_m_init_AED # initiate first cycle of cohort trace 
    m_M_noAED[1, ] <- v_m_init_noAED # initiate first cycle of cohort trace  
    
    ## Initialize transition array which will capture transitions from each state to another over time
    
    # for AED
    a_A_AED <- array(0,
                     dim      = c(n_states, n_states, n_cycles + 1),
                     dimnames = list(v_names_states, v_names_states, 0:n_cycles))
    
    
    # Initialize transition array for strategies no AED
    a_A_noAED <- a_A_AED
    
    # Set first slice of a_A with the initial state vector in its diagonal
    diag(a_A_AED[, , 1])   <- v_m_init_AED
    diag(a_A_noAED[, , 1]) <- v_m_init_noAED
    
    
    # create the transition probability array for AED strategy
    # All transitions to a non-death state are assumed to be conditional on survival
    a_P_AED <- array(0,  # Create 3-D array
                     dim = c(n_states, n_states, n_cycles),
                     dimnames = list(v_names_states, v_names_states, 
                                     v_names_cycles[-length(v_names_cycles)])) # name the dimensions of the array
    
    #From CPC1
    a_P_AED["CPC1", "CPC1", ] <- 1 - v_CPC1_D
    a_P_AED["CPC1", "D", ]    <-     v_CPC1_D
    
    #From CPC2
    a_P_AED["CPC2", "CPC2", ] <- 1 - v_CPC2_D
    a_P_AED["CPC2", "D", ]    <-     v_CPC2_D
    
    #From CPC3
    a_P_AED["CPC3", "CPC3", ] <- 1 - v_CPC3_D
    a_P_AED["CPC3", "D", ]    <-     v_CPC3_D
    
    #From CPC4
    a_P_AED["CPC4", "CPC4", ] <- 1 - v_CPC4_D
    a_P_AED["CPC4", "D", ]    <-     v_CPC4_D
    
    #From Dead
    a_P_AED["D", "D", ] <-1
  
    # Store the probability matrix of AED im the noAED one. 
    a_P_noAED <- a_P_AED
    
    ## 05.1 Run Markov model
    for (t in 1:n_cycles){     # loop through the number of cycles
      ## Fill in cohort trace for the next cycle (t + 1)
      
      # For the AED strategy
      m_M_AED[t + 1,   ] <- m_M_AED[t, ]  %*% a_P_AED[, , t]    
      # For the no AED strategy
      m_M_noAED[t + 1, ] <- m_M_noAED[t, ] %*% a_P_noAED[, , t]    
      
      # For the AED strategy
      a_A_AED[, , t + 1]    <- m_M_AED[t, ]   *  a_P_AED[, , t]
      # For the no AED strategy
      a_A_noAED[, , t + 1]  <- m_M_noAED[t, ] *  a_P_noAED[, , t]  
      
      
      ## Store the cohort traces in a list
      l_m_M        <- list(m_M_AED,
                           m_M_noAED)
      names(l_m_M) <- v_names_str
      
      
      ## Store the cohort traces in a list
      l_a_A        <- list(a_A_AED,
                           a_A_noAED)
      names(l_a_A) <- v_names_str
      
    } # close the loop
    
  
    # create the cost array for AED strategy
    a_C_AED <- array(0,  # Create 3-D array
                     dim = c(n_states, n_states, n_cycles),
                     dimnames = list(v_names_states, v_names_states, 
                                     v_names_cycles[-length(v_names_cycles)])) # name the dimensions of the array 
    
    a_C_noAED <- a_C_AED
    
    #AED Costs
    # Common base cost
    c_base <- c_AED_Arrest + c_emergency_department
    
    # AED Costs per CPC level (Year 1)
    c_AED_CPC1_Y1 <- c_base + (c_hospital_CPC1 + c_medical_CPC1) * c_multiplier
    c_AED_CPC2_Y1 <- c_base + (c_hospital_CPC2 + c_medical_CPC2) * c_multiplier
    c_AED_CPC3_Y1 <- c_base + (c_hospital_CPC3 + c_medical_CPC3) * c_multiplier
    c_AED_CPC4_Y1 <- c_base + (c_hospital_CPC4 + c_medical_CPC4) * c_multiplier
    
    c_AED_D_Y1   <-  c_AED_Arrest + (p_AED_used * p_shockable_AED * p_arrest_D_shockable_AED * (c_emergency_department * c_multiplier * (p_die_emergency_department_shockable_AED + p_die_hospital_shockable_AED) + c_hospital_D * c_multiplier * p_die_hospital_shockable_AED)) + (p_AED_used * p_nonshockable_AED * p_arrest_D_nonshockable_AED * (c_emergency_department * c_multiplier * (p_die_emergency_department_nonshockable + p_die_hospital_nonshockable) + c_hospital_D * c_multiplier * p_die_hospital_nonshockable)) + (p_AED_not_used * p_shockable_noAED * p_arrest_D_shockable_noAED * (c_emergency_department * c_multiplier * (p_die_emergency_department_shockable_noAED + p_die_hospital_shockable_noAED) + c_hospital_D * c_multiplier * p_die_hospital_shockable_noAED)) + (p_AED_not_used * p_nonshockable_noAED * p_arrest_D_nonshockable_noAED * (c_emergency_department * c_multiplier * (p_die_emergency_department_nonshockable + p_die_hospital_nonshockable) + c_hospital_D * c_multiplier * p_die_hospital_nonshockable))
    
    #No AED Costs
    c_noAED_CPC1_Y1 <- (c_emergency_department + c_hospital_CPC1 + c_medical_CPC1) * c_multiplier
    c_noAED_CPC2_Y1 <- (c_emergency_department + c_hospital_CPC2 + c_medical_CPC2) * c_multiplier
    c_noAED_CPC3_Y1 <- (c_emergency_department + c_hospital_CPC3 + c_medical_CPC3) * c_multiplier 
    c_noAED_CPC4_Y1 <- (c_emergency_department + c_hospital_CPC4 + c_medical_CPC4) * c_multiplier 
    
    c_noAED_D_Y1    <- (p_shockable_noAED    * p_arrest_D_shockable_noAED  * 
                          (c_emergency_department * c_multiplier * (p_die_emergency_department_shockable_noAED + p_die_hospital_shockable_noAED) + c_hospital_D * c_multiplier * p_die_hospital_shockable_noAED)
    ) + 
      (p_nonshockable_noAED * p_arrest_D_nonshockable_noAED * 
         (c_emergency_department * c_multiplier * (p_die_emergency_department_nonshockable    + p_die_hospital_nonshockable) + c_hospital_D * c_multiplier * p_die_hospital_nonshockable)
      )
    
    #Costs for subsequent years
    c_annual <- c_medical_yearly * c_multiplier
    
    #Create vectors of costs over time 
    #AED
    v_C_AED_CPC1 <- c(c_AED_CPC1_Y1, rep(c_annual, each = n_cycles - 1))
    v_C_AED_CPC2 <- c(c_AED_CPC2_Y1, rep(c_annual, each = n_cycles - 1))
    v_C_AED_CPC3 <- c(c_AED_CPC3_Y1, rep(c_annual, each = n_cycles - 1))
    v_C_AED_CPC4 <- c(c_AED_CPC4_Y1, rep(c_annual, each = n_cycles - 1))
    v_C_AED_D    <- c(c_AED_D_Y1,    rep(0,        each = n_cycles - 1))
    
    #No AED
    v_C_noAED_CPC1 <- c(c_noAED_CPC1_Y1, rep(c_annual, each = n_cycles - 1))
    v_C_noAED_CPC2 <- c(c_noAED_CPC2_Y1, rep(c_annual, each = n_cycles - 1))
    v_C_noAED_CPC3 <- c(c_noAED_CPC3_Y1, rep(c_annual, each = n_cycles - 1))
    v_C_noAED_CPC4 <- c(c_noAED_CPC4_Y1, rep(c_annual, each = n_cycles - 1))
    v_C_noAED_D    <- c(c_noAED_D_Y1,    rep(0,        each = n_cycles - 1))
    
    
    #Costs applied to the AED array
    a_C_AED["CPC1", "CPC1", ] <- v_C_AED_CPC1
    a_C_AED["CPC2", "CPC2", ] <- v_C_AED_CPC2
    a_C_AED["CPC3", "CPC3", ] <- v_C_AED_CPC3
    a_C_AED["CPC4", "CPC4", ] <- v_C_AED_CPC4
    a_C_AED["D",     "D", ]   <- v_C_AED_D
    
    #Costs applied to the no AED array
    a_C_noAED["CPC1", "CPC1", ] <- v_C_noAED_CPC1
    a_C_noAED["CPC2", "CPC2", ] <- v_C_noAED_CPC2
    a_C_noAED["CPC3", "CPC3", ] <- v_C_noAED_CPC3
    a_C_noAED["CPC4", "CPC4", ] <- v_C_noAED_CPC4
    a_C_noAED["D",     "D", ]   <- v_C_noAED_D
    
    # Vectors with utilities by treatment
    #AED
    v_u_AED  <- c(1 - du_CPC1 + u_change, 1 - du_CPC2 + u_change, 1 - du_CPC3 + u_change, 1 - du_CPC4 + u_change, u_D)
    
    #No AED
    v_u_noAED <- v_u_AED
    
    names(v_u_AED) <- names(v_u_noAED) <- v_names_states
    
    l_c <- list(AED   = a_C_AED, 
                noAED = a_C_noAED)
    
    l_u <- list(AED   = v_u_AED,
                noAED = v_u_noAED)
  
    # Calculate QALYs 
    #Check disutility arrest
    #AED
    v_qaly_AED      <- m_M_AED %*% v_u_AED # calculate the QALYs
    v_qaly_AED[1, ] <- v_qaly_AED[2, ] - sum(m_M_AED[2, - n_states] * du_arrest) # calculate the proportion in each CPC state (all without dead), multiply with the disutility for arrest and subtract that in that cycle
    v_tot_qaly_AED <- t(v_qaly_AED) %*% (v_dwe * v_wcc) 
    
    # noAED
    v_qaly_noAED      <- m_M_noAED %*% v_u_noAED
    v_qaly_noAED[1, ] <- v_qaly_noAED[2, ] - sum(m_M_noAED[2, - n_states] * du_arrest) # calculate the proportion in each CPC 
    v_tot_qaly_noAED  <- t(v_qaly_noAED) %*% (v_dwe * v_wcc) 
    
    # Calculate Costs
    # Create matrix to calculate costs over time
    m_R_AED <- m_R_noAED<-  matrix(0,
                                   nrow = n_cycles+1,
                                   ncol = n_states,
                                   dimnames = list( 0:(n_cycles), v_names_states)) # name the columns and rows of the matrix
    
    for (t in 1:(n_cycles)){
      m_R_AED[t+1, ]   <- colSums(a_A_AED[, , t] * a_C_AED[, , t])
      m_R_noAED[t+1, ] <- colSums(a_A_noAED[, , t] * a_C_noAED[, , t])
      
    }
    
    #AED
    v_costs_AED <- rowSums(m_R_AED)
    v_tot_costs_AED <- t(v_costs_AED) %*% (v_dwc * v_wcc) 
    
    # noAED
    v_costs_noAED <- rowSums(m_R_noAED)
    v_tot_costs_noAED <- t(v_costs_noAED) %*% (v_dwc * v_wcc) 
    
    v_tc <- round(c(v_tot_costs_AED, v_tot_costs_noAED))
    v_tu <- round(c(v_tot_qaly_AED,  v_tot_qaly_noAED),3)
  
  
  # Dataframe with discounted costs and effectiveness in QALYs
  df_cu       <- data.frame(Strategy = v_names_str,
                            Cost     = v_tc,
                            Effect   = v_tu) # effects are in utility (QALYs)
  
  # combine results of utility and effects in one dataframe 
  df_results   <- data.frame(Strategy = v_names_str,
                             Cost   = v_tc,
                             QALYs  = v_tu) 
  
  ## Vector with discounted net monetary benefits (NMB) - QALY
  v_nmb    <- (v_tu * n_wtp) - v_tc
  
  ## Vector with discounted net health benefits (NHB) - QALY
  v_nhb    <- v_tu - (v_tc/n_wtp)
  
  df_ce_combined_icer <- calculate_icers(cost       = v_tc,
                                    effect     = v_tu,
                                    strategies = v_names_str)
  
  df_ce_combined <- cbind(df_ce_combined_icer,
                          NMB = v_nmb,
                          NHB = v_nmb)
  
  
  
  
  # Create a summary data frame  
  l_out_all <- list(df_cu          = df_cu,
                    df_results     = df_results,
                    df_ce_combined = df_ce_combined)
  
  # dataframe with combined CE results
  l_out_ce <- df_ce_combined
  
  # if the Return all is true, report the full information
  if(return_all == TRUE){
    l_out_ce <-  l_out_all
  }
  
  # RETURN OUTPUT #
  return(l_out_ce) 
  
    })
}
