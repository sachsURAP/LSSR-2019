#' @title Compares numerical slope behavior between two IDERS
#'  
#' @description \code{dropIDER} Takes two IDER functions and returns the lesser IDER function and the dosage at which the slope of that IDER approaches a numerical bound.
#' (
#' @details dropIDER numerically approximates the slope of each IDER over a vector of two thousand potential dosages in units of gray.
#' 
#' @param ider1 An IDER function such that takes in a dosage in gray for its first argument and LET as its second argument
#' @param ider2 A different IDER function such that takes in a dosage in gray for its first argument and LET as its second argument
#' @param L1 The LET parameter for the first IDER. Defaults to NULL
#' @param L2 The LET parameter for the second IDER. Defaults to NULL
#' @param bound The slope for which dosage will be approximated. Defaults to 0.0001
#' @param max_dose The maximum dosage to analyze the IDERs over. Defaults to 20 gray
#' @param increment The increment between each dosage to be tested. Defaults to 0.01 gray
#' 
# @example 
# HZEC <- function(Dose, L = 173) {
#     return(1-exp(-0.01*(HZEc[1]*L*Dose*exp(-HZEc[2]*L)+(1-exp(-150*phi*Dose/L))*HZEc[3])))
# }
#
# low_LET_IDER <- function(d, lambda, beta) {
#     return(beta*(1-exp(-lambda*d)))
# }
# lower_IDER, lower_dose <- dropIDER(HZEC, low_LET_IDER(lambda = 0.0012, beta = 0.0038), L1 = 173, L2 = NULL)
#'
#' @return The lower IDER function, and numeric dose at the specified bound
#'
#' @author Edward G. Huang <edwardgh@@berkeley.edu>
#' @export

dropIDER <- function(ider1, ider2, L1 = NULL, L2 = NULL, bound = 0.0001, max_dose = 20, increment = 0.01) {
  dose <- seq(0, max_dose, increment)
  ider1 <- ider1(dose, L1); ider2 <- ider2(dose, L2)
  if (mean(ider1) >= mean(ider2)) { # Calculates mean value of each ider.
    lesser <- ider1
    } else {lesser <- ider2 # Designates the lower mean ider to be the lesser ider.
  } 
  slope <- diff(lesser) / diff(dose) # Numerical approximation of derivative of the lower IDER
  dosage_index <- which(slope <= bound)[1] # Returns indice of first element that is less that the approximate bound
  dosage <- dose[dosage_index] # Finds the corresponding dosage
  return(lesser, dosage)
}

HZEC <- function(Dose, L = 173) {
    return(1-exp(-0.01*(HZEc[1]*L*Dose*exp(-HZEc[2]*L)+(1-exp(-150*phi*Dose/L))*HZEc[3])))
}

low_LET_IDER <- function(Dose, lambda, beta) {
    return(beta*(1-exp(-lambda*Dose)))
}

HZEC=function(Dose,L) 1-exp(-0.01*(HZEc[1]*L*Dose*exp(-HZEc[2]*L)+(1-exp(-150*phi*Dose/L))*HZEc[3]))#calibrated IDER. The above equations imply: HZEC(dose=0)=0; HZEC increases monotonically and apporaches 1 at large doses.  

