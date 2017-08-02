dropIDER <- function(ider1, ider2, L1 = NULL, L2 = NULL, bound = 0.0001) {
  dose <- seq(0, 20, 0.01)
  ider1 <- ider1(dose, L1); ider2 <- ider2(dose, L2)
  if (mean(ider1) >= mean(ider2)) { # Calculates mean value of each ider.
    lesser <- ider1
  } else { 
    lesser <- ider2
  } # Designates the lower mean ider to be the lesser ider.
  slope <- diff(lesser) / diff(dose) # Numerical approximation of derivative.
  dosage_index <- which(slope > bound + bound * 0.1 & 
                        slope < bound - bound * 0.1)[1] # Returns indice of first element that is in the specified interval. 
  dosage <- dose[dosage_index] # Finds the corresponding dosage.
  return(lesser, dosage)
}

HZEC <- function(Dose, L = 173) {
  return(1-exp(-0.01*(HZEc[1]*L*Dose*exp(-HZEc[2]*L)+(1-exp(-150*phi*Dose/L))*HZEc[3])))
}

low_LET_IDER <- function(lambda, beta, d) {
  return(beta*(1-exp(-lambda*d)))
}


