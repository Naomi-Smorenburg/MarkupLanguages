###############################################
# Author: Naomi Smorenburg
# Utility function maximizing the overall life expectancy for the entire population for a combination of values for t1 and t2.
###############################################

# OLE = overall life expectancy 

OLE <- function(N, hE, hC, t1, t2, t_star) { 
  
  ### Function input ---
  
  # N = patient horizon 
  # hE = hazard rate experimental treatment 
  # hC = hazard rate control treatment 
  # t1 = starting point of trial (in years)
  # t2 = endpoint of trial (in years)
  # t_star = time until N new onset  
  
  ### Computations ---
  
  # 1. Calculate the sample size of the different groups
  n1 <- 0.5*N*(t1/t_star)           # Patients in the trial receiving the standard treatment  
  n2 <- 0.5*N*(t1/t_star)           # Patients in the trial receiving the experimental treatment
  n3 <- N*((t2-t1)/t_star)          # Patients not in the trial, with disease onset between t1 and t2
  n4 <- N*((t_star-t2)/t_star)      # Patients not in the trial, with disease onset after t2. 
  
  # 2. Calculating the average life expectancy when E is not superior 
  C.LE1 <- 1/hC                                                                   # group only receives C 
  C.LE2 <- (1/hE^2) * (hE + (((exp(t1*hE)-1)*exp(-t2*hE)*(hE-hC))/(t1*hC)))       # group receives E first, after t' C
  C.LE3 <- 1/hC                                                                   # group only receives C
  C.LE4 <- 1/hC                                                                   # group only receives C
  
  # 3. Calculating the average life expectancy when E is superior
  E.LE1 <- (1/hC^2) * (hC + (((exp(t1*hC)-1)*exp(-t2*hC)*(hC-hE))/(t1*hE)))       # group receives C first, after t' E                                                       
  E.LE2 <- 1/hE                                                                   # group only receives E
  E.LE3 <- (1/hC^2) * (hC + (((exp(t1*hC)-1)*exp(-t2*hC)*(hC-hE))/(t1*hE)))       # group receives C first, after t' E                                                      
  E.LE4 <- 1/hE                                                                   # group only receives E
  
  # 4. Calculating the overall life expectancy for the entire population N for situation E is superior and E is not          superior
  
  #### By multiplying group size with LE for individuals 
  C.OLE <- (n1 * C.LE1) + (n2 * C.LE2) + (n3 * C.LE3) + (n4 * C.LE4) 
  E.OLE <- (n1 * E.LE1) + (n2 * E.LE2) + (n3 * E.LE3) + (n4 * E.LE4) 
  
  # 5. Calculating the power 
  
  ### First calculate the expected number of events for each group 
  D.E <- n2 * (1 + ((exp(-1*(t2-t1)*hE) - exp(-1*t2*hE))/((t2-t1)*hE - t2*hE)))
  D.C <- n1 * (1 + ((exp(-1*(t2-t1)*hC) - exp(-1*t2*hC))/((t2-t1)*hC - t2*hC)))
  
  ### Calculate the total expected number of events by adding  both groups 
  D <- D.E + D.C
  
  #### PropE = proportion in E group, in the case of 1:1 allocation this is 0.5
  PropE <- 0.5
  
  ### Calculate the mean under the alternative hypothesis 
  ### Log-rank test statistic T has distribution (Schoenfeld, 1981 Biometrika)
  muHa <- log(hE/hC) * sqrt(PropE * (1-PropE) * D)
  
  ### Determine the critical value for one sided alpha of 0.05. 
  Cr <- qnorm(0.025, 0, 1) 
  
  #### Determine the power by taking the area left of the critical value 
  pwr <- pnorm(Cr, muHa, 1)
  
  # 6. Weigh the two possible outcomes by the power 
  U <- (pwr*E.OLE) + ((1-pwr)*C.OLE) 
  
  
  ### Output --- 
  
  # 7. Take the results and put in a list 
  result <- list(
    'Overall life expectancy E not superior' = C.OLE,
    'Overall life expectancy E superior' = E.OLE,
    'Utility' = U
  )
  
  # 6. Return a list with the results 
  return(result)
  
  # End of function 
}