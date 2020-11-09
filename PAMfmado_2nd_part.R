library(cluster)
PAMfmado_2ndpart <-function(DD, K){
  # Description 
  #   This function performs the PAM algorithm based on the F-madogram distance
  #
  #
  # Aguments 
  #   DD the matrix distance
  #   K number of clusters
  
  
  #---------- CLUSTERING WITH PAM -------#
  output = pam(DD,K,diss = TRUE,medoids = NULL)
  
  return(output)
}
