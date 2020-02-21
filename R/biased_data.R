
#' @export
SMI_post_biased_data = function( Z, Y, sigma_z, sigma_y, sigma_phi, sigma_theta, sigma_theta_tilde, eta=1 ){
  
  # Calculates the posterior mean and variance of the Semi-Modular Posterior
  # P_{smi,\eta}( \varphi, \theta, \tilde\theta | Z, Y )
  
  # Z, Y : data from the unbiased and biased modules, respectively
  # eta : degree of influence of the unbiased module into phi
  # sigma_z, sigma_y : likelihood std for Z and Y, respectively
  # sigma_phi, sigma_theta : std dev for the priors of phi, theta, theta_tilde respectively
  
  n = length(Z)
  m = length(Y)
  
  A = matrix(0,3,3)
  A[1,1] = n/(sigma_z^2) + (1+eta) * m/(sigma_y^2) - (m/(sigma_y^2+m*sigma_theta^2)) + 1/(sigma_phi^2)
  A[2,2] = m/(sigma_y^2) + 1/(sigma_theta^2)
  A[3,3] = eta * m/(sigma_y^2) + 1/(sigma_theta_tilde^2)
  A[1,2] = A[2,1] = m/(sigma_y^2)
  A[1,3] = A[3,1] = eta*m/(sigma_y^2)
  
  b = matrix(0,3,1)
  b[1,1] = sum(Z)/(sigma_z^2) + (1+eta) * sum(Y)/(sigma_y^2) - sum(Y) / (sigma_y^2+m*sigma_theta^2)
  b[2,1] = sum(Y)/(sigma_y^2)
  b[3,1] = eta*sum(Y)/(sigma_y^2)
  
  cov = solve(A)
  mean = cov %*% b
  
  return( list(mean, cov) )
}

#' @export
SMI_pred_biased_data = function( Z, Y, sigma_z, sigma_y, sigma_phi, sigma_theta, sigma_theta_tilde, eta=1 ){
  
  # Calculate de joint Predictive posterior mean and variance of the Semi-Modular Posterior
  # P_{smi,\eta}( \varphi, \theta, \tilde\theta, Z_0, Y_0 | Z, Y )
  
  # Z, Y : data from the unbiased and biased modules, respectively
  # eta : degree of influence of the unbiased module into phi
  # sigma_z, sigma_y : likelihood std for Z and Y, respectively
  # sigma_phi, sigma_theta : prior std for phi and b, respectively
  
  n = length(Z)
  m = length(Y)
  
  A = matrix(0,5,5)
  A[1,1] = (n+1)/(sigma_z^2) + (m+eta*m+1)/(sigma_y^2) - (m/(sigma_y^2+m*sigma_theta^2)) + 1/(sigma_phi^2)
  A[2,2] = (m+1)/(sigma_y^2) + 1/(sigma_theta^2)
  A[3,3] = eta * m/(sigma_y^2) + 1/(sigma_theta_tilde^2)
  A[4,4] = 1/(sigma_z^2)
  A[5,5] = 1/(sigma_y^2)
  
  A[1,2] = A[2,1] = (m+1)/(sigma_y^2)
  A[1,3] = A[3,1] = eta*m/(sigma_y^2)
  A[1,4] = A[4,1] = -1/(sigma_z^2)
  A[1,5] = A[5,1] = -1/(sigma_y^2)
  A[2,5] = A[5,2] = -1/(sigma_y^2)
  
  b = matrix(0,5,1)
  b[1,1] = sum(Z)/(sigma_z^2) + (1+eta) * sum(Y)/(sigma_y^2) - sum(Y) / (sigma_y^2+m*sigma_theta^2)
  b[2,1] = sum(Y)/(sigma_y^2)
  b[3,1] = eta*sum(Y)/(sigma_y^2)
  
  cov = solve(A)
  mean = cov %*% b
  
  # return( list(mean, cov) )
  return( list(mean[4:5], cov[4:5,4:5]) )
}
