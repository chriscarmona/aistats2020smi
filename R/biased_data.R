
#' @export
SMI_post_biased_data = function( Z, Y, sigma_z, sigma_y, sigma_phi, sigma_theta, sigma_theta_tilde, eta=1 ){
  
  # Calculate de posterior mean and variance of the Semi-Modular Posterior
  # Z, Y : data from the unbiased and biased modules, respectively
  # eta : degree of influence of the unbiased module into phi
  # sigma_z, sigma_y : likelihood std for Z and Y, respectively
  # sigma_phi, sigma_theta : std dev for the priors of phi, theta, theta_tilde respectively
  
  # browser()
  
  n = length(Z)
  m = length(Y)
  
  A = matrix(0,3,3)
  A[1,1] = (n/(sigma_z^2)) + (1+eta) * m/(sigma_y^2) - (m/(sigma_y^2+m*sigma_theta^2)) + 1/(sigma_phi^2)
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
predict_post = function( Z, Y, sigma_z, sigma_y, sigma_phi, sigma_theta, eta=1, stage2=True ) {
  
  # Calculate de Predictive posterior mean and variance of the Semi-Modular Posterior
  
  # Z, Y : data from the unbiased and biased modules, respectively
  # eta : degree of influence of the unbiased module into phi
  # sigma_z, sigma_y : likelihood std for Z and Y, respectively
  # sigma_phi, sigma_theta : prior std for phi and b, respectively
  
  
  n = length(Z)
  m = length(Y)
  
  # Posterior of phi, b, b_tilde
  b = SMI_post_biased_data( Z, Y, sigma_z, sigma_y, sigma_phi, sigma_theta, sigma_theta, eta )
  mean_post = b[[1]]
  cov_post = b[[2]]
  rm(b)
  
  # mean_post, cov_post = phi_b_post( Z, Y, sigma_z, sigma_y, sigma_phi, sigma_theta, eta )
  
  if (stage2){ # keep b
    mean_post = mean_post[c(1,2)]
    cov_post = cov_post[c(1,2),c(1,2)]
  } else{ # keep b tilde
    mean_post = mean_post[c(1,3)]
    cov_post = cov_post[c(1,3),c(1,3)]
  }
  
  cov_post_inv = solve(cov_post)
  
  A = matrix(0,4,4)
  
  A[1,1] = 1/(sigma_z^2)
  A[2,2] = 1/(sigma_y^2)
  A[3,3] = 1/(sigma_z^2)
  A[2:4,2:4] = A[2:4,2:4] + 1/(sigma_y^2)
  A[2:4,2:4] = A[2:4,2:4] + cov_post_inv
  
  A[1,3] = A[3,1] = -1/(sigma_z^2)
  A[2,3] = A[3,2] = -1/(sigma_y^2)
  A[2,4] = A[4,2] = -1/(sigma_y^2)
  
  cov = solve(A)
  b = matrix(0,4,1)
  b[2:4,] = cov_post_inv %*% mean_post
  mean = cov %*% b
  
  return( list( mean[c(0,1),], cov[c(0,1),c(0,1)] ) )
}
