
#' @export
mu_b_bt_post = function( Z, Y, sigma_z, sigma_y, gamma, rho, rho_tilde, eta=1 ){
  
  # Calculate de posterior mean and variance of the Semi-Modular Posterior
  # Z, Y : data from the unbiased and biased modules, respectively
  # eta : degree of influence of the unbiased module into mu
  # sigma_z, sigma_y : likelihood std for Z and Y, respectively
  # gamma, rho : prior std for mu and b, respectively
  
  # browser()
  
  n = length(Z)
  m = length(Y)
  
  # if eta==0:
  # eta = 0.0000001
  # rho_tilde = rho * np.sqrt(1/eta)
  # rho_tilde = rho * 1
  
  A = matrix(0,3,3)
  A[1,1] = (n/(sigma_z^2)) + (m/(sigma_y^2)) * (1+eta) - (m/(sigma_y^2+m*rho^2)) + 1/(gamma^2)
  A[2,2] = m/(sigma_y^2) + 1/(rho^2)
  A[3,3] = eta * m/(sigma_y^2) + 1/(rho_tilde^2)
  A[1,2] = A[2,1] = m/(sigma_y^2)
  A[1,3] = A[3,1] = eta*m/(sigma_y^2)
  
  cov = solve(A)
  
  aux = matrix(0,3,1)
  aux[1,1] = sum(Z)/sigma_z^2 + (sum(Y)/sigma_y^2) * ( (eta*(sigma_y^2 + m*rho^2) + m*rho^2) / (sigma_y^2 + m*rho^2) )
  aux[2,1] = sum(Y)/sigma_y^2
  aux[3,1] = eta*sum(Y)/sigma_y^2
  mean = cov %*% aux
  
  return( list(mean, cov) )
}

#' @export
mu_b_post = function( Z, Y, sigma_z, sigma_y, gamma, rho, eta=1 ){
  # Calculate de posterior mean and variance of the Semi-Modular Posterior
  # Z, Y : data from the unbiased and biased modules, respectively
  # eta : degree of influence of the unbiased module into mu
  # sigma_z, sigma_y : likelihood std for Z and Y, respectively
  # gamma, rho : prior std for mu and b, respectively
  n = length(Z)
  m = length(Y)
  
  A = matrix(0,2,2)
  A[1,1] = (n/(sigma_z^2)) + (m/(sigma_y^2)) * ( ( eta*sigma_y^2 + m * rho^2) / ( sigma_y^2 + m * rho^2) ) + 1 / gamma^2
  A[1,2] = A[2,1] = m/(sigma_y^2)
  A[2,2] = m/(sigma_y^2) + 1/(rho^2)
  cov = solve(A)
  
  aux = matrix(0,2,1)
  aux[1,1] = sum(Z)/sigma_z^2 + (sum(Y)/sigma_y^2) * ( (eta*sigma_y^2 + m * rho^2) / (sigma_y^2 + m * rho^2) )
  aux[2,1] = sum(Y)/sigma_y^2
  mean = cov %*% aux
  return( list(mean, cov) )
}

#' @export
mu_post = function( Z, Y, sigma_z, sigma_y, gamma, rho, eta=1 ) {
  # Calculate de posterior mean and variance of the Semi-Modular Posterior
  # Z, Y : data from the unbiased and biased modules, respectively
  # eta : degree of influence of the unbiased module into mu
  # sigma_z, sigma_y : likelihood std for Z and Y, respectively
  # gamma, rho : prior std for mu and b, respectively
  
  n = length(Z)
  m = length(Y)
  
  var = n / sigma_z^2 + eta * ( m /( sigma_y^2 + m * rho^2) ) + 1 / gamma^2
  var = 1/var
  
  mean = sum(Z)/sigma_z^2 + eta * sum(Y) / (sigma_y^2 + m * rho^2)
  mean = mean * var
  
  return( list(mean, var) )
  
}

#' @export
predict_post = function( Z, Y, sigma_z, sigma_y, gamma, rho, eta=1, stage2=True ) {
  
  # Calculate de Predictive posterior mean and variance of the Semi-Modular Posterior
  
  # Z, Y : data from the unbiased and biased modules, respectively
  # eta : degree of influence of the unbiased module into mu
  # sigma_z, sigma_y : likelihood std for Z and Y, respectively
  # gamma, rho : prior std for mu and b, respectively
  
  
  n = length(Z)
  m = length(Y)
  
  # Posterior of mu, b, b_tilde
  aux = mu_b_bt_post( Z, Y, sigma_z, sigma_y, gamma, rho, rho, eta )
  mean_post = aux[[1]]
  cov_post = aux[[2]]
  rm(aux)
  
  # mean_post, cov_post = mu_b_post( Z, Y, sigma_z, sigma_y, gamma, rho, eta )
  
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
  aux = matrix(0,4,1)
  aux[2:4,] = cov_post_inv %*% mean_post
  mean = cov %*% aux
  
  return( list( mean[c(0,1),], cov[c(0,1),c(0,1)] ) )
}
