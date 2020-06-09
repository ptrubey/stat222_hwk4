# functions.R 

hat_gamma = function(pi, alpha){
  # Computes the value for (gamma_l1,gamma_l2)
  # that minimizes KL for each l
  N = ncol(pi) # Truncation level of VB
  
  sumpi = colSums(pi)[-N] # sum(1{L_i = l})
  csump = c(rev(cumsum(rev(sumpi[-1]))),0) # checked.
  
  return(cbind(1 + sumpi, alpha + csump))
}
hat_xi    = function(beta, pi, y, pri){
  s.ypi = as.numeric(t(pi) %*% y)
  s.pi  = colSums(pi)
  
  b1b2i = beta[1] / beta[2]
  
  xi1 = (pri$m + pri$s2 * b1b2i * s.ypi)/(1 + pri$s2 * b1b2i * s.pi)
  xi2 = 1/(1/pri$s2 + b1b2i * s.pi)
  return(cbind(xi1, xi2))
}
hat_pi    = function(beta, xi, gamma, y){
  # This function needs help.
  # should be cleared up....
  N     = nrow(xi)
  b1b2i = beta[1] / beta[2]
  
  C1 = 0.5 * (digamma(beta[1]) - log(beta[2]))
  C2 = 0.5 * b1b2i * sapply(1:N, function(l){
    y^2 - 2 * y * xi[l,1] + xi[l,2] + xi[l,1]^2
    })
  C3 = c(digamma(gamma[,1]) - digamma(gamma[,1] + gamma[,2]), 0)
  C4 = cumsum(c(0, digamma(gamma[,2]) - digamma(gamma[,1] + gamma[,2])))
  
  WW = t(C1 + t(C2) + C3 + C4) # Fixed, should be good.
  W = WW - apply(WW, 1, min)
  pi = exp(W) / apply(exp(W), 1, sum)
  return(pi)
}
hat_beta  = function(xi, pi, y, pri){
  n = length(y)
  N = nrow(xi)
  C = sapply(1:N, function(l){
    pi[,l] * (y^2 - 2 * y * xi[l,1] + xi[l,2] + xi[l,1]^2)
    })
  return(c(pri$a_phi + 0.5 * n, pri$b_phi + 0.5 * sum(C)))
}

convergence = function(be, ga, xi, pi, max.diff = 1e-6){
  diff = (
    + sum((be[2,] - be[1,])^2)
    + sum((ga[2,,] - ga[1,,])^2)
    + sum((xi[2,,] - xi[1,,])^2)
    + sum((pi[2,,] - pi[1,,])^2)
  )
  return(diff < max.diff)
}
convergence = function(be, ga, xi, pi, max.diff = 1e-6){
  diff = (
    + sum((be[2,] - be[1,])^2)
    + sum((ga[2,,] - ga[1,,])^2)
    + sum((xi[2,,] - xi[1,,])^2)
    + sum((pi[2,,] - pi[1,,])^2)
  )
  return(diff)
}

variational = function(y, N, alpha, pri, max.iter = 10000){
  i.pi = array(NA, dim = c(max.iter, length(y), N))
  i.ga = array(NA, dim = c(max.iter, N - 1, 2))
  i.xi = array(NA, dim = c(max.iter, N, 2))
  i.be = array(NA, dim = c(max.iter, 2))
  
  i.ga[1,,] = rgamma((N - 1) * 2, 2, 1)
  i.xi[1,,1] = rnorm(N, pri$mu, sqrt(pri$s2)); i.xi[1,,2] = pri$s2
  i.be[1,]  = rgamma(2, 2, 1)
  i.pi[1,,] = hat_pi(i.be[1,], i.xi[1,,], i.ga[1,,], y)
  
  converged = FALSE
  i = 1
  
  while(!converged && (i < max.iter)){
    i = i + 1
    print(i)
    
    i.ga[i,,] = hat_gamma(i.pi[i-1,,], alpha)
    i.pi[i,,] = hat_pi(i.be[i-1,], i.xi[i-1,,], i.ga[i,,], y)
    i.xi[i,,] = hat_xi(i.be[i-1,], i.pi[i,,], y, pri)
    i.be[i,]  = hat_beta(i.xi[i,,], i.pi[i,,], y, pri)
    
    converged = convergence(i.be[(i-1):i,],  i.ga[(i-1):i,,], 
                            i.xi[(i-1):i,,], i.pi[(i-1):i,,])
  }
  return(list(
    pi = i.pi[iter.cnt,,],
    ga = i.ga[iter.cnt,,],
    xi = i.xi[iter.cnt,,],
    be = i.be[iter.cnt,]
    ))
}











# EOF