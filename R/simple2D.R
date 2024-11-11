
iterate <- function(func, n, d) {
  index = 1
  if (d == 8) {
    for (i0 in 0:n) { 
      for (i1 in 0:(n - i0)) {
        for (i2 in 0:(n - i0 - i1)) {
          for (i3 in 0:(n - i0 - i1 - i2)) {
            for (i4 in 0:(n - i0 - i1 - i2 - i3)) {
              for (i5 in 0:(n - i0 - i1 - i2 - i3 - i4)) {
                for (i6 in 0:(n - i0 - i1 - i2 - i3 - i4 - i5)) {
                  x = c(i0, i1, i2, i3, i4, i5, i6)
                  x = c(x, n - sum(x))
                  func(x, index)
                  index = index + 1
                }
              }
            }
          }
        }
      }
    }
  } else if (d == 6) {
    for (i0 in 0:n) {
      for (i1 in 0:(n - i0)) {
        for (i2 in 0:(n - i0 - i1)) {
          for (i3 in 0:(n - i0 - i1 - i2)) {
            for (i4 in 0:(n - i0 - i1 - i2 - i3)) {
              x = c(i0, i1, i2, i3, i4)
              x = c(x, n - sum(x))
              func(x, index)
              index = index + 1
            }
          }
        }
      }
    }
  } else if (d == 4) {
    for (i0 in 0:n) {
      for (i1 in 0:(n - i0)) {
        for (i2 in 0:(n - i0 - i1)) {
          func(c(i0, i1, i2, n - i0 - i1 - i2), index)
          index = index + 1
        }
      }
    }
  } else if (d == 3) {
    for (i0 in 0:n) {
      for (i1 in 0:(n - i0)) {
        func(c(i0, i1, n - i0 - i1), index)
        index = index + 1
      }
    }
  } else if (d == 2) {
    for (i0 in 0:n) {
      func(c(i0, n - i0), index)
      index = index + 1
    }
  } else {
    stop("unsupported dimensionality")
  }
}

log <- function(v) {
  result = base::log(v)
  result[is.infinite(result)] = 0 
  return(result)
}

getIndex <- function(x) {
  n = sum(x)
  k = length(x)
  s = 1
  remaining = 0
  for (i in 1:(k - 1)) {
    s = s + x[i]
    remaining = remaining + choose(n + k - i - s, k - i)
  }
  choose(n + k - 1, k - 1) - remaining
}

getTheta <- function(t, m, min_allowed) {
  d = length(t)
  t / m * (1 - d * min_allowed) + min_allowed
}

cardinality <- function(n, d) {
  choose(n + d - 1, d - 1)
}

multiChoose <- function(x) {
  factorial(sum(x)) / prod(factorial(x))
}

getLikelihoodMatrix <- function(n, m, d, min_allowed) {
  pi = matrix(NA, cardinality(n, d), cardinality(m, d))
  iterate(function(x, i) {
    c_x = multiChoose(x)
    iterate(function(t, j) {
      theta = getTheta(t, m, min_allowed)
      pi[i, j] <<- c_x * prod(theta ^ x) 
    }, m, d)
  }, n, d)
  return(pi)
}

getPhiVector <- function(likelihood_matrix) {
  lp = log(rowSums(likelihood_matrix))
  exp(apply(likelihood_matrix, MARGIN = 2, FUN = function(col) sum((log(col) - lp) * col)))
}

getPhiPrecompute <- function(n, m, dim, row_transform, col_transform) {
  cat("likelihood matrix size =", cardinality(n, dim) * cardinality(m, dim) * 8 / 2^20, "MB\n")
  cat("coefficient matrix size =", cardinality(m, dim / 2) * cardinality(m, dim / 2) * 8 / 2^20, "MB\n")
  
  result = list(
    coefficient = matrix(0, cardinality(m, dim / 2), cardinality(m, dim / 2)),
    conditional = c(NA, cardinality(m, dim))
  )
  
  cat("calculating phi vectors...")
  full_phi = getPhiVector(getLikelihoodMatrix(n, m, dim, min_allowed = 0))
  half_phi = getPhiVector(getLikelihoodMatrix(n, m, dim / 2, min_allowed = 0))
  cat("done!\n")
  
  iterate(function(theta, theta_index) {
    i = getIndex(row_transform(theta))
    j = getIndex(col_transform(theta))
    result$coefficient[i, j] <<- result$coefficient[i, j] + full_phi[theta_index]
  }, m, dim)
  
  iterate(function(theta, theta_index) {
    i = getIndex(row_transform(theta))
    j = getIndex(col_transform(theta))
    result$conditional[theta_index] <<- full_phi[theta_index] / result$coefficient[i, j]
  }, m, dim)
  
  for (i in 1:nrow(result$coefficient)) {
    for (j in 1:ncol(result$coefficient)) {
      result$coefficient[i, j] = (result$coefficient[i, j] * half_phi[i] * half_phi[j]) ^ (1 / 3)
    }
  }
  
  return(result)
}

getThetaPrior <- function(m, dim, joint, conditional, row_transform, col_transform) {
  result = c(NA, cardinality(m, dim))
  iterate(function(theta, theta_index) {
    i = getIndex(row_transform(theta))
    j = getIndex(col_transform(theta))
    result[theta_index] <<- conditional[theta_index] * joint[i, j]
  }, m, dim)
  return(result)
}

getProbability <- function(prior, query, m) {
  d = length(query)
  result = 0
  iterate(function(t, i) {
    theta = getTheta(t, m, 0)
    result <<- result + prod(theta ^ query) * prior[i]
  }, m, d)
  return(result)
}

lgTransform <- function(x) {
  if (length(x) == 8) result = c(x[1] + x[2], x[3] + x[4], x[5] + x[6], x[7] + x[8])
  #if (length(x) == 8) result = c(x[1] + x[2] + x[7] + x[8], x[3] + x[4] + x[5] + x[6])
  else if (length(x) == 6) result = c(x[1] + x[2], x[3] + x[4], x[5] + x[6])
  else if (length(x) == 4) result = c(x[1] + x[2], x[3] + x[4])
  else stop("unsupported")
  
  return(result)
}

oeTransform <- function(x) {
  #if (length(x) == 8) result = c(x[1] + x[3] + x[5] + x[7], x[2] + x[4] + x[6] + x[8])
  if (length(x) == 8) result = c(x[1] + x[3], x[2] + x[4], x[5] + x[7], x[6] + x[8])
  #if (length(x) == 8) result = c(x[1] + x[3] + x[6] + x[8], x[2] + x[4] + x[5] + x[7])
  else if (length(x) == 6) result = c(x[1] + x[3], x[5] + x[2], x[4] + x[6])
  else if (length(x) == 4) result = c(x[1] + x[3], x[2] + x[4])
  else stop("unsupported")
  
  return(result)
}

dim = 8
n = 9
m = 7

phi = getPhiPrecompute(n, m, dim, lgTransform, oeTransform)


joint = matrix(runif(cardinality(m, dim / 2)^2), cardinality(m, dim / 2), cardinality(m, dim / 2))
for (iteration in 1:25) {
  row_marginal = rowSums(joint)
  col_marginal = colSums(joint)
  for (i in 1:nrow(joint)) {
    for (j in 1:ncol(joint)) {
      if (joint[i, j] != 0) {
        joint[i, j] = phi$coefficient[i, j] *
          (joint[i, j] / row_marginal[i] * joint[i, j] / col_marginal[j]) ^ (1 / 3)
      }
    }
  }
  
  measure = log(sum(joint)) * 3
  print(c(iteration, measure))
}

prior = getThetaPrior(m, dim, joint, phi$conditional, lgTransform, oeTransform)

# 0.5
k = 1
e = c(1,0,0,0,  0,0,0,1)
w = getProbability(prior, k * e + c(0,0,1,0,  0,0,0,0), m)
l = getProbability(prior, k * e + c(0,0,0,0,  0,0,1,0), m)
w / (w + l)


# These three probabilities should be equal:
k = 1
e = c(0,0,0,1,  1,0,1,0)
w = getProbability(prior, k * e + c(0,1,0,0,  0,0,0,0), m)
l = getProbability(prior, k * e + c(0,0,0,0,  0,1,0,0), m)
w / (w + l)

e = c(0,0,1,0,  1,1,0,0)
w = getProbability(prior, k * e + c(0,0,0,1,  0,0,0,0), m)
l = getProbability(prior, k * e + c(0,0,0,0,  0,0,0,1), m)
w / (w + l)

e = c(1,1,0,0,  0,0,0,1)
w = getProbability(prior, k * e + c(0,0,0,0,  0,0,1,0), m)
l = getProbability(prior, k * e + c(0,0,1,0,  0,0,0,0), m)
w / (w + l)