getPredictor <- function(x, y, alpha, beta, basis) {
  X = t(sapply(x, basis))
  Y = matrix(y)
  m = ncol(X)
  A = alpha * diag(m) + beta * t(X) %*% X
  M = beta * solve(A) %*% t(X) %*% Y
  
  return(list(
    predict = function(x) basis(x) %*% M, 
    log_marginal = logMarginalLikelihood(X, Y, alpha, beta, A, M),
    M = M)
  )
}

logMarginalLikelihood <- function(X, Y, alpha, beta, A, M) {
  m = ncol(X)
  n = nrow(X)
  
  err_M = beta * sum((Y - X %*% M) ^ 2) + alpha * t(M) %*% M
  
  l = m * log(alpha) + n * log(beta) - err_M - log(det(A)) - n * log(2 * pi)
  return(l / 2)
}

getCentralPredictor <- function(predictors, log_marginals) {
  weights = exp(log_marginals)
  weights = weights / sum(weights)
  
  function(x) sum(weights * sapply(predictors, function(pred) pred(x)))
}

evaluate <- function(predictor, target, test_data) {
  err = sapply(test_data, function(x) predictor(x) - target(x))
  return(sqrt(sum(err ^ 2) / length(test_data)))
}

getPlynomialBasis <- function(n) {
  function(x) x ^ (0:n)
}

###########################################################################

alpha = 0.001
beta = 100

train_size = 7
test_size = 300
min_input = -1
max_input = +1

max_power = 6

x = seq(min_input, max_input, length.out = train_size)  
#input = runif(train_size, min_input, max_input)

target = function(x) sin(pi * x)
noise = rnorm(n = length(x), mean = 0, sd = 1 / beta)
y = target(x) + noise

plot(x, y)

predictors = c()
log_marginals = c()
for (i in 1:max_power) {
  model = getPredictor(x, y, alpha , beta, basis = getPlynomialBasis(i))
  predictors = c(predictors, model$predict)
  log_marginals = c(log_marginals, model$log_marginal)
}

log_marginals
sapply(predictors, function(predict) predict(1.2))
target(1.2)

test = runif(test_size, min_input, max_input)

sapply(predictors, evaluate, target = target, test_data = test)

evaluate(predictor = getCentralPredictor(predictors, log_marginals), target, test)


#############################################################################