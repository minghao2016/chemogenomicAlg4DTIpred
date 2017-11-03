###########################################
###########################################
###########################################
# weights <- c(2, 3)
# kernels <- list(a = replicate(2, rnorm(2)), b = replicate(2, rnorm(2)))
# Map("*", weights, kernels)

## combine kernels
combineKernels <- function(weights, kernels) {
  ## INPUT
  # weight: vector
  # kernels: list
  
  ## OUTPUT:
  
  ## operate in two Lists
  res <- Map("*", weights, kernels)
  
  ## operate in one List
  res <- Reduce("+", res)
  
  return(res)
}
