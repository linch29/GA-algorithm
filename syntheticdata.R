gen_linearcombo = function(x){
  # This function takes in a dataset and returns 
  # a vector that is a linear combination of 
  # a subset of the columns of the data
  
  keep = rep(FALSE, ncol(x))
  
  # Make sure we get a non trivial linear combination
  while(sum(keep) < 1){
    keep = ga_rlogical(ncol(x))
  }
  
  betas = rnorm(sum(keep))
  
  # In case we get a combination of only one column
  if(sum(keep) == 1){
    return(betas*x[,keep])
  } else{
    return(rowSums(betas*x[,keep]))
  }
}

# Synthetic dataset creation
n = 500
x = matrix(rnorm(n*3), n, 3)
linearcombos = replicate(7, gen_linearcombo(x))
probes = matrix(rnorm(n*10, 0, 5), n, 10)
betas = rnorm(3)
y = rnorm(n, betas[1]*x + betas[2]*x + betas[3]*x)
x = cbind(x,linearcombos,probes)
x = as.data.frame(x)
pop = ga_compute(dim = ncol(x), p = nrow(x), t = 100, selection_method = 'score', partial_update = FALSE, data = x, fitness = AIC, func = lm, response = y, min = TRUE)