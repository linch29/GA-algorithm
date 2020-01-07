ga_rlogical = function(size, p = 0.5){
  # This function samples from Trues and Falses with 
  # probability 1/2 each by default
  
  return(sample(c(TRUE, FALSE), size,
                replace = TRUE, prob = c(p, 1-p)))
}

ga_crossover = function(parentA, parentB){
  # This function take in two parents and uses a simple split
  # to create children. It then picks a random child to pass on 
  # in the algorithm
  
  n = length(parentA)
  splitpt = sample(1:(n-1),1)
  
  # Create two children
  childA = parentA
  childA[(splitpt+1):n] = parentB[(splitpt+1):n]
  childB = parentB
  childB[(splitpt+1):n] = parentA[(splitpt+1):n]
  
  # Return both children
  return(list(childA, childB))  
}

ga_mutate = function(child, mprob = 0.03){
  # This function takes in a child and randomly mutates each
  # gene with probability mprob.
  
  # Create vector of positions to mutate
  changes = ga_rlogical(length(child), mprob)

  # For the subset of child that needs to be mutated, flip values
  child[changes] = !child[changes]  
  
  return(child)
}

n = 20

parents = matrix(ga_rlogical(n^2), n, n) 
child_a = ga_crossover(parents[,1], parents[,2])
child_b = ga_mutate(child_a)
