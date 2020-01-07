library(assertthat)
library(rlist)
library(data.table)
library(testthat)
library(devtools)

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

ga_fitness_score <- function(list_of_gene, data, fitness = AIC, func = lm, response, min = FALSE) {
  #The new fitness function.
  #list_of_gene is a list of gene likes
  #list(c(TRUE,FALSE,FALSE), c(TRUE,FALSE,TRUE), c(TRUE,FALSE,TRUE)).
  #data is a dataframe containing several x columus.
  #fitness is the fitness function, default in AIC.
  #func is the regression method, likes lm or glm, default in lm.
  #response is the response value (y).
  #min: TRUE is for those fitness function the smaller the fitness value the better the model
  #like AIC,while FALSE is for those the larger fitness value the better the model.

  assert_that(is.data.frame(data), msg = "the data must be a dataframe.")
  assert_that(length(response) == nrow(data),
              msg = "the dimension of observed vectors should be the same as that of response vector.")
  fitness_value <- vector()

  for (i in 1:length(list_of_gene)){
    gene <- list_of_gene[[i]]

    regression_data <- data[,gene]
    regression_data <- data.frame(response,regression_data)
    model <- func(response~., data = regression_data)
    fitness_value <- c(fitness_value,fitness(model))
  }
  if (min == TRUE){
    return (-fitness_value)
  }
  else{
    return (fitness_value)
  }
}

ga_initialize = function(dim, p = 20){
  # This function takes dimension and the number of individuals
  # and returns a list of initialized individuals

  individuals = data.frame(matrix(runif(n = dim*p) > 0.5, dim, p))

  return(as.list(individuals))
}




#' GENETIC ALGORITHMS
#'
#' @param dim A number, indicating the dimension of right-hand-side variables.
#' @param p A number, indicating the number of individuals in your population.
#' @param t A number, indicating the time of iterations.
#' @param selection_method A string, either 'rank' or 'score'. For 'rank', the algorithm will select parents according to their fitness ranks; for 'score', the algorithm will select parents according to their fitness raw scores
#' @param partial_update A logical variable, deciding whether to retain some parents of previous generation when updating population
#' @param  parent_ratio A number between 0 and 1, deciding the ratio of parents retained from previous generation, if selection_method is TRUE
#' @param m_prob A number between 0 and 1, indicating the probability of mutation
#' @param ... A set of parameters to decide how to calculate fitness function, explained below
#' @param data A data.frame of right-hand-variables, which should be (n by dim)
#' @param response A vector of response in your regression, which should be n by 1
#' @param fitness A function, to calculalte the fitness, default AIC
#' @param func A function of linear regression, either 'glm' or 'lm', default lm
#' @param min A logical variable, if TRUE, it will return the negative value of fitness scores
#' @return A list of three :a list of latest generation, a vector of highest fitness score of each iteration, the best individual appeared(result & fitness score)
#' @examples
#'
#' \dontrun{
#' data = read.table('madelon_train.data')
#' y = unlist(read.table('madelon_train.labels'))
#' pop = ga_compute(dim = ncol(data), p = nrow(data), t = 100, selection_method = 'score', partial_update = FALSE, data = data, fitness = AIC, func = glm, response = y, min = TRUE)
#' }
#'
#' \dontrun{
#' set.seed(1)
#' n = 500
#' x = matrix(rnorm(n*3), n, 3)
#' linearcombos = replicate(7, gen_linearcombo(x))
#' probes = matrix(rnorm(n*10, 0, 5), n, 10)
#' betas = rnorm(3)
#' y = rnorm(n, betas[1]*x + betas[2]*x + betas[3]*x)
#' x = cbind(x,linearcombos,probes)
#' x = as.data.frame(x)
#' pop1 = ga_compute(dim = ncol(x), p = nrow(x), t = 100, selection_method = 'score', partial_update = FALSE, data = x, fitness = AIC, func = lm, response = y, min = TRUE)
#' pop2 = ga_compute(dim = ncol(x), p = nrow(x), t = 300, selection_method = 'score', partial_update = FALSE, data = x, fitness = AIC, func = lm, response = y, min = TRUE)
#' }
#'


select = function(dim, p, t = 100, selection_method = 'rank', partial_update = FALSE, parent_ratio = 0.5, m_prob = 0.03, ...){
  # This function computes the GA results
  # dim is the dimension of genes
  # p is the number of individuals in population
  # t is the time of iterating

  assert_that(parent_ratio >= 0 & parent_ratio <= 1, msg = 'Ratio of parents should be between 0 and 1')
  assert_that(nrow(data)>=dim, msg = 'The dimenstion exceeds the length of observed data vector.')
  pop = ga_initialize(dim, p)
  highest_fitness = numeric()
  best_individual = 0
  best_score = -1e8
  for(i in 1:t){
    # Find fitness
    ga_fitness_scores = ga_fitness_score(pop, ...)
    max_fit_score = max(ga_fitness_scores)

    # Find the best individual
    highest_fitness = c(highest_fitness, max_fit_score)
    if(max_fit_score > best_score){
      best_individual = pop[ga_fitness_scores == max_fit_score][1]
      best_score = max_fit_score
    }

    #UPDATE(pop)
    parents = pop[ga_select_index(ga_fitness_scores, method = selection_method)]
    p = length(parents)
    n_cross = floor(p/2)
    children = list(n_cross*2)

    for(i in 1:n_cross){
      childs = ga_crossover(parentA = parents[[i]], parentB = parents[[p-i+1]])
      children[[(i*2-1)]] = childs[[1]]
      children[[(i*2)]] = childs[[2]]
    }

    if(partial_update){
      ind = sample(1:p, floor(parent_ratio*p), replace = FALSE)
      children[ind] = pop[ind]
    }

    children = lapply(children, ga_mutate, mprob = m_prob)
    pop = children
  }

  # Latest Update of best_individual
  ga_fitness_scores = ga_fitness_score(pop, ...)
  if(max_fit_score > best_score){
    best_individual = pop[ga_fitness_scores == max_fit_score][1]
    best_score = max_fit_score
  }
  return(list(pop,highest_fitness, list(best_individual, best_score)))
}

ga_select_index = function(ga_fitness_scores, method = 'rank'){
  assert_that(method %in% c('score', 'rank'), msg = "method should be {score, rank}")
  p = length(ga_fitness_scores)
  if(method == 'rank'){
    fitness_rank = frankv(ga_fitness_scores, order = 1)
    selection_prob = 2*fitness_rank / (p^2 + p)
    parents = sample(1:p, size = p, replace = TRUE, prob = selection_prob)
  }
  if(method == 'score'){
    parents = (1:p)[ga_fitness_scores >= median(ga_fitness_scores)]
    parents = c(parents, sample(1:p, size = p - length(parents), replace = TRUE))
  }
  return(parents)
}

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
