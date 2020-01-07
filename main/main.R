library(rlist)
library(data.table)
source('fitness/project.R')
source('crossover/crossover.R')




## GA_initialize() function initializes the 
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
#' GA_compute(dim = 500, p = 20, t = 100, selection_method = 'rank', partial_update = TRUE, data = data, fitness = AIC, func = glm, response = y)
 



# This function computes the GA results
# dim is the dimension of genes
# p is the number of individuals in population
# t is the time of iterating
ga_compute = function(dim, p, t = 100, selection_method = 'rank', partial_update = FALSE, parent_ratio = 0.5, m_prob = 0.03, ...){
    assert_that(parent_ratio >= 0 & parent_ratio <= 1, msg = 'Ratio of parents should be between 0 and 1')
    assert_that(nrow(data)>=dim, msg = 'The dimenstion exceeds the length of observed data vector.')
    pop = ga_initialize(dim, p)
    highest_fitness = numeric()
    best_individual = list(0,0)
    for(i in 1:t){
        # Find fitness
        ga_fitness_scores = ga_fitness_score(pop, ...)
        max_fit_score = max(ga_fitness_scores)
        
        # Find the best individual
        highest_fitness = c(highest_fitness, max_fit_score)
        if(max_fit_score > best_individual[[2]]){
            best_individual[[1]] = pop[ga_fitness_scores == max_fit_score]
            best_individual[[2]] = max_fit_score
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
    if(max_fit_score > best_individual[[2]]){
        best_individual[[1]] = pop[ga_fitness_scores == max_fit_score]
        best_individual[[2]] = max_fit_score
    }
    return(list(pop,highest_fitness, best_individual))
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



#y = c(1,3,5,7,9)
#data <- data.frame(x1 = c(10,9,5,7,6), x2 = c(7,6,5,4,3), x3 = c(1,2,3,4,5),
#                   x4 =c (5,4,6,2,4),x5 = c(100,200,300,400,500))
data = read.table('madelon_train.data')
y = unlist(read.table('madelon_train.labels'))
pop = ga_compute(dim = 5, p = 20, t = 100, selection_method = 'score', partial_update = FALSE, data = data, fitness = AIC, func = glm, response = y, min = TRUE)

ft_score = ga_fitness_score(pop[[1]], data, fitness = AIC, func = glm, response = y, min = TRUE)
result = pop[[1]][ft_score == max(ft_score)]
plot(x = 1:100, y = pop[[2]], xlab = 'iteration', ylab = 'highest fitness')
