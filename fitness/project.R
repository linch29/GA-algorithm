library(assertthat)
library(testthat)
# The new fitness function.
#list_of_gene is a list of gene likes
#list(c(TRUE,FALSE,FALSE), c(TRUE,FALSE,TRUE), c(TRUE,FALSE,TRUE)).
#data is a dataframe containing several x columus.
#fitness is the fitness function, default in AIC.
#func is the regression method, likes lm or glm, default in lm.
#response is the response value (y).
#min: TRUE is for those fitness function the smaller the fitness value the better the model
#like AIC,while FALSE is for those the larger fitness value the better the model.
ga_fitness_score <- function(list_of_gene, data, fitness = AIC, func = lm, response, min = FALSE) {

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



