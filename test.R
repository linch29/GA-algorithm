library(assertthat)
library(testthat)
library(prodlim)

# Test main function select and ga_fitness_score
x = 1:100
data = cbind(x, cos(1:100), sin(1:100))
y = x + rnorm(n = 100)

gene = list()
gene[[1]] = c(TRUE, TRUE, TRUE)
gene[[2]] = c(TRUE, TRUE, FALSE)
gene[[3]] = c(TRUE, FALSE, TRUE)
gene[[4]] = c(FALSE, TRUE, TRUE)
gene[[5]] = c(FALSE, FALSE, TRUE)
gene[[6]] = c(FALSE, TRUE, FALSE)
gene[[7]] = c(TRUE, FALSE, FALSE)

fitness = numeric()
for(i in 1:7){
  mod = lm(y~data[,gene[[i]]])
  fitness = c(fitness, AIC(mod))
}

data = as.data.frame(data)
fitness_scores = ga_fitness_score(list_of_gene = gene, data = data, fitness = AIC, func = lm, response = y, min = FALSE)
test_that("ga_fitness_score works", {
  expect_identical(fitness,fitness_scores)
})

pop_result = select(dim = 3, p = 25, t = 50, m_prob = 0.01, data = data, fitness = AIC, func = lm, response = y, min = TRUE)
names(pop_result[[3]][[1]]) = NULL
test_that("select works",{
  expect_identical(unlist(pop_result[[3]][[1]]), c(TRUE, FALSE, FALSE))
})

# Test fitness_score again
y = c(1,3,5,7,9)
data2 = cbind(c(10,9,5,7,6),  c(7,6,5,4,3), c(1,2,3,4,5))
data3 = data.frame(x1 = c(10,9,5,7,6), x2 = c(7,6,5,4,3), x3 = c(1,2,3,4,5))
gene = list()
gene[[1]] = c(TRUE, TRUE, TRUE)
gene[[2]] = c(TRUE, TRUE, FALSE)
gene[[3]] = c(TRUE, FALSE, TRUE)
gene[[4]] = c(FALSE, TRUE, TRUE)
gene[[5]] = c(FALSE, FALSE, TRUE)
gene[[6]] = c(FALSE, TRUE, FALSE)
gene[[7]] = c(TRUE, FALSE, FALSE)

fitness = numeric()
for(i in 1:7){
  mod = lm(y~data2[,gene[[i]]])
  fitness = c(fitness, AIC(mod))
}

fitness_score = ga_fitness_score(list_of_gene = gene, response = y, data = data3)
test_that('ga_fitness_score works', {
  expect_equal(fitness, fitness_score)
})

# Tests for ga_mutate and ga_crossover
test_that("ga_mutate works",{
  all_combos = unique(expand.grid(c(TRUE,FALSE,FALSE),
                                  c(TRUE,FALSE,FALSE),
                                  c(TRUE,FALSE,FALSE)))
  result = ga_mutate(c(TRUE,FALSE,FALSE))
  check = sum(row.match(all_combos, matrix(result, nrow = 1)), na.rm = T)
  expect_equal(check,1L)
})

test_that("ga_crossover works",{
  parentA = c(TRUE, FALSE, FALSE)
  parentB = c(FALSE, TRUE, TRUE)
  child1 = as.data.frame(matrix(
    c(TRUE,TRUE,TRUE,
      TRUE,FALSE,TRUE,
      TRUE,FALSE,FALSE), byrow = T, nrow = 3))
  child2 = as.data.frame(matrix(
    c(FALSE,FALSE,FALSE,
      FALSE,TRUE,FALSE,
      FALSE,TRUE,TRUE), byrow = T, nrow = 3))
  result = ga_crossover(parentA, parentB)
  check1 = row.match(child1, matrix(result[[1]], nrow = 1))
  check2 = row.match(child2, matrix(result[[2]], nrow = 1))
  expect_identical(check1,check2)
})

