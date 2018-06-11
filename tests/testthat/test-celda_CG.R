# celda_CG.R
# library(celda)
library(testthat)
library(Rtsne)
context("Testing celda_CG")

source('R/celda_functions.R')
source('R/celda_list.R')
source('R/celda.R')
source('R/celda_C.R')
source('R/celda_G.R')
source('R/celda_CG.R')
source('R/diffExp.R')
source('R/semi_pheatmap.R')
source('R/split_clusters.R')
source('R/s3_generics.R')
source('R/StateHeatmap.R')
source('R/plot_dr.R')
source('R/feature_selection.R')
source('R/celda_heatmap.R')

library(doParallel)


#Loading pre-made simulatedcells/celda objects
load("tests/celdaCGsim.rda")
load("tests/celdaCG.rda")
load("tests/celdaCG.non.stochastic.res.rda")
model_CG = getModel(celdaCG.res, K = 5, L = 3)[[1]]
factorized <- factorizeMatrix(model_CG, celdaCG.sim$counts)
counts.matrix <- celdaCG.sim$counts

print("starting CG tests.")

test_that(desc = "celda_CG runs without crashing", {
  # celdaCG.test.res <- celda(counts = celdaCG.sim$counts, model = "celda_CG",
  #                                nchains = 1, K = celdaCG.sim$K, L = celdaCG.sim$L, max.iter = 15)

   celdaCG.test.res <- celda_CG(counts = celdaCG.sim$counts,
                                    K = celdaCG.sim$K, L = celdaCG.sim$L, max.iter = 15,
                                    random.state.order = FALSE)

  # Simple check: do the cell/gene cluster lengths match the provided counts matrix dims?
  #expect_equal(length(celdaCG.test.res$res.list[[1]]$z), ncol(celdaCG.sim$counts))
  #expect_equal(length(celdaCG.test.res$res.list[[1]]$y), nrow(celdaCG.sim$counts))
  expect_equal(length(celdaCG.test.res$z), ncol(celdaCG.sim$counts))
  expect_equal(length(celdaCG.test.res$y), nrow(celdaCG.sim$counts))

})


test_that(desc = "Ensure celda_CG always returns same output in 'non-stochastic' mode", {

  # celdaCG.test.res <- celda(counts = celdaCG.sim$counts, model = "celda_CG",
  #                          nchains = 1, K = celdaCG.sim$K, L = celdaCG.sim$L, max.iter = 15,
  #                          random.state.order = FALSE)

  celdaCG.test.res <- celda_CG(counts = celdaCG.sim$counts,
                            K = celdaCG.sim$K, L = celdaCG.sim$L, max.iter = 15,
                            random.state.order = FALSE)
  
  # Non-stochasticity: make sure cluster labels are consistent with reference celda run

  # test.model <- celdaCG.test.res$res.list[[1]]
  # expect.model <- celdaCG.non.stochastic.res$res.list[[1]]
  test.model <- celdaCG.test.res
  expect.model <- celdaCG.non.stochastic.res$res.list[[1]]

  expect_equal(test.model$z, expect.model$z)
  expect_equal(test.model$y, expect.model$y)
  expect_equal(test.model$finalLogLik, expect.model$finalLogLik)
})



# test_that(desc = "Cell simulation works", {
#   celdacg <- simulateCells(K = 5, L = 3, model = "celda_CG")
#   expect_equal(celdacg$K, 5)
#   expect_equal(celdacg$L, 3)
# })


#Making sure getModel if functioning correctly
test_that(desc = "Sanity checking getModel",{
  expect_equal(celdaCG.res$content.type, class(model_CG))
})

#Making sure relationship of counts vs proportions is correct in factorize matrix
test_that(desc = "Checking factorize matrix, counts vs proportions",{
  expect_equal(TRUE,all(factorized$counts$sample.states/sum(factorized$counts$sample.states) 
                        == factorized$proportions$sample.states))
})

#Checking dimension of factorize matrix
test_that(desc = "Checking factorize matrix dimension size",{
  expect_equal(5, ncol(factorized$proportions$population.states))
  expect_equal(3, nrow(factorized$proportions$population.states))
})


# Ensure calculateLoglikFromVariables calculates the expected values
test_that(desc = "calculateLoglikFromVariables.celda_CG returns correct output for various params", {
  expect_lt(calculateLoglikFromVariables(y = celdaCG.sim$y, z = celdaCG.sim$z, 
                                            delta = 1, gamma = 1,  beta = 1, 
                                            alpha = 1, K = 5, L = 3, model="celda_CG", 
                                            s = celdaCG.sim$sample.label, 
                                            counts=celdaCG.sim$counts),
               0)
})


#normalizeCounts
test_that(desc = "Making sure normalizeCounts doesn't change dimensions",{
  norm.counts <- normalizeCounts(counts.matrix)
  expect_equal(dim(norm.counts),dim(counts.matrix))
  expect_equal(rownames(norm.counts),rownames(counts.matrix))
  expect_equal(colnames(norm.counts),colnames(counts.matrix))
})

#recodeClusterY
test_that(desc = "Checking to see if recodeClusterY gives/doesn't give error",{
  expect_error(recodeClusterY(celda.mod = model_CG, from = NULL, to = ))
  expect_error(recodeClusterY(celda.mod = model_CG, from = c(1,2,3), to = c(1,2,4)))
  new.recoded <- recodeClusterY(celda.mod = model_CG, from = c(1,2,3), to = c(3,2,1))
  expect_equal(model_CG$y == 1,new.recoded$y == 3)
})

#recodeClusterZ
test_that(desc = "Checking to see if recodeClusterZ gives/doesn't give error",{
  expect_error(recodeClusterZ(celda.mod = model_CG, from = NULL, to = ))
  expect_error(recodeClusterZ(celda.mod = model_CG, from = c(1,2,3,4,5), to = c(1,2,3,4,6)))
  new.recoded <- recodeClusterZ(celda.mod = model_CG, from = c(1,2,3,4,5), to = c(5,4,3,2,1))
  expect_equal(model_CG$z == 1,new.recoded$z == 5)
})

#compareCountMatrix
test_that(desc = "Checking CompareCountMatrix",{
  expect_true(compareCountMatrix(count.matrix = celdaCG.sim$counts, celda.obj = model_CG))
})

#distinct_colors
test_that(desc = "Making sure distinct_colors gives expected output",{
  expect_equal(distinct_colors(2), c("#FF9999","#99FFFF"))
  expect_equal(distinct_colors(4), c("#FF9999","#99FFFF","#FFDB99","#9999FF"))
})


##renderCeldaHeatmap###
# test_that(desc = "Checking renderCeldaHeatmap",{
#   expect_equal(names(renderCeldaHeatmap(counts = celdaCG.sim$counts, z = model_CG$z, y = model_CG$y)),
#                c("tree_row","tree_col","kmeans","gtable"))
# })


##feature_selection.R##
#topRank
test_that(desc = "Checking topRank",{

  top.rank <- topRank(fm = factorized$proportions$gene.states, n = 1000)
  expect_equal(nrow(counts.matrix),
               sum(sapply(top.rank$names,FUN = length)))
  expect_equal(names(top.rank),
               c("index","names"))
})

#GiniPlot
test_that(desc = "Checking GiniPlot to see if it runs",{
  gini.plot <- GiniPlot(counts = celdaCG.sim$counts, celda.mod = model_CG)
  expect_equal(class(gini.plot),
               c("gg","ggplot"))
})


# #stateHeatmap
# test_that(desc = "Checking stateHeatmap to see if it runs",{
#   expect_equal(names(stateHeatmap(celdaCG.sim$counts, celda.mod = model_CG)),
#                c("tree_row","tree_col","kmeans","gtable"))
# })


#diffExp
test_that(desc = "Checking diffExp",{
 expect_equal(class(diffexp_K1 <- diffExp(counts = counts.matrix, celda.mod = model_CG, c1 = 1)),
		c("data.table","data.frame"))
})


# #plotDrCluster,State
# test_that(desc = "Checking plotDrCluster to see if it runs",{
#   celda.tsne <- celdaTsne(counts = celdaCG.sim$counts, max.iter = 50,celda.mod = model_CG)
#   expect_equal(names(plotDrCluster(dim1 = celda.tsne[,1], dim2 = celda.tsne[,2],cluster = as.factor(model_CG$z))),
#                c("data","layers","scales","mapping","theme","coordinates","facet","plot_env","labels","guides"))
#   expect_equal(names(plotDrState(dim1 = celda.tsne[,1], dim2 = celda.tsne[,2],matrix = factorized$proportions$cell.states)),
#                c("data","layers","scales","mapping","theme","coordinates","facet","plot_env","labels"))
# })

print("Done testing Celda CG.")