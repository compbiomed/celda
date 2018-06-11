#celda_G
library(testthat)
library(Rtsne)
library(doParallel)

context("Testing celda_G")

source('R/celda_functions.R')
source('R/celda_list.R')
source('R/celda_G.R')
source('R/diffExp.R')
source('R/split_clusters.R')
source('R/s3_generics.R')
source('R/feature_selection.R')

load("tests/celdaGsim.rda")
load("tests/celdaG.rda")

model_G = getModel(celdaG.res, L = 5)[[1]]
#factorized_G <- factorizeMatrix(celda.mod = model_G, counts = celdaG.sim$counts)
factorized_G <- factorizeMatrix.celda_G(celda.mod = model_G, counts = celdaG.sim$counts)
counts.matrix_G <- celdaG.sim$counts

#Making sure getModel if functioning correctly
test_that(desc = "Sanity checking getModel",{
  expect_equal(celdaG.res$content.type, class(model_G))
})

#Making sure relationship of counts vs proportions is correct in factorize matrix
test_that(desc = "Checking factorize matrix, counts vs proportions",{
  expect_equal(TRUE,all(factorized_G$counts$sample.states/sum(factorized_G$counts$sample.states)
                        == factorized_G$proportions$sample.states))
})

#Checking dimension of factorize matrix
test_that(desc = "Checking factorize matrix dimension size",{
  expect_equal(5, ncol(factorized_G$proportions$gene.states))
})


# Ensure calculateLoglikFromVariables calculates the expected values
test_that(desc = "calculateLoglikFromVariables.celda_G returns correct output for various params", {
  expect_lt(calculateLoglikFromVariables(counts = celdaG.sim$counts, y = celdaG.sim$y,
                                         L = celdaG.sim$L, delta = 1, gamma = 1, beta = 1, 
                                         model="celda_G"), 0)
})

#normalizeCounts
test_that(desc = "Making sure normalizeCounts doesn't change dimensions of counts matrix_G",{
  norm.counts <- normalizeCounts(counts.matrix_G)
  expect_equal(dim(norm.counts),dim(counts.matrix_G))
  expect_equal(rownames(norm.counts),rownames(counts.matrix_G))
  expect_equal(colnames(norm.counts),colnames(counts.matrix_G))
})

#recodeClusterY
test_that(desc = "Checking recodeClusterY gives/doesn't give error",{
  expect_error(recodeClusterY(celda.mod = model_G, from = NULL, to = ))
  expect_error(recodeClusterY(celda.mod = model_G, from=c(1,2,3,4,5), to = c(1,2,3,4,6)))
  new.recoded = recodeClusterY(celda.mod = model_G, from=c(1,2,3,4,5), to = c(5,4,3,2,1))
  expect_equal(model_G$y == 1, new.recoded$y == 5)
})

#compareCountMatrix
test_that(desc = "Checking CompareCountMatrix",{
  expect_true(compareCountMatrix(count.matrix = celdaG.sim$counts, celda.obj = model_G))
})

#distinct_colors
test_that(desc = "Checking distinct_colors",{
  expect_equal(distinct_colors(2), c("#FF9999","#99FFFF"))
})


###renderCeldaHeatmap###
# test_that(desc = "Checking renderCeldaHeatmap output",{
#   expect_equal(names(renderCeldaHeatmap(counts = celdaG.sim$counts, z = model_G$z, y = model_G$y)),
#                c("tree_row","tree_col","kmeans","gtable"))
# })

##feature_selection.R##
#topRank
test_that(desc = "Checking topRank function",{
  expect_equal(names(topRank(fm = factorized_G$proportions$gene.states)),
               c("index","names"))
})

#stateHeatmap
# test_that(desc = "Checking stateHeatmap to see if it runs",{
#   expect_equal(names(stateHeatmap(celdaG.sim$counts, celda.mod = model_G)),
#                c("tree_row","tree_col","kmeans","gtable"))
# })

##celda_G.R##
test_that(desc = "Making sure celda_G runs without errors",{
  #celdaG.res <- celda_G(counts = celdaG.sim$counts, model = "celda_G", nchains = 2, L = 5, max.iter = 15)
  #expect_true(class(celdaG.res) == "celda_list")  # Only best chain returned by default

  celdaG.res <- celda_G(counts = celdaG.sim$counts, L = 5, max.iter = 15, random.state.order = FALSE)
  expect_true(class(celdaG.res) == "celda_G")

})

#plotDrState
# test_that(desc = "Checking plotDrState",{
#   celda.tsne <- celdaTsne(counts = celdaG.sim$counts,max.iter = 50,celda.mod=model_G)
#   expect_equal(names(plotDrState(dim1 = celda.tsne[,1], dim2 = celda.tsne[,2],matrix = factorized_G$proportions$cell.states)),
#                c("data","layers","scales","mapping","theme","coordinates","facet","plot_env","labels"))
# })

print("Done testing Celda G.")