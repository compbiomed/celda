# Title     : TODO
# Objective : TODO
# Created by: sarahleinicke
# Created on: 3/14/18

# Reference cell clustering script
source('celda_C.R')
source('celda_G.R')
source('celda_CG.R')
source('celda_functions.R')
source('split_clusters.R')
source('s3_generics.R')
library(profvis)


#################################
#    PROFVIS CELDA CG SimCells  #
#################################

# sim_counts_CG = simulateCells.celda_CG("celda_CG", S=10,
#     C.Range=c(100, 500), N.Range=c(500,5000),
#     G=1000, K=5, L=9)
#
# print('Dimensions - simulateCells_celda_CG')
# print(dim(sim_counts_CG$counts))
#
# print('Profiling EDITED CELDA CG')
#
# p <- profvis({
#     celda_CG(counts = sim_counts_CG$counts,
#     K=5, L=9, max.iter = 15, random.state.order=FALSE)
# })
#
# print(p)


##################################
#    PROFVIS CELDA CG RDA DATA   #
##################################

# load("../tests/celdaCGsim.rda")
# p <- profvis(celda_CG(counts = celdaCG.sim$counts,
#     K = celdaCG.sim$K, L=celdaCG.sim$L, max.iter = 15,
#     random.state.order = FALSE))
# print(p)



#################################
#    PROFVIS CELDA G SimCells  #
#################################

sim_counts_G = simulateCells.celda_G()
print('Dimensions - simulateCells.celda_G()')
print(dim(sim_counts_G$counts))

print('Profiling EDITED CELDA G')
p <- profvis(celda_G(counts = sim_counts_G$counts,
    L = 5, beta=1, gamma = 1,
    stop.iter = 10, max.iter= 200, split.on.iter=10, split.on.last=TRUE,
    random.state.order=FALSE, count.checksum=NULL, seed=12345,
    y.init = NULL, process.counts=TRUE, logfile=NULL))

print(p)