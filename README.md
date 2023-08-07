# wrappedpcm
Package wrappedpcm

A series of functions to allow for phylogenetic comparative methods on circular data

wrappedPIC (Wrapped Phylogenetically Independent Contrasts) is a modified version of Phylogenetically Independent Contrasts, but in circular space.

wrappedASE (Wrapped Ancestral State Estimation) is a modified version of Ancestral State Estimation in circular space.

fitWrappedBM (Fit Wrapped Brownian Motion) allows one to fit a wrapped Brownian Motion model to circular data.

##### Installing wrappedpcm

library(devtools)

install_github("markjuhn/wrappedpcm")

library(wrappedpcm)

##### Run examples

tree = ape::rcoal(10)

trait = runif(n = 10, min = 0, max = 1)

names(trait) = tree$tip.label


wrappedASE(trait = trait, phy = tree, lower_bound = 0, upper_bound = 1)


wrappedPIC(trait = trait, phy = tree, lower_bound = 0, upper_bound = 1)


fitWrappedBM(trait = trait, phy = tree, bounds = c(0, 1))


###### To delete package

remove.packages("wrappedpcm")
