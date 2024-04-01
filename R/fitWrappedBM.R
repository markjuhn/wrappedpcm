#' Model fitting for continuous comparative data in circular space
#'
#' Fitting wrapped Brownian Motion model to phylogenetic trees
#'
#' @param phy a phylogenetic tree of class phylo
#' @param trait data vector for a single trait, with names matching tips in phy
#' @param bounds range defining bounds of circular space
#' @param bin_number number of discrete bins used to subdivide the continuous space
#' @export
#' @author Florian Boucher and Mark Juhn
#' @examples
#' tree = ape::rcoal(10)
#' trait = runif(n = 10, min = 0, max = 1)
#' names(trait) = tree$tip.label
#' fitWrappedBM(phy = tree, trait = trait, bounds = c(0, 1))
#'
fitWrappedBM = function(phy, trait, bounds, bin_number = 50){
  wrapped_likelihood=lnL_BBMV_circular_trait(tree = phy,trait,
                                             bounds = bounds,
                                             a = 0,b = 0,c = 0, Npts = bin_number)
  wrapped_fit = find.mle_FPK_circular_trait(model=wrapped_likelihood,method='Nelder-Mead',init.optim=NULL,safe=F)
  n_tips = length(phy$tip.label)
  correct_aic = 2*(4) * (n_tips/(n_tips - 4 - 1)) - 2*wrapped_fit$lnL
  correct_k = 4
  wrapped_fit$aic <- correct_aic
  wrapped_fit$k <- correct_k

  print_info_wbm = function(x){
    # Get Log Likelihood
    cat(paste("lnL =", round(x$lnL, 4)))
    # Get Sigma Squared
    sigsq = x$par$sigsq
    cat("\n sigsq =",formatC(sigsq, format = "e", digits = 4), "\n")
    # Get Root
    a_state = x$root
    z0 = a_state[a_state[,2] == max(a_state[,2]),1]
    cat("x0 =", z0, "\n")
  }
  a_state = wrapped_fit$root
  z0 = a_state[a_state[,2] == max(a_state[,2]),1]
  wrapped_fit$z0 = z0

  print_info_wbm(x = wrapped_fit)
  return(wrapped_fit)
}
