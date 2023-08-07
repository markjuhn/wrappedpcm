#' Model fitting for continuous comparative data in circular space
#'
#' Fitting wrapped Brownian Motion model to phylogenetic trees
#'
#' @param phy a phylogenetic tree of class phylo
#' @param trait data vector for a single trait, with names matching tips in phy
#' @param bounds range defining bounds of circular space
#'
#' @author Florian Boucher and Mark Juhn
#' @import BBMV
fitWrappedBM = function(phy, trait, bounds){
  wrapped_likelihood=lnL_BBMV_circular_trait(tree = phy,trait,
                                             bounds = bounds,
                                             a = 0,b = 0,c = 0, Npts = 50)
  wrapped_fit = find.mle_FPK_circular_trait(model=wrapped_likelihood,method='Nelder-Mead',init.optim=NULL,safe=F)
  n_tips = length(phy$tip.label)
  correct_aic = 2*(4) * (n_tips/(n_tips - 4 - 1)) - 2*wrapped_fit$lnL
  correct_k = 4
  wrapped_fit$aic <- correct_aic
  wrapped_fit$k <- correct_k
  return(wrapped_fit)
}
