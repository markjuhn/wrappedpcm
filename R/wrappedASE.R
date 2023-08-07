#' Ancestral State Estimation in circular space
#'
#' Estimates ancestral states using wrapped phylogenetic independent contrasts
#'
#' @param trait data vector for a single trait, with names matching tips in phy
#' @param phy phylogenetic tree of class phylo
#' @param lower_bound lower bound for trait in circular space
#' @param upper_bound upper bound for trait in circular space
#'
#' @author Mark Juhn
#' @import ape
wrappedASE = function(trait, phy, lower_bound, upper_bound){
  # NOTE TRAITS NEED TO BE LABELED
  if (is.null(names(trait)) == TRUE){
    print("NOTE: TRAIT VECTOR NOT NAMED, WILL SUPPLY NAMES FROM PHYLOGENY USING DEFAULT TIP ORDER")
    names(trait) = phy$tip.label
  }
  # Internal Functions

  # Find the minimum distance in pacman space
  find_min_dist_pac = function(x1, x2, lower_bound, upper_bound){
    non_pac_dist = abs(x1 - x2)
    pac_dist_lower = min(x1, x2) - lower_bound
    pac_dist_upper = upper_bound - max(x1, x2)
    pac_dist = pac_dist_lower + pac_dist_upper
    if (x1 > x2){
      if (non_pac_dist < pac_dist){
        # Return non_pac_dist (+)
        return(non_pac_dist)
      } else{
        # Return pac_dist (-)
        return(pac_dist * -1)
      }
    } else{
      if (non_pac_dist < pac_dist){
        # Return non_pac_dist (-)
        return(non_pac_dist * -1)
      } else {
        # Return pac_dist (+)
        return(pac_dist)
      }
    }

  }


  # Gets the correct sign for the ancestral node estimation for
  correct_sign = function(node_value, lb, ub){
    if (node_value > 0){
      return(lb + node_value)
    } else {
      return(ub - abs(node_value))
    }
  }

  # Start Regular Function
  trait = trait[phy$tip.label] # Reorder Traits
  nb.tip <- length(phy$tip.label)
  nb.node <- phy$Nnode
  # Set number of contrasts
  contr <- numeric(nb.node)
  # Set phenotype
  phenotype = c(trait, numeric(nb.node))
  # Reorder the length matrix
  reordered_matrix = phy$edge[order(phy$edge[,1], decreasing = T),]
  # Reorder the branch lengths
  reordered_bl = phy$edge.length[order(phy$edge[,1], decreasing = T)]
  ## Adapted from the ape code
  for (i in seq(from = 1, by = 2, length.out = nb.node)) {
    j <- i + 1
    anc <- reordered_matrix[i, 1]
    des1 <- reordered_matrix[i, 2]
    des2 <- reordered_matrix[j, 2]
    sumbl <- reordered_bl[i] + reordered_bl[j]
    ic <- anc - nb.tip
    distance_value <- find_min_dist_pac(x1 = phenotype[des1], x2 = phenotype[des2],
                                        lower_bound = lower_bound, upper_bound = upper_bound) # Find the minimum distance
    #print(abs(distance_value) == abs(phenotype[des1] - phenotype[des2]))

    contr[ic] <- distance_value
    # contr[ic] = phenotype[des1] - phenotype[des2] Normal non-pacman
    contr[ic] <- contr[ic]/sqrt(sumbl) # Standardizes contrast
    #phenotype[anc] <- (phenotype[des1]*reordered_bl[j] + phenotype[des2]*reordered_bl[i])/sumbl # Calculate trait value for ancestor
    # If you're using normal non-pacman space
    if (abs(distance_value) == abs(phenotype[des1] - phenotype[des2])){
      phenotype[anc] <- (phenotype[des1]*reordered_bl[j] + phenotype[des2]*reordered_bl[i])/sumbl # Calculate trait value for ancestor
    } else{
      #print("pacman space used")
      #print(c(abs(contr[ic]), abs(phenotype[des1] - phenotype[des2])))
      #print(distance_value)
      # Here use the pacman space estimation for ancestor node
      if (phenotype[des1] < phenotype[des2]){
        # then des1 is smaller then des2
        pacman_lower <- (phenotype[des1] - lower_bound) # get pacman_lower
        pacman_upper <- (upper_bound - phenotype[des2]) # get pacman_upper
        pacman_node <- (pacman_lower*reordered_bl[j] - pacman_upper*reordered_bl[i])/sumbl
        phenotype[anc] <- correct_sign(node_value = pacman_node, lb = lower_bound, ub = upper_bound) # if  positive, just use this, if negative, then subtract from upper bound
      } else {
        # Here, des2 is xlower
        pacman_lower <-  (phenotype[des2] - lower_bound) # get pacman_lower
        pacman_upper <- (upper_bound - phenotype[des1]) # get pacman_upper
        pacman_node <-  (phenotype[des2]*reordered_bl[i] - phenotype[des1]*reordered_bl[j])/sumbl
        phenotype[anc] <- correct_sign(node_value = pacman_node, lb = lower_bound, ub = upper_bound) # if  positive, just use this, if negative, then subtract from upper bound
      }
    }

    k <- which(reordered_matrix[, 2] == anc)
    reordered_bl[k] <- reordered_bl[k] + reordered_bl[i]*reordered_bl[j]/sumbl # What is important, rescales the branch lengths
  }
  # Name the contrasts
  names(contr) <- (nb.tip + 1):(nb.node + nb.tip)
  # Get the variance matric from just the pic function
  #pic_variance = pic(x = trait, phy = phy, var.contrasts = T, scaled = T)
  #new_pic_matrix = cbind(contrasts = contr, variance = pic_variance[,2], phenotype)
  new_pic_matrix = list(contrasts = contr, phenotype)
  asr = phenotype[(length(phy$tip.label)+1):length(phenotype)]
  return(as.numeric(asr))
}

