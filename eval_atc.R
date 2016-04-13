library(reshape2)
library(MESS)


get_pos <- function(atc4) {
  #calculates number of unique drug pairs sharing at least 1 level 4 atc code
  pos <- 0
  drugs <- names(atc4)
  
  for (i in seq_along(drugs)) {
    dr1 <- drugs[i]
    k1 <- atc4[[dr1]]
    
    for (dr2 in drugs[-(1:i)]) {
      k2 <- atc4[[dr2]]
      if (length(intersect(k1, k2)) > 0) pos <- pos + 1
    }
  }
  return(pos)
}

#----------------
get_neg <- function(atc4, pos){
  #calculated number of unique drug pairs not sharing at least 1 level 4 atc code
  
  pairs <- length(atc4) * (length(atc4) - 1) / 2
  neg <- pairs - pos
  return(neg)
}

#----------------

add_rate <- function(rates, denom, val=1) {
  #used by get_rates to add (val/denom) to tp and fp rates
  rate <- tail(rates, 1) + (val/denom)
  rates <- c(rates, rate)
  return(rates)
}

#---------------

get_rates <- function(sims, atc4, fp_max = 0.1) {
  #returns a list of tp and fp rates up to fp_max
  
  #only use sims with atc4 code
  sims <- sims[names(atc4), names(atc4)]
  sims[lower.tri(sims, diag = TRUE)] <- NA
  sims <- melt(sims, na.rm = TRUE)
  sims <- sims[order(sims$value, decreasing=T), ]
  
  #obtain number of unique drug pairs with a shared atc4
  tp_rates <- c(0)
  fp_rates <- c(0)
  pos <- get_pos(atc4)
  neg <- get_neg(atc4, pos)
  
  i <- 1
  while(tail(fp_rates, 1) < fp_max) {
    
    #check if next most similar pair share an atc4
    sub_sim <- sims[i, ]
    k1 <- atc4[[as.character(sub_sim$Var1)]]
    k2 <- atc4[[as.character(sub_sim$Var2)]]
    
    if (length(intersect(k1, k2)) > 0) {
      #add to tp (fp stays the same)
      tp_rates <- add_rate(tp_rates, pos, val=1)
      fp_rates <- add_rate(fp_rates, neg, val=0)
    } else {
      #add to fp (tp stays the same)
      tp_rates <- add_rate(tp_rates, pos, val=0)
      fp_rates <- add_rate(fp_rates, neg, val=1)
    }
    i <- i + 1
  }
  auc <- auc(fp_rates, tp_fates)
  rates <- list(fp=fp_rates, tp=tp_rates, auc=auc)
  return(rates)
}