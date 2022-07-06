allocate_dhondt <- function(votes, n) {
  seats <- integer(length(votes))
  score <- votes
  for (m in 1:n) {
    nextparty <- which.max(score)
    seats[nextparty] <- seats[nextparty] + 1
    score[nextparty] <- votes[nextparty] / (seats[nextparty] + 1)
  }
  names(seats) <- names(votes)
  seats
}

apply_threshold <- function(x, threshold, rescale = TRUE, all_under_action = c("ignore_threshold","zero","error")) {
  
  over_threshold <- (x >= threshold)
  y <- x
  
  if (!any(over_threshold)) {
    if (all_under_action == "ignore_threshold") {
      over_threshold[] <- TRUE
    }
    else if (all_under_action == "error") {
      stop("All ",length(x)," vote counts were below threshold ",thres,".")
    } 
    else if (all_under_action == "zero") {
      y[] <- 0
      return(y)
    }
  }
  
  y[!over_threshold] <- 0
  if (rescale) {
    y <- y * sum(x) / sum(y)
  }
  y
}

allocate_largest_remainders <- function(x, ns) {
  p <- x/sum(x)
  seats <- floor(p*ns)
  remainders <- p*ns - seats
  seats_remaining <- ns - sum(seats)
  has_largest_remainders <- (rank(-remainders, ties.method = "random") <= seats_remaining)
  seats[has_largest_remainders] <- seats[has_largest_remainders] + 1
  seats
}

allocate_largest_remainders_threshold <- function(x, ns, threshold, ...) {
  allocate_largest_remainders(apply_threshold(x, threshold, ...), ns)
}
