library(plyr)

rnsp <- nsp_robust <- function(x, M = 1000, thresh = NULL, alpha = 0.1, thresh.type = c("bern", "gauss", "sim"), thresh.sim.times = 1000, deg = 0, max.length = 3000, overlap = FALSE, zeros = TRUE, mixing = c("butterfly", "all.pairs", "fixed.sep")[1], max.bandwidth = 10) {
	
	# Main function. Some of the parameters should be ignored as they relate to the piecewise-linear version only (currently not backed up by anything written up).
	# The only parameters relevant to the piecewise-constant version are:
	#
	# x - the data
	# M - number of intervals
	# thresh - threshold for the deviation measure; best left as NULL but specified indirectly by the significance parameter alpha and thresh.type
	# alpha - significance parameter
	# thresh.type - how to calculate the threshold; "bern" is the symmetric Bernoulli threshold from Kabluchko and Wang (2014), "gauss" is 0.9 * corresponding Gaussian threshold
	#   (a very good approximation of "bern") and "sim" determines the threshold by simulation but there is no need to use this option
	# thresh.sim.times - if thresh.type == "sim" then how many sample paths to simulate
	# deg - leave as 0; 1 refers to the piecewise-linear version (undocumented for now)
	# max.length - maximum length of an interval of significance returned
	# overlap - FALSE for no overlap; TRUE for overlap as in the paper
	# zeros - leave as TRUE
	# the remaining parameters - please ignore
	
	n <- length(x)

	if (is.null(thresh)) {
		
		if (is.null(alpha)) alpha <- 0.1
		if (length(thresh.type) == 3) thresh.type <- thresh.type[1]
		if (thresh.type == "bern") thresh <- thresh_kab_bern(n, alpha)
		else if (thresh.type == "gauss") thresh <- 0.9 * thresh_kab(n, alpha)
			else {
				
				msn <- sim_multi_sup_norm(n, thresh.sim.times)
				thresh <- as.numeric(quantile(msn, 1-alpha))
				
			}
		
	}
	
	res <- iter_random_checks_scan_signed(c(1, n), x, M, thresh, max.length, overlap, 0, deg, zeros, mixing, max.bandwidth)
	
	intervals <- data.frame(t(order_chron(res)))
	colnames(intervals) <- c("starts", "ends", "values")
	
	list(x=x, intervals=intervals, threshold.used=thresh)
	
}


thresh_kab <- function(n, alpha = 0.1, method = "asymp") 
{
  an <- sqrt(2 * log(n)) + (1/2 * log(log(n)) + log(0.8197466/2/sqrt(pi)))/sqrt(2 * log(n))
  bn <- 1/sqrt(2 * log(n))
  if (method == "bound") 
    beta <- alpha/2
  else if (method == "asymp") 
    beta <- 1 - sqrt(1 - alpha)
  an + bn * log(1/log(1/(1 - beta)))
}


thresh_kab_bern <- function(n, alpha = 0.1) {
	
	an <- sqrt(2 * log(n*log(n)^(-1/2)))
	
	tau <- -log(-1/(2*0.2740311) * log(1 - alpha))
	
	an + tau / an
	
}


iter_random_checks_scan_signed <- function(ind, x, M, thresh, max.length, overlap = FALSE, buffer = 0, deg = 0, zeros = TRUE, mixing = c("butterfly", "all.pairs", "fixed.sep")[1], max.bandwidth = 10) {
	
	s <- ind[1]
	e <- ind[2]
	n <- e - s + 1
	x <- x[s:e]
		
	if (n > 1) {
		
		next.int <- random_checks_scan_signed_2stage(c(1,n), x, M, thresh, max.length, deg, zeros, mixing, max.bandwidth)
		
		if (!is.na(next.int$selected.ind))  {
		
			if (!overlap) {
			
			if (next.int$selected.val[1,1]-buffer >= 2) left <- iter_random_checks_scan_signed(c(1, next.int$selected.val[1,1]-buffer), x, M, thresh, max.length, overlap, buffer, deg, zeros, mixing, max.bandwidth) else left <- matrix(NA, 3, 0)
			if (n - next.int$selected.val[1,2]-buffer >= 1) {
				
				right <- iter_random_checks_scan_signed(c(next.int$selected.val[1,2]+buffer, n), x, M, thresh, max.length, overlap, buffer, deg, zeros, mixing, max.bandwidth)
				if (dim(right)[2]) right <- right + c(rep(next.int$selected.val[1,2]-1+buffer, 2), 0)
				
			}
			else right <- matrix(NA, 3, 0)
			}
			
			else {
				

			if (floor(mean(next.int$selected.val[1,1:2]))-buffer >= 2) left <- iter_random_checks_scan_signed(c(1, floor(mean(next.int$selected.val[1,1:2]))-buffer), x, M, thresh, max.length, overlap, buffer, deg, zeros, mixing, max.bandwidth) else left <- matrix(NA, 3, 0)
			if (n - floor(mean(next.int$selected.val[1,1:2]))-buffer >= 2) {
				right <- iter_random_checks_scan_signed(c(floor(mean(next.int$selected.val[1,1:2]))+1+buffer, n), x, M, thresh, max.length, overlap, buffer, deg, zeros, mixing, max.bandwidth)
				if (dim(right)[2]) right <- right + c(rep(floor(mean(next.int$selected.val[1,1:2]))+buffer, 2), 0)
				
			}
			else right <- matrix(NA, 3, 0)
				
				
			}
			
			
			return(cbind(t(next.int$selected.val), left, right))
			
			
		}
		
		else(return(matrix(NA, 3, 0)))
		
		
	}

	else(return(matrix(NA, 3, 0)))
	
}


random_checks_scan_signed_2stage <- function(ind, x, M, thresh, max.length, deg = 0, zeros = TRUE, mixing = c("butterfly", "all.pairs", "fixed.sep")[1], max.bandwidth = 10) {
		
	s1 <- random_checks_scan_signed(ind, x, M, thresh, max.length, deg, zeros, mixing, max.bandwidth)
		
	if (!is.na(s1$selected.ind)) {
		
		s <- s1$selected.val[1,1] + ind[1] - 1
		e <- s1$selected.val[1,2] + ind[1] - 1
		
		s2 <- random_checks_scan_signed(c(s,e), x, M, thresh, max.length, deg, zeros, mixing, max.bandwidth)

		if (!is.na(s2$selected.ind)) {
			
			replacement <- s2$selected.val
			replacement[1,1:2] <- replacement[1,1:2] + s - ind[1]
			s1$selected.val <- replacement
			
		}

	}
	
	s1	
	
}


random_checks_scan_signed <- function(ind, x, M, thresh, max.length, deg = 0, zeros = TRUE, mixing = c("butterfly", "all.pairs", "fixed.sep")[1], max.bandwidth = 10) {

	s <- ind[1]
	e <- ind[2]
	n <- e - s + 1
	
	if (n > 1) {
		
		x <- x[s:e]

		M <- min(M, (n-1)*n/2)

		ind <- grid_intervals_sorted(n, M)      # record where this function can be found

		M <- dim(ind)[2]

		res <- matrix(0, M, 3)

		res[,1:2] <- t(ind)

		zero.check <- TRUE
		j <- 1
		
		while (zero.check && (j <= M)) {
			
			res[j,3] <- check_interval_robust_sup_norm(res[j,1:2], x, thresh, max.length, deg, zeros, mixing, max.bandwidth)
			zero.check <- (res[j,3] == 0)
			j <- j + 1
			
		}

		if (zero.check) {
			
			selected.ind <- NA
			selected.val <- matrix(0, 0, 3)
			
		}

		else {
			
			selected.ind <- j-1
			selected.val <- res[selected.ind,,drop=FALSE]
			
		}


	}

	else {
		
		selected.val <- matrix(0, 0, 3)
		selected.ind <- NA
		M <- 0
		
	}
	
	list(selected.ind = selected.ind, selected.val = selected.val, M.eff=max(1, M))

}




check_interval_robust_sup_norm <- function(ind, x, thresh, max.length, deg = 0, zeros = TRUE, mixing = c("butterfly", "all.pairs", "fixed.sep")[1], max.bandwidth = 10) {
	
	if (ind[2] - ind[1] + 1 > max.length) return(0) else {
		
		if (deg == 0) sup.nm <- robust_sup_norm_rv_wwozr(x[ind[1]:ind[2]], zeros)$sup.norm
			else if (deg == 1) sup.nm <- robust_sup_norm_lin_combined(x[ind[1]:ind[2]], mixing, max.bandwidth)$sup.norm
		
		sup.nm * (sup.nm > thresh)
	
	}

}



robust_sup_norm_rv_wwozr <- function(data, zeros = TRUE) {
	
	# rv stands for "rank version"
	# wwozr stands for "with or without zeros"
	# this version uses rank rather than order
	# the user can set zeros=FALSE if they are sure there are no masses at any local medians (e.g. if the distribution of the data is continuous)
	# requires plyr::mapvalues !!
	

	n <- length(data)
	
	if (n <= 1) {
		
		sup.norm <- 0
		d <- (zeros + 1) * n + 1
		left.max <- right.max <- comb.max <- rep(0, d)
		centre <- data
		sig <- data		
		
	} else {	

		data.r <- rank(data, ties.method = "min")
		data.r.s <- sort(unique(data.r))
		how.many <- length(data.r.s)
		data.sr <- plyr::mapvalues(data.r, data.r.s, 1:how.many)
	
		if (zeros) {

			pseudo.data <- matrix(1, n, 2*how.many+1)
	
			for (i in 1:n) {
				pseudo.data[i, (2*data.sr[i]+1):(2*how.many+1)] <- -1
				pseudo.data[i, 2*data.sr[i]] <- 0
			}

		}
		else {
			
			pseudo.data <- matrix(1, n, how.many+1)
	
			for (i in 1:n) pseudo.data[i, (data.sr[i]+1):(how.many+1)] <- -1
			
		}
		
		left.cumsum	<- apply(pseudo.data, 2, cumsum)
		left.max <- apply(abs(left.cumsum) / sqrt(1:n), 2, max)

		right.cumsum <- apply(pseudo.data[n:1,], 2, cumsum)
		right.max <- apply(abs(right.cumsum) / sqrt(1:n), 2, max) 

		comb.max <- pmax(left.max, right.max)
	
		which.comb.min <- which.min(comb.max)
	
		sup.norm <- comb.max[which.comb.min]

		if (zeros) {
			
			which.comb.trim.min <- which.min(comb.max[2:(2*how.many)])
			centre.ind.up <- data.r.s[floor((which.comb.trim.min + 1)/2)]
			centre.ind.down <- data.r.s[ceiling((which.comb.trim.min + 1)/2)]
			centre <- 1/2 * (data[which(data.r == centre.ind.up)][1] + data[which(data.r == centre.ind.down)][1])
			
		}
		else {
			
			if (how.many > 1) {
				which.comb.trim.min <- which.min(comb.max[2:how.many])
				centre.ind.up <- data.r.s[floor(which.comb.trim.min + 1/2)]
				centre.ind.down <- data.r.s[ceiling(which.comb.trim.min + 1/2)]
				centre <- 1/2 * (data[which(data.r == centre.ind.up)][1] + data[which(data.r == centre.ind.down)][1])
				
			}
			else centre <- data[1]
			
			
		}

		sig <- rep(centre, n)
		
	}
	
	return(list(sup.norm = sup.norm, left.max = left.max, right.max = right.max, comb.max = comb.max, centre = centre, sig = sig))

}


robust_sup_norm_lin_combined <- function(data, mixing = c("butterfly", "all.pairs", "fixed.sep")[1], max.bandwidth = 10) {
  
  #	This function combines robust_sup_norm_lin_gen2p (for mixing == "all.pairs" or "fixed.sep") and robust_sup_norm_lin_2p (for mixing == "butterfly")
  #	To be used in check_interval_robust_sup_norm
  #	max.bandwidth only relevant if mixing == "all.pairs" or "fixed.sep"
  
  n <- length(data)
  
  if (n <= 1) {
    
    sup.norm <- 0
    d <- 1
    middle.max <- left.max <- right.max <- comb.max <- rep(0, d)
    ret <- list(sup.norm = sup.norm, middle.max = middle.max, left.max = left.max, right.max = right.max, comb.max = comb.max)
    
  } else if (mixing == "all.pairs") {
    
    bandwidth = min(max.bandwidth, floor(length(data)/2))		
    g <- matrix(0, bandwidth^2, 2)
    g[,1] <- rep(1:bandwidth, bandwidth)
    g[,2] <- rep((n-bandwidth+1):n, each = bandwidth)
    ret <- robust_sup_norm_lin_gen2p(data, g)
    
  } else if (mixing == "fixed.sep") {
    
    bandwidth = min(max.bandwidth, floor(length(data)/2))				
    g <- matrix(0, bandwidth, 2)
    g[,1] <- 1:bandwidth
    g[,2] <- (n-bandwidth+1):n
    ret <- robust_sup_norm_lin_gen2p(data, g)
    
  } else if (mixing == "butterfly")
    ret <- robust_sup_norm_lin_2p(data)
  
  ret	
  
}

robust_sup_norm_lin_gen2p <- function(data, point.pairs) {
  
  #	This function checks all linear fits through point pairs specified in point.pairs, therefore is more general than robust_sup_norm_lin_2p
  #	point.pairs is an how.many x 2 matrix
  
  n <- length(data)
  
  if (n <= 1) {
    
    sup.norm <- 0
    d <- 1
    middle.max <- left.max <- right.max <- comb.max <- rep(0, d)
    
  } else {
    
    how.many <- dim(point.pairs)[1]
    
    pseudo.data <- matrix(1, n, how.many)
    
    for (i in 1:how.many) {
      
      try.fit <- straight_line(point.pairs[i,1], data[point.pairs[i,1]], point.pairs[i,2], data[point.pairs[i,2]], 1:n)
      pseudo.data[,i] <- sign(data - try.fit)
      
    }
    
    
    reorder.x <- order(abs((1:n) - (n+1)/2))
    
    middle.cumsum <- apply(pseudo.data[reorder.x,,drop=FALSE], 2, cumsum)
    middle.max <- apply(abs(middle.cumsum) / sqrt(1:n), 2, max)
    
    left.cumsum <- apply(pseudo.data, 2, cumsum)
    left.max <- apply(abs(left.cumsum) / sqrt(1:n), 2, max)
    
    right.cumsum <- apply(pseudo.data[n:1,,drop=FALSE], 2, cumsum)
    right.max <- apply(abs(right.cumsum) / sqrt(1:n), 2, max)
    
    
    comb.max <- pmax(middle.max, left.max, right.max)
    
    which.comb.min <- which.min(comb.max)
    
    sup.norm <- comb.max[which.comb.min]
    
  }
  
  return(list(sup.norm = sup.norm, middle.max = middle.max, left.max = left.max, right.max = right.max, comb.max = comb.max))
  
}

robust_sup_norm_lin_2p <- function(data) {
  
  n <- length(data)
  
  if (n <= 1) {
    
    sup.norm <- 0
    d <- 1
    middle.max <- left.max <- right.max <- comb.max <- rep(0, d)
    
  } else {
    
    how.many <- floor(n/2)
    
    pseudo.data <- matrix(1, n, how.many)
    
    for (i in 1:how.many) {
      
      try.fit <- straight_line(i, data[i], n-i+1, data[n-i+1], 1:n)
      pseudo.data[,i] <- sign(data - try.fit)
      
    }
    
    
    reorder.x <- order(abs((1:n) - (n+1)/2))
    
    middle.cumsum <- apply(pseudo.data[reorder.x,,drop=FALSE], 2, cumsum)
    middle.max <- apply(abs(middle.cumsum) / sqrt(1:n), 2, max)
    
    left.cumsum <- apply(pseudo.data, 2, cumsum)
    left.max <- apply(abs(left.cumsum) / sqrt(1:n), 2, max)
    
    right.cumsum <- apply(pseudo.data[n:1,,drop=FALSE], 2, cumsum)
    right.max <- apply(abs(right.cumsum) / sqrt(1:n), 2, max)
    
    
    comb.max <- pmax(middle.max, left.max, right.max)
    
    which.comb.min <- which.min(comb.max)
    
    sup.norm <- comb.max[which.comb.min]
    
  }
  
  return(list(sup.norm = sup.norm, middle.max = middle.max, left.max = left.max, right.max = right.max, comb.max = comb.max))
  
}



straight_line <- function(x1, y1, x2, y2, range) {
  
  a <- (y2 - y1) / (x2 - x1)
  
  a * (range - x1) + y1
  
}




sim_multi_sup_norm <- function(n, N, prob0 = 0, prob1 = 1/2) {
	
	max.use.sample <- max.abs.use.sample <- rep(0, N)
	
	for (i in 1:N) {
		
		eps <- sample(c(-1, 0, 1), n, TRUE, c(1-prob0-prob1, prob0, prob1))
				
		cur.max.use <- multi_sup_norm(eps)
		
		max.use.sample[i] <- cur.max.use
		
	}
	
	max.use.sample	
		
}


multi_sup_norm <- function(eps) {
	
	
	n <- length(eps)
	
	eps.cum <- c(0, cumsum(eps))
		
	max.use <- max.abs.use <- 0
	
	for (i in 0:(n-1)) for (j in (i+1):n) {
		
		scan.stat <- (eps.cum[j+1] - eps.cum[i+1]) / sqrt(j-i)
						
		if (abs(scan.stat) > max.abs.use) max.abs.use <- abs(scan.stat)
				
	}
	
	max.abs.use
	
}



grid_intervals_sorted <- function(n, M) {
  
  if ((n==2) || (M == 0)) ind <- matrix(c(1, n), 2, 1)
  
  else if (M >= (n-1)*n/2) ind <- all_intervals_sorted(n)
  
  else {
    k <- 1
    while (k*(k-1)/2 < M) k <- k+1
    ind2 <- all_intervals_sorted(k)
    ind2.mx <- max(ind2)
    ind <- round((ind2 - 1) * ((n-1) / (ind2.mx-1)) + 1)
  }	
  
  ind	
}


all_intervals_sorted <- function(n) {
  
  d <- all_intervals_flat(n)
  d.ord <- order(d[2,] - d[1,])
  d[,d.ord, drop=FALSE]
  
}


all_intervals_flat <- function(n) {
  
  if (n == 2) ind <- matrix(1:2, 2, 1) else {
    M <- (n-1)*n/2	
    ind <- matrix(0, 2, M)
    ind[1,] <- rep(1:(n-1), (n-1):1)
    ind[2,] <- 2:(M+1) - rep(cumsum(c(0, (n-2):1)), (n-1):1)
  }
  ind
  
}


order_chron <- function(nsp.obj) {
  
  # Order intervals of significance chronologically.
  # nsp.obj - quantity returned by one of the nsp_* functions.
  
  d <- dim(nsp.obj)
  if (d[2]) {
    nsp.obj.ord <- order(nsp.obj[1,])
    nsp.obj <- nsp.obj[,nsp.obj.ord]
  }	
  
  nsp.obj
  
}

rnsp_pointwise <- function(rnsp_obj) {
  
  # Finds pointwise change-point location estimates given an object returned by rnsp.
  
  int <- subset(rnsp_obj$intervals, select = c(starts, ends, values))
  
  d <- dim(int)
  
  midpoints <- floor((int$starts + int$ends)/2)
  
  max_sign_cusums <- rep(0, d[1])
  
  if (d[1]) for (i in 1:d[1]) {
    
    max_sign_cusums[i] <- argmax_sign_cusum(rnsp_obj$x[int$starts[i]:int$ends[i]]) + int$starts[i] - 1			
    
  }	
  
  int <- cbind(int, midpoints, max_sign_cusums)
  
  list(x = rnsp_obj$x, intervals = int, threshold.used = rnsp_obj$threshold.used)
  
}

sign_cusum <- function(x) {
  
  y <- sign(x - median(x))
  inner_prod_iter(y)
  
} 

argmax_sign_cusum <- function(x) {
  
  sc <- sign_cusum(x)
  max.ind <- which(sc == max(sc))
  median(max.ind)
  
}

inner_prod_iter <- function(x) {
  
  m <- length(x)
  z <- cumsum(x)
  
  ip <- sqrt(((m-1):1) / m / (1:(m-1))) * z[1:(m-1)] - sqrt((1:(m-1)) / m / ((m-1):1)) * (z[m] - z[1:(m-1)])
  
  abs(ip)	
}


draw_rects <- function(nsp.obj, yrange, density = 10, col = "red", x.axis.start = 1) {
  
  # Draw intervals of significance, as shaded rectangular areas, on the current plot.
  # nsp.obj - quantity returned by one of the nsp_* functions.
  # yrange - vector of length two specifying the (lower, upper) vertical limit of the rectangles.
  # density - density of the shading; try using 10 or 20.
  # col - colour of the shading.
  # x.axis.start - time index the x axis stars from.
  
  d <- dim(nsp.obj$intervals)
  if (d[1]) for (i in 1:d[1]) {
    
    rect(nsp.obj$intervals[i,1]+x.axis.start-1, yrange[1], nsp.obj$intervals[i,2]+x.axis.start-1, yrange[2], density=density, col=col)
    
    
  }
  
}


###########################################################
### Computing multiscale norms for all sub-intervals [s,e]
###########################################################

all_multiscale_norms_combined <- function(x) {
  
  # Computes sign-multiresolution sup-norms for all subintervals of [1,T], where
  # T is the length of x. The norm for interval [s,e] is in the [s,e]th element of the output
  # matrix. Warning: this is an O(T^3) procedure. Execution times in seconds on our machine for a
  # rnorm(n) input, for n in the left column:
  #       [,1]    [,2]
  # [1,]  100   0.077
  # [2,]  200   0.575
  # [3,]  300   1.819
  # [4,]  400   5.258
  # [5,]  500   7.755
  # [6,]  600  12.820
  # [7,]  700  20.530
  # [8,]  800  28.992
  # [9,]  900  40.442
  # [10,] 1000  63.058
  # [11,] 1100  75.089
  # [12,] 1200  94.862
  # [13,] 1300 131.075
  # [14,] 1400 155.171
  # [15,] 1500 179.712
  # [16,] 1600 218.301
  # [17,] 1700 260.205
  # [18,] 1800 308.756
  # [19,] 1900 358.544
  # [20,] 2000 420.255
  # [21,] 2100 606.606
  # [22,] 2200 708.727
  # [23,] 2300 779.498
  # [24,] 2400 872.683
  # [25,] 2500 980.730
  
  # NOTE: if norms up to a certain max value of e-s are required, please use all_multiscale_norms_v2 instead.
  
  if (length(x) <= 2050) all_multiscale_norms_v2(x) else all_multiscale_norms(x)
  
}


all_multiscale_norms <- function(x) {
  
  x_rol <- rank_order_levels(x)
  n <- length(x)
  multipliers <- partial_sums <- matrix(0, n, n)
  for (i in 0:(n-1)) {
    B <- matrix(c(1:(n-i), (i+1):n), n-i, 2)
    multipliers[B] <- 1/sqrt(i+1)
    partial_sums[B] <- sqrt(i+1)
  }
  
  norms <- final_norms <- partial_sums
  
  start_of_next_level <- 1
  
  for (i in 1:x_rol$no_of_levels) {
    
    locations <- x_rol$x_o[start_of_next_level:(start_of_next_level + x_rol$level_lengths[i] - 1)]
    
    for (l in 1:2) {
      
      for (j in 1:x_rol$level_lengths[i]) {
        partial_sums[1:locations[j], locations[j]:n] <- partial_sums[1:locations[j], locations[j]:n] - multipliers[1:locations[j], locations[j]:n]
      }
      for (j in 1:x_rol$level_lengths[i]) {
        norms[1:locations[j], locations[j]:n] <- abs(partial_sums[1:locations[j], locations[j]:n])
      }
      for (k in 1:(n-1)) {
        start_x <- max(1, locations[1]-k)
        start_y <- start_x + k
        end_y <- min(locations[x_rol$level_lengths[i]]+k, n)
        end_x <- end_y - k
        B <- matrix(c(start_x:end_x, start_y:end_y), end_x - start_x + 1, 2)
        B_left <- matrix(c(start_x:end_x, (start_y-1):(end_y-1)), end_x - start_x + 1, 2)
        B_down <- matrix(c((start_x+1):(end_x+1), start_y:end_y), end_x - start_x + 1, 2)
        norms[B] <- pmax(norms[B], norms[B_left], norms[B_down])
      }
      for (j in 1:x_rol$level_lengths[i]) {
        final_norms[1:locations[j], locations[j]:n] <- pmin(final_norms[1:locations[j], locations[j]:n], norms[1:locations[j], locations[j]:n])
      }
      
    }
    
    start_of_next_level <- start_of_next_level + x_rol$level_lengths[i]
    
  }
  
  final_norms
  
}

rank_order_levels <- function(x) {
  
  x_r <- rank(x, ties = "min")
  x_o <- order(x)
  ordered_ranks <- x_r[x_o]
  level_lengths <- rle(ordered_ranks)$lengths
  no_of_levels <- length(level_lengths)
  list(x_r = x_r, x_o = x_o, level_lengths = level_lengths, no_of_levels = no_of_levels)
  
}

all_multiscale_norms_v2 <- function(x, max_length = length(x)) {
  
  x_rol <- rank_order_levels(x)
  n <- length(x)
  data <- matrix(1, n, 2 * x_rol$no_of_levels + 1)
  
  start_of_next_level <- 1
  
  for (i in 1:x_rol$no_of_levels) {
    
    locations <- x_rol$x_o[start_of_next_level:(start_of_next_level + x_rol$level_lengths[i] - 1)]
    
    data[, 2 * i] <- data[, 2 * i + 1] <- data[, 2 * i - 1]
    
    data[locations, 2 * i] <- 0
    data[locations, 2 * i + 1] <- -1
    
    start_of_next_level <- start_of_next_level + x_rol$level_lengths[i]
    
  }
  
  data_cumsum <- rbind(rep(0, 2 * x_rol$no_of_levels + 1), apply(data, 2, cumsum))
  
  current_stage_norms <- abs(data)
  
  final_norms <- matrix(0, n, n)
  
  for (i in 3:min(n+1, max_length+1)) {
    
    current_partial_sums <- (data_cumsum[i:(n+1),,drop=FALSE] - data_cumsum[1:(n-i+2),,drop=FALSE])/sqrt(i-1)
    
    current_stage_norms[1:(n-i+2),] <- pmax(current_stage_norms[1:(n-i+2),], current_stage_norms[2:(n-i+3),], abs(current_partial_sums))
    
    B <- matrix(c(1:(n-i+2), (i-1):n), n-i+2, 2)
    
    final_norms[B] <- apply(current_stage_norms[1:(n-i+2),,drop=FALSE], 1, min)
    
  }
  
  final_norms
  
}
