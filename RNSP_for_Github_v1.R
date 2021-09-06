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
	# overlap - FALSE for no overlap; TRUE for overlap as in Sections 4 and 5 of the paper
	# zeros - leave as TRUE
	# the remaining parameters - please ignore
	
	# NOTE 1: functions: thresh_kab, order_chron, grid_intervals_sorted are defined in the NSP code (the same Github repo)
	# NOTE 2: code requires R package plyr

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
	
	list(intervals=intervals, threshold.used=thresh)
	
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



sim_multi_sup_norm <- function(n, N, prob0 = 0, prob1 = 1/2) {
	
	max.use.sample <- max.abs.use.sample <- rep(0, N)
	
#	an <- sqrt(2 * log(n)) + (1/2*log(log(n)) + log(0.8197466 / 2 /sqrt(pi))) / sqrt(2 * log(n))
	
#	bn <- 1 / sqrt(2 * log(n))
	
	for (i in 1:N) {
		
#		eps <- 2 * rbinom(n, 1, 1/2) - 1

		eps <- sample(c(-1, 0, 1), n, TRUE, c(1-prob0-prob1, prob0, prob1))
				
		cur.max.use <- multi_sup_norm(eps)
		
		max.use.sample[i] <- cur.max.use
		
#		max.abs.use.sample[i] <- cur.max.use$max.abs.use
		
#		max.use.sample.scaled <- (max.use.sample - an) / bn
		
#		max.abs.use.sample.scaled <- (max.abs.use.sample - an) / bn
		
	}
	
#	list(max.use.sample = max.use.sample, max.abs.use.sample = max.abs.use.sample, max.use.sample.scaled = max.use.sample.scaled, max.abs.use.sample.scaled = max.abs.use.sample.scaled)

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
