mean.from.cpt <- function(x, cpt) {



	n <- length(x)

	len.cpt <- length(cpt)

	if (len.cpt) cpt <- sort(cpt)

	beg <- endd <- rep(0, len.cpt+1)

	beg[1] <- 1

	endd[len.cpt+1] <- n

	if (len.cpt) {

		beg[2:(len.cpt+1)] <- cpt+1

		endd[1:len.cpt] <- cpt

	}

	means <- rep(0, len.cpt+1)

	for (i in 1:(len.cpt+1)) means[i] <- mean(x[beg[i]:endd[i]])

	rep(means, endd-beg+1)

}



universal.M.th.v3 <- function(n, lambda = 0.9) {

		

	mat.90 <- matrix(0, 24, 3)

	mat.90[,1] <- c(10, 50, 100, 150, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1500, 2000, 2500, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000)

	mat.90[,2] <- c(1.420, 1.310, 1.280, 1.270, 1.250, 1.220, 1.205, 1.205, 1.200, 1.200, 1.200, 1.185, 1.185, 1.170, 1.170, 1.160, 1.150, 1.150, 1.150, 1.150, 1.145, 1.145, 1.135, 1.135)

	mat.90[,3] <- rep(100, 24)

	

	mat.95 <- matrix(0, 24, 3)

	mat.95[,1] <- c(10, 50, 100, 150, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1500, 2000, 2500, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000)

	mat.95[,2] <- c(1.550, 1.370, 1.340, 1.320, 1.300, 1.290, 1.265, 1.265, 1.247, 1.247, 1.247, 1.225, 1.225, 1.220, 1.210, 1.190, 1.190, 1.190, 1.190, 1.190, 1.190, 1.180, 1.170, 1.170)

	mat.95[,3] <- rep(100, 24)



	if (lambda == 0.9) A <- mat.90 else A <- mat.95



	d <- dim(A)

	if (n < A[1,1]) {

		th <- A[1,2]

		M <- A[1,3]

	}

	else if (n > A[d[1],1]) {

		th <- A[d[1],2]

		M <- A[d[1],3]

	}

	else {

		ind <- order(abs(n - A[,1]))[1:2]

		s <- min(ind)

		e <- max(ind)

		th <- A[s,2] * (A[e,1] - n)/(A[e,1] - A[s,1]) + A[e,2] * (n - A[s,1])/(A[e,1] - A[s,1])

		M <- A[s,3] * (A[e,1] - n)/(A[e,1] - A[s,1]) + A[e,3] * (n - A[s,1])/(A[e,1] - A[s,1])

	}



	list(th.const=th, M=M)

}


all.intervals.flat <- function(n) {
	
	if (n == 2) ind <- matrix(1:2, 2, 1) else {
		M <- (n-1)*n/2	
		ind <- matrix(0, 2, M)
		ind[1,] <- rep(1:(n-1), (n-1):1)
		ind[2,] <- 2:(M+1) - rep(cumsum(c(0, (n-2):1)), (n-1):1)
	}
	ind

}

grid.intervals <- function(n, M) {
		
	if (n==2) ind <- matrix(c(1, 2), 2, 1)
	
	else if (M >= (n-1)*n/2) ind <- all.intervals.flat(n)
	
	else {
		k <- 1
		while (k*(k-1)/2 < M) k <- k+1
		ind2 <- all.intervals.flat(k)
		ind2.mx <- max(ind2)
		ind <- round((ind2 - 1) * ((n-1) / (ind2.mx-1)) + 1)
	}	
	
	ind	
}
