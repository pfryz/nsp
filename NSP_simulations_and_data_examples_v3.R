smuce_coverage <- function(truth, est) {
	
	# Returns "TRUE" if the SMUCE estimator "est" is such that each confidence interval covers a true change-point of "truth", and "FALSE" otherwise. See "sim_coverage" for an example of use.
	
	res <- TRUE
	
	j <- jumpint(est)
	d <- dim(j)

	cpts <- which(abs(diff(truth)) > 0)
		
	if (d[1]-1) for (i in 1:(d[1]-1)) {
		
		cur_int <- j[i,3]:j[i,4]
		res <- res && any(cpts %in% cur_int)
	}
	
	res	
	
}


nsp_coverage <- function(truth, est) {

	# Returns "TRUE" if the NSP estimator "est" is such that each interval of significance covers a true change-point of "truth", and "FALSE" otherwise. See "sim_coverage" for an example of use.
	
	res <- TRUE
	
	d <- dim(est$intervals)

	cpts <- which(abs(diff(truth)) > 0)
	
	if (d[1]) for (i in 1:(d[1])) {
		
		cur_int <- est$intervals[i,1]:(est$intervals[i,2]-1)
		res <- res && any(cpts %in% cur_int)
	}
	
	res	
	
}


nsp_coverage_extended <- function(truth, est) {
	
	# Returns as in "nsp_coverage", plus the number of intervals of significance in "est".	
	
	res <- TRUE
	
	d <- dim(est$intervals)

	cpts <- which(abs(diff(truth)) > 0)
	
	if (d[1]) for (i in 1:(d[1])) {
		
		cur_int <- est$intervals[i,1]:(est$intervals[i,2]-1)
		res <- res && any(cpts %in% cur_int)
	}
	
	list(coverage = res, how.many = d[1])
	
}



sim_coverage <- function(truth, sigma, M = 100, N = 100) {
	
	# Simulates N realisations of truth + sigma * rnorm(length(truth)) and checks SMUCE and NSP coverage for each sample path, via "nsp_coverage" and "smuce_coverage".
	
	res.nsp <- res.smuce <- rep(0, N)
	
	for (i in 1:N) {
		print(i)
		x <- truth + sigma * rnorm(length(truth))

		x.n <- nsp_poly(x, M = M)
		print(x.n)
		res.nsp[i] <- nsp_coverage(truth, x.n)
		
		x.s <- stepFit(x, alpha = 0.1, confband = T)
		print(jumpint(x.s))
		res.smuce[i] <- smuce_coverage(truth, x.s)
		
		print(paste("NSP:", res.nsp[i]))
		print(paste("SMUCE:", res.smuce[i]))
				
	}
	
	list(res.nsp=res.nsp, res.smuce=res.smuce)
	
}

sim_coverage_extended <- function(truth, ar.coef, div.factor, M = 100, N = 100) {

	# As in "sim_coverage", but only for NSP and also returns the numbers of intervals of significance in each sample path.
	
	coverage <- how.many <- rep(0, N)
	
	for (i in 1:N) {
		print(i)
		x <- truth + arima.sim(list(ar=ar.coef), n=length(truth)) / div.factor
		
		x.n <- nsp_poly_ar(x, 1, M = M)
		cvg <- nsp_coverage_extended(truth, x.n)
		
		coverage[i] <- cvg$coverage
		how.many[i] <- cvg$how.many
				
	}
	
	list(coverage = coverage, how.many = how.many)
	
}




#####################################
### Data examples
#####################################

# Notes of caution:

# (1) Examples 1 and 2 in the paper used a version of NSP which used a random interval sampling mechanism. Code for this implementation is not and will not be provided - we have replaced it with
# a deterministic interval sampling mechanism, for improved reproducibility of results (with the random sampling mechanism reproducibility is still possible via set.seed, but this is frequently
# forgotten and so not very reliable). 
# The results for Examples 1 and 2 reported in the paper may therefore differ *slightly* from the below.

# (2) Since releasing the September 2020 version of the paper, we also significantly accelerated the NSP algorithms. The key component of the acceleration is that interval sampling is now performed
# from the shortest to the longest interval. This means the results of the scripts below, even notwithstanding issue (1) above, may differ *slightly* from those reported in the paper.



# Example 1 - The US ex-post real interest rate.
# real_dat_sc is the scaled time series, denoted in the paper (Section 5.1) by \tilde{Y}_t
# real_dat is the original time series Y_t

sig_const_orig <- mean.from.cpt(real_dat, c(47, 82))
ts.plot(real_dat, xlab="Time (quarters)", ylab="")
lines(sig_const_orig, col="red", lwd = 3)

#set.seed(1)
dd <- nsp_poly(real_dat)
draw_rects(dd, c(-10, 15), col="grey", density = 20)

sig_const <- mean.from.cpt(real_dat_sc, c(47, 82))
ts.plot(real_dat_sc, xlab="Time (quarters)", ylab="")
lines(sig_const, col="red", lwd = 3)

#set.seed(1)
ddd <- nsp_poly(real_dat_sc)
draw_rects(ddd, c(-3, 5), col="grey", density = 20)


#set.seed(1)
nsp_poly(real_dat_sc, deg=1)
sig_lin <- rep(0, 103)
real_dat_sc_11 <- real_dat_sc[1:73]
time_1 <- 1:73
sig_1 <- lm(real_dat_sc_11 ~ time_1)$fitted
sig_lin[1:73] <- sig_1

mu <- sig_1[73]
real_dat_sc_12 <- real_dat_sc[74:103]
which.max(abs(inner.prod.iter(real_dat_sc_12))) -> d
sig_2 <- mean.from.cpt(real_dat_sc_12, d)

sig_lin[74:103] <- sig_2

ts.plot(sig_lin)


real_dat_sc_1alt <- real_dat_sc[1:76]
time_1alt <- 1:76
sig_1alt <- lm(real_dat_sc_1alt ~ time_1alt)$fitted
sig_2alt <- rep(mean(real_dat_sc[77:103]), 103-77+1)
sig_linalt <- c(sig_1alt, sig_2alt)

# BIC for linear+const
103/2 * log(var(real_dat_sc - sig_linalt)) + 4/2 * log(103)
# BIC for piecewise-constant
103/2 * log(var(real_dat_sc - sig_const)) + 5/2 * log(103)

ts.plot(real_dat_sc, xlab="Time (quarters)", ylab="")
lines(sig_linalt, col="red", lwd = 3)


sig1orig <- real_dat[1:76]
tim <- 1:76
tim2 <- tim^2 / 76
dd <- lm(sig1orig ~ tim + tim2)$fitted
sig_quad <- c(dd, rep(mean(real_dat[77:103]), 103-77+1))
ts.plot(real_dat, xlab="Time (quarters)", ylab="")
lines(sig_quad, col="red", lwd = 3)




# Example 2 - UK Covid-19 deaths time series.


covid.new <- read.csv("data_2021-Jan-31.csv")
Dt <- rev(covid.new[[5]])
At <- sqrt(2 * Dt + 3/8)
ts.plot(At, xlab="Time (days starting from 29th February 2020)", ylab = "")
nsp_poly(At, deg = 1) -> At.n
draw_rects(At.n, c(0, 51))





# Example 3 - House prices in the London Borough of Newham



nm <- read.csv("n.csv")
nm[[7]] -> nmp
n <- 131

covv <- matrix(1, 130, 2)
covv[,2] <- log(nmp[1:130])
resp <- log(nmp[2:131])
nsp_tvreg(resp, covv, 1000, overlap=TRUE)

summary(lm(resp[1:60]~covv[1:60,]-1))
summary(lm(resp[61:130]~covv[61:130,]-1))

ts.plot(log(nmp), ylab="", xlab="Time (months) starting January 2010")
abline(v = 60, col="red")







# Simulated examples
# NOTE: the below is not "completely complete" - some of the quantities below are missing, but they are hopefully described precisely enough in the paper.
# If in doubt, please ask: p.fryzlewicz@lse.ac.uk .

set.seed(1)
x <- blocks + 10 * rnorm(2048)
x.n <- nsp_poly(x)
ts.plot(x, col="grey", ylab="")
draw_rects(x.n, c(-40, 50), 20, "red")
truecpt <- which(abs(diff(blocks)) > 0)
abline(v = truecpt, col="blue")

nsp_poly(x, overlap=TRUE) -> x.n.o
ts.plot(x, col="grey", ylab="")
draw_rects(x.n.o, c(-40, 50), 20, "red")
truecpt <- which(abs(diff(blocks)) > 0)
abline(v = truecpt, col="blue")

plot.ts(837:1303, x.sect, type="lines", xlab="Time", ylab="")



ts.plot(x.easy, ylab="")
cpt_importance(x.easy.n)
cpt_importance(x.easy.n.1stage)



set.seed(1)
sim.cov.res.100 <- sim_coverage(blocks, 10, 100, 100)

set.seed(1)
sim.cov.res.1000 <- sim_coverage(blocks, 10, 1000, 100)


xw <- wave2sect + rnorm(450)/2
ts.plot(xw, col="grey", ylab="")
lines(wave2sect)
draw_rects(xw.n.0, c(-1, 5), 20, "brown")


ts.plot(xw, col="grey", ylab="")
lines(wave2sect)
draw_rects(xw.n.1, c(-1, 5), 20, "red")

ts.plot(xw, col="grey", ylab="")
lines(wave2sect)
draw_rects(xw.n.2, c(-1, 5), 20, "brown")


set.seed(1)
x.rt.hard <- squarewave + rt(800, 4) * seq(from = 2, to = 8, length = 800)
x.rt.hard.sn <- nsp_poly_selfnorm(x.rt.hard)
ts.plot(x.rt.hard, ylab="")
draw_rects(x.rt.hard.sn, c(-70, 60), 20, "red")
abline(v = c(200, 400, 600), col="blue")

cvg.holger.1 <- sim_coverage_extended(sig.holger.1, 0.9, 5)

sig.holger.1 <- c(rep(0, 100), rep(1, 200), rep(0, 200), rep(2, 50), rep(0, 200), rep(-1, 250))
noise <- arima.sim(list(ar = 0.9), n = 2048)
x.dep <- sig.holger.1 + noise[1:1000]/5
nsp_poly_ar(x.dep, M=100) -> x.dep.n
