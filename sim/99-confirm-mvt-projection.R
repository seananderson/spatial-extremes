library(cluster)
library(mvtnorm)
library(MASS)

# Simulate data and grid
grid = as.matrix(expand.grid("lon" = seq(5,15,1), "lat" = seq(5,15,1)))
nLocs = dim(grid)[1]
nKnots = 20
knots = jitter(pam(grid,nKnots)$medoids)
distKnots = as.matrix(dist(knots))
distKnotsSq = distKnots^2 # squared distances
# note: shape parameter scaled to distance matrix
gp_theta = 0.001
sigma.norm = 0.001
corKnots = exp(-gp_theta*distKnotsSq)
Sigma.normal = corKnots * sigma.norm * sigma.norm

# Plot distance vs covariance
plot(distKnots, Sigma.normal, xlab="Distance between knots", ylab="Covariance")

# We'll replicate the knots a 2nd time, and use the projection equation to project to them.
# This is a good comparison because the relative distances between them stays the same.
# Calculate distance from knots to grid
distAll = as.matrix(dist(rbind(knots, knots)))^2
distKnots21Sq = t(distAll[1:nKnots, -c(1:nKnots)])
Sigma21.normal = exp(-gp_theta*distKnots21Sq) * sigma.norm * sigma.norm

re.norm = rmvt(1, sigma = Sigma.normal, df = 6)
re.norm = re.norm-mean(re.norm)# Scale

# Generate a lot of random samples, and show that empirically,
# we can recover the covariance
sim.norm = rmvt(10000, sigma = Sigma.normal, df = 6)
estCov = cov.trob(sim.norm, cor = FALSE, center = TRUE, nu = 6)

plot(Sigma.normal, estCov$cov, xlab="True values", ylab = "Estimated")
abline(0,1,col="red")

# take inverse for projection
invSigmaKnots.norm = solve(Sigma.normal)

# Project
proj.norm = t((Sigma21.normal %*% invSigmaKnots.norm) %*% t(sim.norm))

# Try to estimate covariance with cov.trob
sim.norm.proj = cov.trob(proj.norm, cor = FALSE, center = TRUE, nu = 6)
points(Sigma.normal, sim.norm.proj$cov, col="blue",cex=0.5)

