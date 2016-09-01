###math, statistics & probability
* nonLinSolve.m - solves a set of non-linear equations (gradient free)
* isUnitCirclePolynomial.m - test if the roots of a given polynomial lie inside the unit circle
* stabilizePoly.m - reflect polynomial unstable roots inside the unit circle (make a polynomial stable)
* danChol.m - performs the cholesky factorization for positive semi definite matrices
* danQR.m - performs a generalized ("partial") QR factorization
* danDetrend.m - removes any order of trend from a given observation (make it wide sense stationary)
* polygonMean.m - calculate the mean polygon of several 2D polygons
* nJacob.m - numericaly calculate the forward / centeral jacobian of a vectorial function
             (see \systems & simulation\danLinMod.m for usage of complex step jacobian calculation)
* nGrad.m - numericaly calculate the forward / centeral gradient of a vectorial function
* chebyQuad.m - 1D function integration utilizing chebyshev polynomial approximation
* PCA.m - full principle component analysis (with appropriate Z-score, covariance eigenvalues and t-squared statistics)
* cloudAngle.m - estimate the main axis of scattered 2D data using principle component analysis
* toBand.m - expand a given matrix to band matrix (where the given matrix is its operator)
* envelope.m - retrieve the envelope of a given 2D data set
* erfir.m - calculate the error function of pure imaginary\real argument
* partialFrac.m - factorize c(x)/(a(x)*b(x)) polynomial to s(x)/a(x) + t(x)/b(x) form
* polyDiv.m - given a(x) & b(x), returns q(x) & r(x) which fulfil the equation a(x) = q(x) * b(x) + r(x)
* secondOrderNonHom.m - numerically solve X''(t) = p(x) * X'(t) + q(x) * X(t) + r(t) with given boundary conditions
* quadMin.m - solves a constrained quadratic eigenvalue problem
* confCorr.m - calculates the confidence interval of the correlation coefficient
* isGaussianLinear.m - tests a time series for gaussian (i.e - no skewness) and linear
* kolmogorovSmirnov.m - perform the two-sample kolmogorov-smirnov goodnes of fit hypothesis test
* ECDF.m - calculate a given vector empirical cumulative distribution function
* EPDF.m - estimate a given vector probability density
* jack.m - determine the bias and standard error of any statistical calculation performed on a given vector
* bootstrap.m - determine the bootstrap percentile interval of any statistical calculation performed
                on a given vector and in a given confidence interval
* mutualInformtaion.m - estimate the mutual information between two vectors
* distGen.m - generate random numbers according to a defined discrete probability distribution
* distSample.m - generate a random sample from an arbitrary discrete/finite distribution
* quat.m - returns the sampled quantiles of a given vector
* logProb.m - returns the log(P(lo < Z < hi)), where Z~N(0,1), for the region [lo, hi]
* danDct.m -  performs the direct or inverse discrete cosine transform on a given array
