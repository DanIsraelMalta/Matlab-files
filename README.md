#Matlab files backup repository

###estimation, processing & identification
* TFidentify.m - identify a transfer function of a process (frequency domain method)
* ARidentify.m - identify an autoregressive process (time domain method)
* MAidentify.m - identify a moving-average process (time domain method)
* multiARestimate.m - perform a multivariate adaptive autoregressive process estimation
                     (time domain method; output is in kalman filter form)
* markovIdentification.m - given static model and an output signal, estimate model parameters
                          (markov chain monte carlo method)
* yuleWalkerAR.m - identify an autoregressive process (time domain method) using Yule-walker method.
* OptimalFilter.m - designs an optimal direct form low pass filter
* sTrans.m - implements stockwell frequency transform
* wpsd.m - estimate a signal power spectral density using Welch's method
* cpsd.m - estimate two signals cross power spectral density, coherence & phase
* modelOrder.m - estimate a given signal subspace order
* adaptWiener.m - implements an adaptive wiener filter
* burg.m - implements burg's method of linear prediction
* autocorrExt.m - extends an autocorrelation sequence using levinson-durbin recrusion
* firWiener.m - design a finite impulse response wiener filter using levinson-durbin
* levinsonDurbin.m - performs the levinson-durbin recrusion
* kmeans.m - perform the k-means clustering algorithm (using euclidean distance)
* chirpZ.m - perform the chirp z trnasform in a given frequency span
* danGaussfir.m - returns the coeeficients of gaussian filter (good for pulse shaping)
* stepFinder.m - extracts a piecewise constant signal from a noisy observation
                 using a maximum likelihood based penalization schema

###smoothers
* bayesianSmooth.m - smooths an observation by minimizing a quadratic cost function using the bayesian approach
* directSmooth.m - smooths an observation by minimizing a quadratic cost using the maximum likelihood approach
* fftSmooth.m - smooths an observation by minimizing the general cross validation function, but in the frequency domain!
* totalVarFilter.m - smooths an observation using the total-variance schema
* termSmooth.m - smooths an observation such that its output will either have a non-negative impulse
                 response or a perfect response to gradual curvature
* impulseSmooth.m - spike smoothing using first & second order moments
* RTSsmooth.m - rauch-tung-striebel smoother
* estimateVar.m - estimate the variance of an additive noise to a given signal (then a perfect smoother can be created)
* splineSmooth.m - smooths an observation using a spline
* danSgolay.m - fast savitsky-golay smoothing for a uniformely spaced data set
* bandPass.m - perfect band pass filter
* loess2D.m - 2D LOESS method smoothing

###systems & simulation
* ODE87.m - 8'th order dorman prince explicit integrator
* SDF.m - vairable-mass rigid-body six-degrees-of-freedom equation motion solver
* danLinMod.m - linmod for any sort of non-simulink non-linear dynamic system
* modelReduction.m - reduces system order by removing fast dynamics (not pole-zero overlap!) while maintaining static gain
* danDLQR.m - enhanced (more features) version of matlab's control toolbox function 'dlqr'
* lowFreqPart.m - express a given transfer function as the sum of its low frequency dominant component
                  and an additive term, both in the form of partial fraction expression
* danImpulse.m - matlab's control toolbox 'impulse' function for non LTI objects

###curve fitting & regression
* optimalPoly.m - automatically find's a polynomial optimal degree (for polynomial fitting; maximum likelihood style)
* orthogonalPolyFit.m - fit's an orthogonal polynomial (not linear least square polynomial)
* loessFit.m - LOESS method regression
* scatterReduction.m - scatterplot interpolant using grouped first order difference smoothing
* danPchip.m - enhanced version of matlab's builtin function pchip.

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

###general
* runFunc.m - implements a general moving window inlined scalar function upon a given vector
* corrAlign.m - align two numeric vectors such that their inner product (correlation) is maximal
* linDependent.m - finds linearly dependent columns in a given matrix
* danBuffer.m - partitions a vector into matrix whos columns are non-overlapping data segments ("frames") of a given length
               (zero padded if necessary)
* logSpan.m - logarithmically grid's a given 1D region with required number of points
* space.m - spacially distibuted grid of a given 1D region, relative to its edges or center.
* oneLinersLib.m - a library of usefull function handlers
                       
###dedicated graphical functions:
* qqplot.m - draws a quantile-quantile plot (theoretical quantiles are plotted against sample order statistics)
* pvd.m - progressive vector plot (using quiver)
* modelCompare.m - a comparison diagran of various test models (compared against a given reference model)
