###general mathematical & statistical functions:
###math, statistics & probability
* nonLinSolve.m - solves a set of non-linear equations (gradient free)
* isUnitCirclePolynomial.m - test if the roots of a given polynomial lie inside the unit circle
\math, statistics & probability\stabilizePoly.m - reflect polynomial unstable roots inside the unit circle (make a polynomial stable)
\math, statistics & probability\danChol.m - performs the cholesky factorization for positive semi definite matrices
\math, statistics & probability\danQR.m - performs a generalized ("partial") QR factorization
\math, statistics & probability\danDetrend.m - removes any order of trend from a given observation (make it wide sense stationary)
\math, statistics & probability\polygonMean.m - calculate the mean polygon of several 2D polygons
\math, statistics & probability\nJacob.m - numericaly calculate the forward / centeral jacobian of a vectorial function
                                (see \systems & simulation\danLinMod.m for usage of complex step jacobian calculation)
\math, statistics & probability\nGrad.m - numericaly calculate the forward / centeral gradient of a vectorial function
\math, statistics & probability\cloudAngle.m - estimate the main axis of scattered 2D data using principle component analysis
\math, statistics & probability\toBand.m - expand a given matrix to band matrix (where the given matrix is its operator)
\math, statistics & probability\envelope.m - retrieve the envelope of a given 2D data set
\math, statistics & probability\erfir.m - calculate the error function of pure imaginary\real argument
\math, statistics & probability\partialFrac.m - factorize c(x)/(a(x)*b(x)) polynomial to s(x)/a(x) + t(x)/b(x) form
\math, statistics & probability\polyDiv.m - given a(x) & b(x), returns q(x) & r(x) which fulfil the equation a(x) = q(x) * b(x) + r(x)
\math, statistics & probability\confCorr.m - calculates the confidence interval of the correlation coefficient
\math, statistics & probability\isGaussianLinear.m - tests a time series for gaussian (i.e - no skewness) and linear
\math, statistics & probability\kolmogorovSmirnov.m - perform the two-sample kolmogorov-smirnov goodnes of fit hypothesis test
\math, statistics & probability\ECDF.m - calculate a given vector empirical cumulative distribution function
\math, statistics & probability\EPDF.m - estimate a given vector probability density
\math, statistics & probability\jack.m - determine the bias and standard error of any statistical calculation performed on a given vector
\math, statistics & probability\bootstrap.m - determine the bootstrap percentile interval of any statistical calculation performed
                                              on a given vector and in a given confidence interval
\math, statistics & probability\distGen.m - generate random numbers according to a defined discrete probability distribution
\math, statistics & probability\distSample.m - generate a random sample from an arbitrary discrete/finite distribution
\math, statistics & probability\quat.m - returns the sampled quantiles of a given vector

#functions related to system identification, signal processing & parameter estimation:
-------------------------------------------------------------------------------------
\estimation, processing & identification\TFidentify.m - identify a transfer function of a process (frequency domain method)
\estimation, processing & identification\ARidentify.m - identify an autoregressive process (time domain method)
\estimation, processing & identification\MAidentify.m - identify a moving-average process (time domain method)
\estimation, processing & identification\multiARestimate.m - perform a multivariate adaptive autoregressive process estimation
                                                             (time domain method; output is in kalman filter form)
\estimation, processing & identification\markovIdentification.m - given static model and an output signal, estimate model parameters
                                                                  (markov chain monte carlo method)
\estimation, processing & identification\OptimalFilter.m - designs an optimal direct form low pass filter
\estimation, processing & identification\sTrans.m - implements stockwell frequency transform
\estimation, processing & identification\wpsd.m - estimate a signal power spectral density using Welch's method
\estimation, processing & identification\cpsd.m - estimate two signals cross power spectral density, coherence & phase
\estimation, processing & identification\modelOrder.m - estimate a given signal subspace order
\estimation, processing & identification\adaptWiener.m - implements an adaptive wiener filter
\estimation, processing & identification\burg.m - implements burg's method of linear prediction
\estimation, processing & identification\autocorrExt.m - extends an autocorrelation sequence using levinson-durbin recrusion
\estimation, processing & identification\firWiener.m - design a finite impulse response wiener filter using levinson-durbin
\estimation, processing & identification\levinsonDurbin.m - performs the levinson-durbin recrusion

#smoothers desighed especially for certain cases or to maintain certain properties:
----------------------------------------------------------------------------------
\smoothers (specially designed)\bayesianSmooth.m - smooths an observation by minimizing a quadratic cost function using the bayesian approach
\smoothers (specially designed)\directSmooth.m - smooths an observation by minimizing a quadratic cost using the maximum likelihood approach
\smoothers (specially designed)\fftSmooth.m - smooths an observation by minimizing the general cross validation function, but in the frequency domain!
\smoothers (specially designed)\totalVarFilter.m - smooths an observation using the total-variance schema
\smoothers (specially designed)\termSmooth.m - smooths an observation such that its output will either have a non-negative impulse
                                               response or a perfect response to gradual curvature
\smoothers (specially designed)\impulseSmooth.m - spike smoothing using first & second order moments
\smoothers (specially designed)\RTSsmooth.m - rauch-tung-striebel smoother
\smoothers (specially designed)\estimateVar.m - estimate the variance of an additive noise to a given signal (then a perfect smoother can be created)
\smoothers (specially designed)\splineSmooth.m - smooths an observation using a spline
\smoothers (specially designed)\bandPass.m - perfect band pass filter
\smoothers (specially designed)\loess2D.m - 2D LOESS method smoothing

#functions related to system design & simulation:
------------------------------------------------
\systems & simulation\ODE87.m - 8'th order dorman prince explicit integrator
\systems & simulation\SDF.m - vairable-mass rigid-body six-degrees-of-freedom equation motion solver
\systems & simulation\danLinMod.m - linmod for any sort of non-simulink non-linear dynamic system
\systems & simulation\modelReduction.m - reduces system order by removing fast dynamics (not pole-zero overlap!) while maintaining static gain
\systems & simulation\danDLQR.m - enhanced (more features) version of matlab's control toolbox function 'dlqr'
\systems & simulation\lowFreqPart.m - express a given transfer function as the sum of its low frequency dominant component
                                      and an additive term, both in the form of partial fraction expression

#general functions:
------------------
\general\runFunc.m - implements a general moving window inlined scalar function upon a given vector
\general\corrAlign.m - align two numeric vectors such that their inner product (correlation) is maximal
\general\linDependent.m - finds linearly dependent columns in a given matrix
\general\danBuffer.m - partitions a vector into matrix whos columns are non-overlapping data segments ("frames") of a given length
                       (zero padded if necessary)
                       
#curve fitting & regression functions:
-------------------------------------
\curve fitting & regression\optimalPoly.m - automatically find's a polynomial optimal degree (for polynomial fitting; maximum likelihood style)
\curve fitting & regression\orthogonalPolyFit.m - fit's an orthogonal polynomial (not linear least square polynomial)
\curve fitting & regression\loessFit.m - LOESS method regression

#dedicated graphical functions:
------------------------------
\visualization\qqplot.m - draws a quantile-quantile plot (theoretical quantiles are plotted against sample order statistics)
\visualization\pvd.m - progressive vector plot (using quiver)                       
