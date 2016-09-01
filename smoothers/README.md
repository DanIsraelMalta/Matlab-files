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
