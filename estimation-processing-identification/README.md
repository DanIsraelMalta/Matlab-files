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
