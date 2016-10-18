###systems & simulation
* ODE87.m - 8'th order dorman prince explicit integrator
* SDF.m - vairable-mass rigid-body six-degrees-of-freedom equation motion solver
* danLinMod.m - linmod for any sort of non-simulink non-linear dynamic system
* modelReduction.m - reduces system order by removing fast dynamics (not pole-zero overlap!) while maintaining static gain
* danDLQR.m - enhanced (more features) version of matlab's control toolbox function 'dlqr'
* lowFreqPart.m - express a given transfer function as the sum of its low frequency dominant component
                  and an additive term, both in the form of partial fraction expression
* danImpulse.m - matlab's control toolbox 'impulse' function for non LTI objects