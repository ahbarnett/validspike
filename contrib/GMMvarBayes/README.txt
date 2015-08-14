Description 	

This is the variational Bayesian procedure (also called mean field) for inference of Gaussian mixture model. This is the Bayesian treatment of Gaussian mixture model.

Unlike the EM algorithm (Maximum likelihood estimation), it can automatically determine the number of the mixture components k.

Example code:
load data;
label=vbgm(x,10);
spread(x,label);

The data set is of 3 Gaussian. You only need set a number (say 10) which is larger than the intrinsic number of components. The algorithm will automatically find the right k.

Detail description of the algorithm can be found in the reference.

Reference: Pattern Recognition and Machine Learning by Christopher M. Bishop (P.474)
Acknowledgements 	

This file inspired Em Algorithm For Gaussian Mixture Model.
MATLAB release 	MATLAB 7.13 (R2011b)

