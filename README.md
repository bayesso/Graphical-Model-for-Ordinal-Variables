# Graphical-Model-for-Ordinal-Variables
This R function implements the algorithms proposed in paper "Graphical model for ordinal variables"

	## Input:
	## --X: a matrix recording n observations (n rows), each with p features (p columns)
  ## --lambda: a scalar tuning parameter balancing the likelihood function and the penalty function
  ## Output: a list object with the following components
  ## --Omega_value: sample inverse covarince matrix (p by p matrix) calculated from input data matrix X
  ## --Sigma_value: sample covarince matrix (p by p matrix) calculated from input data matrix X
  ## --Theta_value: estimated inverse covariance matrix (diagnal scaled to 1) of the latent Gaussian graphical model
  
  
