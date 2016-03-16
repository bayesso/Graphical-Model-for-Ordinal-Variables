
## Remove the junk variables in the memory
rm(list=ls(all=TRUE))



truncate_component <- function(Y, theta=NA, mu=0, sigma=1)
{
  	delta1 <- (theta[Y] - mu) / sigma
  	delta2 <- (theta[Y+1] - mu) / sigma
  	tmp1 <- (dnorm(delta1) - dnorm(delta2)) / (pnorm(delta2) - pnorm(delta1))
  	EX <- mu + tmp1 * sigma
  	
  	delta1[delta1 < -1e+10] <- -1e+10
  	delta2[delta2 > 1e+10] <- 1e+10
  	tmp2 <- (delta1*dnorm(delta1) - delta2*dnorm(delta2)) / (pnorm(delta2) - pnorm(delta1))
  	EXX <- sigma^2 + mu^2 + sigma^2 * tmp2 + 2*mu*sigma*tmp1
  
	return(list(EX=EX, EXX=EXX))
}


impute_S <- function(X, Z, Theta_value, Sigma_value)
{            
	p <- ncol(X)
	n <- nrow(X) 
	Z_new <- matrix(0, n, p)
	diag_element <- rep(0, p)
  
	for (j in seq(1, p))
	{
		## mean and variance for conditional Gaussian
		tmp <- matrix(Sigma_value[j, -j], 1, p-1)
		tmp1 <- solve(Sigma_value[-j, -j])

		mu <- tmp %*% tmp1 %*% t(Z[, -j])				
		mu <- as.vector(mu)		 
		sigma <- Sigma_value[j, j] - tmp %*% tmp1 %*% t(tmp)
		sigma <- sqrt(sigma)  		
      
		obj <- truncate_component(Y=X[, j], theta=Theta_value[[j]]$cutoff, mu=mu, sigma=sigma)      

		Z_new[, j] <- obj$EX
		diag_element[j] <- mean(obj$EXX)
	}
	ES <- t(Z_new) %*% Z_new / n
	diag(ES) <- diag_element

	output <- list()
	output$ES <- ES
	output$Z <- Z_new
	return(output)         
}


rating_graphical_model <- function(X, lambda_value)
{
	## Input:
	## --X: a matrix recording n observations (n rows), each with p features (p columns)
    ## --lambda: a scalar tuning parameter balancing the likelihood function and the penalty function
    ## Output: a list object with the following components
    ## --Omega_value: sample inverse covarince matrix (p by p matrix) calculated from input data matrix X
    ## --Sigma_value: sample covarince matrix (p by p matrix) calculated from input data matrix X
    ## --Theta_value: estimated inverse covariance matrix (diagnal scaled to 1) of the latent Gaussian graphical model

	## This function requires the R package "glasso"
	require(glasso)

    p <- ncol(X)
    n <- nrow(X)
    
    num_iter <- 0
    max_iter <- 25
    tol_value <- 1e-3
    diff_value <- 1e+10
    
  	## Estimate theta
	Theta_value <- vector('list', length=p)
    for (j in seq(1, p))
	{
        K <- length(unique(as.vector(X[, j])))
        tmpvec <- rep(0, K+1)
        tmpvec[1] <- -Inf
        tmpvec[K+1] <- Inf
        for (k in seq(2, K))
  	    {
  		    tmp = mean(X[, j] < k)
  		    tmpvec[k] <- qnorm(tmp)
  	    }
        Theta_value[[j]]$cutoff <- tmpvec
        #Theta_value[[j]]$ifordinal <- T
    }
	
    
    ## Initialize S
	Z <- matrix(0, n, p)
	diag_element <- rep(0, p)
	for (j in seq(1, p))
	{
	    tmp2 <- truncate_component(Y=X[, j], theta=Theta_value[[j]]$cutoff)

	    Z[, j] <- tmp2$EX
		diag_element[j] <- mean(tmp2$EXX)
	}
	ES <- t(Z) %*% Z / n
	diag(ES) <- diag_element

	
    ## Initalize Omega
	obj <- glasso(s=ES, rho=lambda_value, maxit=1000, penalize.diagonal=F)
	Omega_value <- (t(obj$wi) + obj$wi) / 2
    Sigma_value <- (t(obj$w) + obj$w) / 2
	sd_marginal <- sqrt(diag(Sigma_value))
	sd_marginal[abs(sd_marginal) < 1e-10] <- 1e-10
	Sigma_value <- diag(1/sd_marginal) %*% Sigma_value %*% diag(1/sd_marginal)
	Omega_value <- diag(sd_marginal) %*% Omega_value %*% diag(sd_marginal)

    ## EM iteration
    while((num_iter < max_iter) & (diff_value > tol_value))
    {
		## E-step: impute Z
		obj_Estep <- impute_S(X=X, Z=Z, Theta_value=Theta_value, Sigma_value=Sigma_value)

		ES <- obj_Estep$ES
		Z <- obj_Estep$Z
		
		## M-step: estimate Omega
		obj <- glasso(s=ES, rho=lambda_value, maxit=1000, penalize.diagonal=F)
		Omega_new <- (t(obj$wi) + obj$wi) / 2
		Sigma_value <- (t(obj$w) + obj$w) / 2

		sd_marginal <- sqrt(diag(Sigma_value))
		sd_marginal[abs(sd_marginal) < 1e-10] <- 1e-10
		Sigma_value <- diag(1/sd_marginal) %*% Sigma_value %*% diag(1/sd_marginal)
		Omega_new <- diag(sd_marginal) %*% Omega_new %*% diag(sd_marginal)

		## Check stopping criterion
		diff_value <- sum(abs(Omega_new - Omega_value))
		num_iter <- num_iter + 1
		Omega_value <- Omega_new
		# cat(sprintf('iter=%d,    diff=%f\n', num_iter, diff_value))
	}

	output <- list()
	output$Omega_value <- Omega_value
	output$Sigma_value <- Sigma_value
	# output$Z <- Z
	# output$ES <- ES
	# output$num_iter <- num_iter
	output$Theta_value <- Theta_value

	return(output)
}





