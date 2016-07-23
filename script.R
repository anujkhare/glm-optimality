works_with_R("3.2.3",
             data.table="1.9.7",
             ggplot2="1.0.1",
             directlabels="2015.6.17",
             microbenchmark="1.4.2.1")

require(iregnet)
library(ggplot2)
require(survival)

## The spams package can be downloaded from
## http://spams-devel.gforge.inria.fr/hitcounter2.php?file=33815/spams-R-v2.5-svn2014-07-04.tar.gz
get_data <- function()
{
	data(prostate,package="ElemStatLearn")
	pros <- subset(prostate,select=-train,train==TRUE)
	ycol <- which(names(pros)=="lpsa")
	X.unscaled <- as.matrix(pros[-ycol])
	y.unscaled <- pros[[ycol]]
	M <- matrix(
	  colMeans(X.unscaled), nrow(X.unscaled), ncol(X.unscaled), byrow=TRUE)
	X.centered <- X.unscaled - M
	sd.vec <- apply(X.unscaled, 2, sd)
	S <- diag(1/sd.vec)
	X.scaled <- X.centered %*% S
	dimnames(X.scaled) <- dimnames(X.unscaled)
	m <- mean(y.unscaled)
	sigma <- sd(y.unscaled)
	y.scaled <- (y.unscaled - m)/sigma
	
	## X and y will be used in the various solvers.
	X <- X.scaled
	y <- y.scaled
	y_censored <- cbind(y, y)
	list("x" = X, "y" = y, "y_censored" = y_censored)
}

compute_gradient <- function(X, y_censored, beta, scale, dist) {
  stopifnot(is.numeric(X))
  stopifnot(is.numeric(y))
  stopifnot(is.matrix(X))
  stopifnot(is.matrix(y_censored))
  stopifnot(nrow(X) == nrow(y_censored))
  stopifnot(ncol(y_censored) == 2)
  stopifnot(is.numeric(beta))
  stopifnot(ncol(X)+1 == length(beta))
  stopifnot(is.numeric(scale))

  X_1 <- cbind(rep(1, nrow(X)), X)
  eta <- X_1  %*% beta
  # gives us a list ("gradients" = grad of 1/n*LL wrt (intercept, weights, scale)
  #                  "mu" = grad of LL wrt eta)
  gradients <- iregnet_compute_gradients(X_1, y_censored, eta, scale, dist)
  # we need wrt NLL (cost)
  gradients$gradients = gradients$gradients * -1
  gradients$mu = gradients$mu * -1

  return (gradients)
}

subdifferentialCriteria_iregnet <- function
### Compute subdifferential optimality criteria for the elastic net problem.
(gradient.vec,
### Numeric vector of the gradients (intercept + p weights + scale)
 gradient.intercept,
 gradient.scale,
 w,
### Numeric vector of weights (intercept + p features).
 lambda,
### Numeric lambda regularization parameter (scalar).
 alpha
### Numeric elastic net parameter (between 0 and 1).
 ){
  stopifnot(is.numeric(gradient.vec))
  stopifnot(is.numeric(w))
  stopifnot(is.numeric(lambda))
  stopifnot(length(lambda) == 1)
  stopifnot(length(gradient.vec) == length(w))
  stopifnot(is.numeric(alpha))
  stopifnot(length(alpha)==1)
  stopifnot(0 <= alpha && alpha <= 1)

  positive.part <- function(x)ifelse(x<0, 0, x)

  common.term <- gradient.vec + lambda * (1-alpha) * w
  subdiff <- ifelse(w==0,
             positive.part(abs(common.term) - lambda*alpha),
             abs(common.term + lambda*alpha*sign(w)))

  return (cbind(abs(gradient.intercept), subdiff, abs(gradient.scale)))
### p+2-vector of subdifferential optimality criteria. Lower values mean
### closer to the global optimum of the elastic net problem.
}

plot_optimality_criteria <- function(coef.mat.list, scale.mat, coef.mat.lambda, title, fname, alpha,
																	 	 X=X, y_censored=y_censored, iregnet.dist=iregnet.dist)
{
	convergence.list <- list()
	for(pkg in names(coef.mat.list)){
	  print(pkg)
	  coef.mat <- coef.mat.list[[pkg]]
		# FIXME: CHANGED
		# coef.mat.lambda <- lambda.path.list[[pkg]]
	  scale.list <- scale.mat[[pkg]]
	  for(lambda.i in seq_along(coef.mat.lambda)){
	    weight.vec <- coef.mat[lambda.i, ]
	    lambda <- coef.mat.lambda[[lambda.i]]
	    scale <- scale.list[[lambda.i]]
	
	    p <- ncol(X)
      print(p)
	    # print (scale)
	    gradients <- compute_gradient(X, y_censored, weight.vec, scale, iregnet.dist)
	    mu <- gradients$mu
	    gradients <- gradients$gradients * scale ** 2 # FIXME: scale factor in calculated values! :P
	    # gradients <- gradients$gradients
	    param.criteria <- subdifferentialCriteria_iregnet(gradients[2:(p+1)], gradients[1], gradients[p+2],
	                                                      weight.vec[2:(p+1)], lambda, alpha)
	
	    param.criteria <- param.criteria[2:(p+1)]
	    criterion.value <- c(
	      # dualityGap=lassoDualityGap(y, X, weight.vec, lambda),
	      subdifferentialL1=sum(abs(param.criteria)),
	      subdifferentialLInf=max(abs(param.criteria)),
	      subdifferentialL2=sqrt(sum(param.criteria * param.criteria)))
	
			# FIXME: CHANGED!
	    # convergence.list[[paste(pkg, lambda_path[[lambda.i]])]] <- data.table(
	    convergence.list[[paste(pkg, lambda)]] <- data.table(
	      pkg, lambda, criterion.value, criterion.name=names(criterion.value))
		}
	}
	convergence <- do.call(rbind, convergence.list)
	# print(convergence)
  print("done there!")
	
	## Accuracy of different solvers using the four different criteria.
	with.legend <- ggplot()+
	  ggtitle(title)+
	  theme_bw()+
	  theme(panel.margin=grid::unit(0, "lines"))+
	  facet_grid(criterion.name ~ ., scales="free")+
	  geom_point(aes(log10(lambda), log10(criterion.value), color=pkg),
	             shape=1,
	             data=convergence)
	(with.labels <- direct.label(with.legend+xlim(min(log10(lambda)), 0), "first.polygons"))
	pdf(fname)
	print(with.labels)
	dev.off()
}

# get_var_path <- function(fit.iregnet, lambda_path)
# {
#   iregnet.path.list <- list()
#   for(i in 1:length(lambda_path)) {
#   	iregnet.path.list[[paste(i)]] <- data.table(
#   			lambda = fit.iregnet$lambda[[i]],
#   			coef = fit.iregnet$beta,
#   			arclength = sum(abs(fit.iregnet$beta[, i])),
#   			variable = c("intercept", colnames(X))
#   		)
#   }
#   # iregnet.path <- do.call(rbind, iregnet.path.list)
# }

fit_and_plot <- function(thresholds, iregnet.dist, iregnet.alpha, X, y,
                        iregnet.scale_init, iregnet.estimate_scale, iregnet.maxiter, iregnet.unreg)
{
  thresholds <- sort(thresholds) # so that the lambda path is as accurate as possible
  
  coef.mat.list <- list()
  scale.mat <- list()
  coef.mat.lambda <- NA
  
  for (threshold in thresholds) {
    # fixing scale at 1 gives better opt. criteria than estimating scale
    name <- as.character(threshold)
    
    print(y)
    fit.iregnet <- iregnet(X, y, alpha=iregnet.alpha, iregnet.dist, maxiter=iregnet.maxiter, threshold=threshold,
                           scale_init=iregnet.scale_init, estimate_scale=iregnet.estimate_scale, unreg_sol=iregnet.unreg)
    print(fit.iregnet)
    
    # only store the lambdas for the first fit, use the same for other fits
    if (is.na(coef.mat.lambda)) {
      lambda_path <- fit.iregnet$lambda * (fit.iregnet$scale ** 2) # LAMBDA PATH TO USE
      coef.mat.lambda <- lambda_path
    }
    
    # NOTE: get_var_path has the code for the path of the vars!
    coef.mat.list[[name]] = t(fit.iregnet$beta)
    lapply(coef.mat.list, head)
    
    scale.mat[[name]] <- fit.iregnet$scale
  }

  cen <- ifelse(any(is.na(y)), "censored", "uncensored")
  title <- sprintf("%s, %s, alpha %.2f: subdifferntial optimality criteria", iregnet.dist, cen, iregnet.alpha)
  fname <- sprintf("plots/iregnet_%s_%s_%s.pdf", iregnet.dist, cen, as.character(iregnet.alpha))
  ## Compute the subdifferentials and plot the norms
  print("HEY HEY HEY!")
  plot_optimality_criteria(coef.mat.list, scale.mat, coef.mat.lambda, title, fname, iregnet.alpha,
  												 X=X, y_censored=y, iregnet.dist=iregnet.dist)
}

get_xy <- function(n_obs, n_vars, type = c('none', 'right', 'left', 'interval'), standardize=std, positive=F) {
  type <- match.arg(type)

  x <- matrix(rnorm(n_obs * n_vars), n_obs, n_vars)
  y <- rnorm(n_obs)
  y_u <- rnorm(n_obs)

  # standardize x and y
  if (standardize == T) {
    for (i in 1:ncol(x)) {
      x[, i] <- (x[, i] - mean(x[, i])) / sd(x[, i]);
    }

    y <- (y - mean(y))
    y <- y / sd(y)

    if (type == 'interval') {
      y_u <- (y_u - mean(y_u))
      y_u <-  y_u / sd(y_u)
    }
  }

  if (positive) {
    y <- abs(y)
    y_u <- abs(y_u)
  }

  if (type == "none") {
    status = rep(1, length(y))
    y_surv <- Surv(time = y, event = status, type = "right")

  } else if (type == 'interval') {
    status <- sample(c(0, 1, 2, 3), size=n_obs, replace=T)
    y <- cbind(y, y_u)
    y <- t(apply(y, 1, sort))    # make sure intervals are increasing in time
    y_surv <- Surv(time=y[, 1], time2=y[, 2], event = status, type = 'interval')

  } else {    # left or right censored
    status <- sample(c(0, 1), size=n_obs, replace=T)
    y_surv <- Surv(time = y, event = status, type = type)
  }

  # get the y matrix
  y <- cbind(y, y)
  if (type=="right") {
    y[status == 0, 2] = NA
  } else if (type == 'interval') {
    y <- NA  # NOTE: Not implemented, not needed, use the Surv object!
  } else {
    y[status == 0, 1] = NA
  }

  return (list("x" = x, "y" = y, "surv" = y_surv))
}
####################################################################################################
####################################################################################################
## Do the plotting!
xy <- get_data()
X <- xy$x
y <- xy$y
y_censored <- xy$y_censored


# iregnet fit
iregnet.scale_init <- NA
iregnet.estimate_scale <- T
iregnet.maxiter <- 1e3
iregnet.unreg <- F

thresholds <- c(1e-17, 1e-4, 1e-7)

##################### ElemStatsLearn::p (uncensored)
# fit_and_plot(thresholds, iregnet.dist = "gaussian", iregnet.alpha = 1, X, y_censored,
#              iregnet.scale_init, iregnet.estimate_scale, iregnet.maxiter, iregnet.unreg)
# fit_and_plot(thresholds, iregnet.dist = "gaussian", iregnet.alpha = 0.4, X, y_censored,
#              iregnet.scale_init, iregnet.estimate_scale, iregnet.maxiter, iregnet.unreg)
# fit_and_plot(thresholds, iregnet.dist = "logistic", iregnet.alpha = 1, X, y_censored,
#              iregnet.scale_init, iregnet.estimate_scale, iregnet.maxiter, iregnet.unreg)
# fit_and_plot(thresholds, iregnet.dist = "logistic", iregnet.alpha = 0.4, X, y_censored,
#              iregnet.scale_init, iregnet.estimate_scale, iregnet.maxiter, iregnet.unreg)

###################### Using survival::ovarian data (right censored)
data(ovarian)
X <- cbind(ovarian$ecog.ps, ovarian$rx)
y_l <- ovarian$futime
y_r <- ovarian$futime
y_r[ovarian$fustat == 0] <- NA
y_surv <- cbind(y_l, y_r)
y_surv <- abs(y_surv)

# fit_and_plot(thresholds, iregnet.dist = "gaussian", iregnet.alpha = 1, X, y_surv,
#              iregnet.scale_init, iregnet.estimate_scale, iregnet.maxiter, iregnet.unreg)
# fit_and_plot(thresholds, iregnet.dist = "gaussian", iregnet.alpha = 0.4, X, y_surv,
#              iregnet.scale_init, iregnet.estimate_scale, iregnet.maxiter, iregnet.unreg)
# fit_and_plot(thresholds, iregnet.dist = "logistic", iregnet.alpha = 1, X, y_surv,
#              iregnet.scale_init, iregnet.estimate_scale, iregnet.maxiter, iregnet.unreg)
# fit_and_plot(thresholds, iregnet.dist = "logistic", iregnet.alpha = 0.4, X, y_surv,
#              iregnet.scale_init, iregnet.estimate_scale, iregnet.maxiter, iregnet.unreg)

y_surv <- log(y_surv)
# fit_and_plot(thresholds, iregnet.dist = "gaussian", iregnet.alpha = 1, X, y_surv,
#              iregnet.scale_init, iregnet.estimate_scale, iregnet.maxiter, iregnet.unreg)
# fit_and_plot(thresholds, iregnet.dist = "gaussian", iregnet.alpha = 0.4, X, y_surv,
#              iregnet.scale_init, iregnet.estimate_scale, iregnet.maxiter, iregnet.unreg)
fit_and_plot(thresholds, iregnet.dist = "logistic", iregnet.alpha = 1, X, y_surv,
             iregnet.scale_init, iregnet.estimate_scale, iregnet.maxiter, iregnet.unreg)
fit_and_plot(thresholds, iregnet.dist = "logistic", iregnet.alpha = 0.1, X, y_surv,
             iregnet.scale_init, iregnet.estimate_scale, iregnet.maxiter, iregnet.unreg)
