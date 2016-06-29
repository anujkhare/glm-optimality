works_with_R("3.2.3",
             data.table="1.9.7",
             ggplot2="1.0.1",
             directlabels="2015.6.17",
             microbenchmark="1.4.2.1",
             penalized="0.9.42",
             glmnet="1.9.5",
             spams="2.5",
             lars="1.2")

## The spams package can be downloaded from
## http://spams-devel.gforge.inria.fr/hitcounter2.php?file=33815/spams-R-v2.5-svn2014-07-04.tar.gz
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

fit.lars <- lars(X, y, type="lasso", normalize=FALSE)
beta.scaled <- coef(fit.lars)
pred.mat <- predict(fit.lars, X)$fit
pred.manual <- with(fit.lars, {
  ## Even when normalize=FALSE, the mean of the variables is
  ## subtracted away (but the variance stays the same).
  scale(X, meanx, FALSE) %*% t(beta.scaled) + mu
})
rbind(pred.mat[1,], pred.manual[1,])
stopifnot(pred.mat == pred.manual)
## If the prediction function for an input vector x is f(x) = b +
## beta'x, then the lars intercept b = meanx'beta + mu.
intercept.scaled <- with(fit.lars, -beta.scaled %*% meanx + mu)

lars.path.list <- list()
for(step.i in 1:nrow(beta.scaled)){
  coef.vec <- fit.lars$beta[step.i,]
  ## The loss in lars is half of the total squared error (0.5 * RSS)
  ## but the cost function in glmnet is half of the mean squared error
  ## (0.5 * RSS / N), so we divide the lars lambda value by N so that
  ## it is comparable to the glmnet lambda values.
  lambda <- fit.lars$lambda[step.i] / nrow(X)
  if(is.na(lambda))lambda <- 0
  lars.path.list[[paste(step.i)]] <- data.table(
    coef=c(intercept.scaled[[step.i]], coef.vec),
    arclength=sum(abs(coef.vec)),
    lambda, variable=c("(Intercept)", names(coef.vec)))
}
lars.path <- do.call(rbind, lars.path.list)

# fixing scale at 1 gives better opt. criteria than estimating scale
fit.iregnet <- iregnet(X, cbind(y, y), "gaussian")
print (fit.iregnet)
# fit.iregnet <- iregnet(X, cbind(y, y), "gaussian")
iregnet.path.list <- list()
for(i in 1:length(fit.iregnet$lambda)) {
	iregnet.path.list[[paste(i)]] <- data.table(
			lambda = fit.iregnet$lambda[[i]],
			coef = fit.iregnet$beta[2:9, i],
			arclength = sum(abs(fit.iregnet$beta[, i])),
			variable = colnames(X)
		)
}

iregnet.path <- do.call(rbind, iregnet.path.list)

# iregnet also uses 1/2 (MSE) cost function
lars.at.iregnet.lambda <- predict(
  fit.lars,
  s=fit.iregnet$lambda * nrow(X),
  type="coefficients", mode="lambda")

coef.mat.list <- list(
  iregnet=t(fit.iregnet$beta[2:(ncol(X)+1), ]),
  lars=lars.at.iregnet.lambda$coefficients)

lapply(coef.mat.list, head)
coef.mat.lambda <- fit.iregnet$lambda

# GLMNET
glmnet.thresh <- 1e-07
fit.glmnet <- glmnet(X, y, standardize=FALSE, thresh=glmnet.thresh, lambda=fit.iregnet$lambda)
glmnet.path.list <- list()
for(lambda.i in 1:ncol(fit.glmnet$beta)){
  coef.vec <- coef(fit.glmnet)[, lambda.i]
  glmnet.path.list[[paste(lambda.i)]] <- data.table(
    lambda=fit.glmnet$lambda[[lambda.i]],
    variable=names(coef.vec),
    coef=coef.vec,
    arclength=sum(abs(coef.vec[-1]))
    )
}
glmnet.path <- do.call(rbind, glmnet.path.list)

# lars.at.glmnet.lambda <- predict(
#   fit.lars,
#   s=fit.glmnet$lambda * nrow(X),
#   type="coefficients", mode="lambda")
# coef.mat.list <- list(
#   glmnet=t(fit.glmnet$beta),
#   lars=lars.at.glmnet.lambda$coefficients)
# lapply(coef.mat.list, head)
# coef.mat.lambda <- fit.glmnet$lambda
# 
coef.mat.list$glmnet <- t(fit.glmnet$beta)


## penalized pacakge.
fit.penalized.list <- with(fit.glmnet, penalized(
  y, X, lambda1=min(lambda), steps=length(lambda),
  standardize=FALSE))
fit.penalized <- penalized(
  y, X,
  lambda1=1,
  model="linear",
  standardize=FALSE)
ppred <- predict(fit.penalized, X)
residuals <- ppred[,1] - y
n <- length(residuals)
ss <- sum(residuals * residuals)
## The equation for the loglik below "is essentially the likelihood,
## but with sigma=1 plugged in" and comes from
## https://github.com/cran/penalized/blob/master/R/linear.R#L18
(loglik <- (-n/2) * (log(2*pi/n) + 1 + log(ss + .Machine$double.xmin)))
fit.penalized@loglik
# stopifnot(fit.penalized@loglik == loglik)
## without adding inside the log, still equal.
# stopifnot(fit.penalized@loglik == (-n/2) * (log(2*pi/n) + 1 + log(ss)))

## some checks using un-regularized glm solver.
fit.unregularized <- glm(y ~ X)
upred <- predict(fit.unregularized)
sum(dnorm(y, ppred[,1], ppred[,2], log=TRUE)) #not the loglik!
sum(dnorm(y, ppred[,1], 1, log=TRUE)) #not the loglik!

lassoDualityGap <- function
### Computes "Duality Gaps with Fenchel Duality" for the lasso
### problem, as explained in section D.2.3 of Julien Mairal's PHD
### thesis. The lasso problem is defined as min_w ||Xw-y||_2^2 / (2n)
### + lambda*||w||_1.
(y,
### Numeric vector of output labels (n observations)
 X,
### Numeric matrix of input features (n observations x p features).
 w,
### Numeric vector of weights (p features).
 lambda
### Numeric regularization parameter (scalar).
 ){
  stopifnot(is.matrix(X))
  stopifnot(is.numeric(X))
  stopifnot(is.numeric(y))
  stopifnot(is.numeric(w))
  stopifnot(is.numeric(lambda))
  stopifnot(length(lambda) == 1)
  stopifnot(nrow(X) == length(y))
  stopifnot(ncol(X) == length(w))
  dual.norm <- function(x)max(abs(x))
  primal.norm <- function(x)sum(abs(x))
  pred.vec <- X %*% w
  gradient <- (pred.vec - y)/nrow(X)
  inside.dual.norm <- t(X) %*% gradient
  dual.norm.value <- dual.norm(inside.dual.norm)
  ## Julien's thesis does the derivation for the square loss
  ## \tilde f(z) = ||y-z||_2^2 / (2n)
  ## which has the Fenchel/convex conjugate function
  ## \tilde f^*(k) = n ||k||_2^2 / 2 + k'y.
  ## so the total cost function is
  ## g(w) = ||y - Xw||_2^2 / (2n) + lambda * ||w||_1
  ## Note that is lambda not lambda1!
  grad.coef <- min(1, lambda/dual.norm.value)
  dual.vec <- grad.coef * gradient
  residual.vec <- y-pred.vec
  primal.cost <- sum(residual.vec * residual.vec) / 2 / nrow(X) + lambda * primal.norm(w)
  loss.conjugate <- nrow(X)/2 * sum(dual.vec * dual.vec) + sum(y * dual.vec)
  ## Theoretically the dual variables above (dual.vec) should never
  ## produce a penalty.conjugate.arg below which has an infinity norm
  ## greater than one, so penalty.conjugate should always be zero, and
  ## we do not need to include its term in the duality gap (lambda *
  ## penalty.conjugate). The penalty function is psi(w) = ||w||_1 and
  ## the penalty conjugate function is psi^*(gamma) = {0 if
  ## ||gamma||_Inf <= 0, and Inf otherwise}. The term which is always
  ## zero in the duality gap is lambda * psi^*(-Xk/lambda).
  ## penalty.conjugate.arg <- - t(X) %*% dual.vec / lambda
  ## penalty.conjugate <- dual.norm(penalty.conjugate.arg)
  ## stopifnot(penalty.conjugate <= 1)
  primal.cost + loss.conjugate
### The duality gap, a numeric scalar which gives an upper bound on
### the difference between this value of the cost function, and the
### best possible value: 0 <= cost(w) - cost(w^*) <= duality gap.
}
W0 <- cbind(rep(0, ncol(X)))
y.spams <- cbind(y)
spams.path.list <- list()
coef.mat.list$spams <- matrix(
  NA, length(coef.mat.lambda), ncol(X), dimnames=list(NULL, colnames(X)))
spams.tol <- 0.000001
for(lambda.i in seq_along(coef.mat.lambda)){
  lambda <- coef.mat.lambda[[lambda.i]]
  lambda1 <- lambda * nrow(X)
  W0 <- spams.fistaFlat(y.spams, X, W0, loss="square", regul="l1", lambda1=lambda1, tol=spams.tol)
  coef.mat.list$spams[lambda.i,] <- W0
  ## From the docs, this solves argmin_w 0.5 || y - Xw ||_2^2 + lambda1 || w ||_1.
  duality.gap <- lassoDualityGap(y, X, W0, lambda)
  cat(sprintf("lambda=%f, duality gap=%e\n",
              lambda, duality.gap))
  spams.path.list[[paste(lambda)]] <- data.table(
    lambda, variable=colnames(X), coef=as.numeric(W0), arclength=sum(abs(W0)))
}
spams.path <- do.call(rbind, spams.path.list)

timing.coef.list <- list()
timing.df <- microbenchmark(glmnet={
  timing.coef.list$glmnet <- glmnet(X, y, standardize=FALSE, thresh=glmnet.thresh)$beta
}, spams={
  timing.coef.list$spams <- matrix(
    NA, length(coef.mat.lambda), ncol(X), dimnames=list(NULL, colnames(X)))
  W0 <- cbind(rep(0, ncol(X)))
  for(lambda.i in seq_along(coef.mat.lambda)){
    lambda <- coef.mat.lambda[[lambda.i]]
    lambda1 <- lambda * nrow(X)
    W0 <- spams.fistaFlat(y.spams, X, W0, loss="square", regul="l1", lambda1=lambda1, tol=spams.tol)
    timing.coef.list$spams[lambda.i,] <- W0
  }
}, lars={
  fit.lars <- lars(X, y, type="lasso", normalize=FALSE)
  lars.at.glmnet.lambda <- predict(
    fit.lars,
    s=fit.glmnet$lambda * nrow(X),
    type="coefficients", mode="lambda")
  timing.coef.list$lars <- lars.at.glmnet.lambda$coefficients
})
timing.df

subdifferentialCriteria <- function
### Compute subdifferential optimality criteria for the elastic net problem.
(y,
### Numeric vector of output labels (n observations)
 X,
### Numeric matrix of input features (n observations x p features).
 w,
### Numeric vector of weights (p features).
 lambda,
### Numeric lambda regularization parameter (scalar).
 alpha=1
### Numeric elastic net parameter (between 0 and 1).
 ){
  stopifnot(is.matrix(X))
  stopifnot(is.numeric(X))
  stopifnot(is.numeric(y))
  stopifnot(is.numeric(w))
  stopifnot(is.numeric(lambda))
  stopifnot(length(lambda) == 1)
  stopifnot(nrow(X) == length(y))
  stopifnot(ncol(X) == length(w))
  stopifnot(is.numeric(alpha))
  stopifnot(length(alpha)==1)
  stopifnot(0 <= alpha && alpha <= 1)
  pred.vec <- X %*% w
  residual.vec <- pred.vec - y
  gradient.vec <- t(X) %*% residual.vec / nrow(X)
  positive.part <- function(x)ifelse(x<0, 0, x)
  common.term <- gradient.vec + lambda * (1-alpha) * weight.vec
  ifelse(weight.vec==0,
         positive.part(abs(common.term) - lambda*alpha),
         abs(common.term + lambda*alpha*sign(weight.vec)))
### p-vector of subdifferential optimality criteria. Lower values mean
### closer to the global optimum of the elastic net problem.
}

library(ggplot2)
addColumn <- function(dt, pkg){
  data.table(dt, pkg=factor(pkg, c("spams", "glmnet", "penalized", "iregnet")))
}
ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(pkg ~ ., scales="free")+
  geom_point(
    aes(lambda, coef, colour=variable),
    addColumn(glmnet.path, "glmnet"),
    shape=1
    )+
  geom_point(
    aes(lambda, coef, colour=variable),
    addColumn(spams.path, "spams"),
    shape=1
    )+
  geom_point(
    aes(lambda, coef, colour=variable),
    addColumn(spams.path, "iregnet"),
    shape=1
    )+
  geom_line(
    aes(lambda, coef, colour=variable),
    lars.path)

fig.coef.path <- ggplot()+
  ggtitle("lars (lines) is consistent with other solvers (dots)")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(pkg ~ ., scales="free")+
  geom_point(
    aes(arclength, coef, colour=variable),
    addColumn(glmnet.path, "glmnet"),
    shape=1
    )+
  geom_point(
    aes(arclength, coef, colour=variable),
    addColumn(spams.path, "spams"),
    shape=1
    )+
  geom_point(
    aes(arclength, coef, colour=variable),
    addColumn(iregnet.path, "iregnet"),
    shape=1
    )+
  geom_line(
    aes(arclength, coef, colour=variable),
    lars.path)
pdf("figure-lasso-criteria-path.pdf")
print(fig.coef.path)
dev.off()

convergence.list <- list()
for(pkg in names(coef.mat.list)){
  coef.mat <- coef.mat.list[[pkg]]
  for(lambda.i in seq_along(coef.mat.lambda)){
    weight.vec <- coef.mat[lambda.i, ]
    lambda <- coef.mat.lambda[[lambda.i]]
    param.criteria <- subdifferentialCriteria(y, X, weight.vec, lambda)
    criterion.value <- c(
      dualityGap=lassoDualityGap(y, X, weight.vec, lambda),
      subdifferentialL1=sum(abs(param.criteria)),
      subdifferentialLInf=max(abs(param.criteria)),
      subdifferentialL2=sqrt(sum(param.criteria * param.criteria)))
    convergence.list[[paste(pkg, lambda)]] <- data.table(
      pkg, lambda, criterion.value, criterion.name=names(criterion.value))
  }
}
convergence <- do.call(rbind, convergence.list)

## Accuracy of different solvers using the four different criteria.
with.legend <- ggplot()+
  ggtitle("optimality criteria using default solver thresholds")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(criterion.name ~ ., scales="free")+
  geom_point(aes(log10(lambda), log10(criterion.value), color=pkg),
             shape=1,
             data=convergence)
(with.labels <- direct.label(with.legend+xlim(-3.3, 0), "first.polygons"))
pdf("figure-lasso-criteria-all.pdf")
print(with.labels)
dev.off()

## Subdifferential optimality criterion is pretty much the same as the
## duality gap.
scatter.dt <- dcast(convergence, pkg + lambda ~ criterion.name, value.var="criterion.value")
abline.dt <- data.table(slope=1, intercept=0)
fig.scatter.leg <- ggplot()+
  ggtitle("subdifferential is consistent with duality gap")+
  coord_equal()+
  geom_abline(aes(slope=slope, intercept=intercept),
              data=abline.dt,
              color="grey")+
  geom_point(aes(log10(dualityGap), log10(subdifferentialL1), color=pkg),
             shape=1,
             data=scatter.dt)
fig.scatter <- direct.label(fig.scatter.leg)
pdf("figure-lasso-criteria.pdf")
print(fig.scatter)
dev.off()
