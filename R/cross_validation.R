#' Title
#'
#' @param A nxp heterogeneous observations (matrix)
#' @param groups vector of size n of group memberships (factor)
#' @param var.type vector of size p indicating types of columns (currently "gaussian", "binomial" or "poisson")
#' @param nb.boot number of folds for cross validation (integer)
#' @param thresh convergence tolerance (positive number)
#' @param maxit maximum number of iterations (positive integer)
#' @param lambda1.max maximum regularization parameter for nuclear norm penalty (positive number)
#' @param lambda2.max maximum regularization parameter for ell_1 norm penalty (positive number)
#' @param lambda1.min minimum regularization parameter for nuclear norm penalty (positive number)
#' @param lambda2.min minimum regularization parameter for ell_1 norm penalty (positive number)
#' @param length size of cross-validation grid (integer)
#' @param alpha0 initial value for the main effects (matrix of size (number of groups)xp)
#' @param theta0 initial value for the interactions (matrix of size nxp)
#' @param trace.it whether messages about convergence should be printed (boolean)
#'
#' @return
#' @export
#' @import softImpute
#'
#' @examples
crossval <- function(A, epsilon=0.1, nb.boot = 5, thresh = 1e-5,
                     maxit = 100, lambda1.max = NULL, lambda2.max = NULL,
                     lambda1.min = NULL, lambda2.min = NULL, length = 10,
                     S0 = NULL, L0 = NULL, trace.it = T){
  prob <- 0.2
  n <- nrow(A)
  omega <- !is.na(A)
  if(is.null(lambda2.max)) lambda2.max <- 2*max(sqrt(colSums(A^2, na.rm=T)))
  if(is.null(lambda1.max)) lambda1.max <- 2*max(softImpute(A)$d)
  if(is.null(lambda1.min)) lambda1.min <- lambda1.max / (100)
  if(is.null(lambda2.min)) lambda2.min <- lambda2.max / (100)
  lambda1.grid.log <- seq(log(lambda1.max), log(lambda1.min), length.out = length)
  lambda2.grid.log <- seq(log(lambda2.max), log(lambda2.min), length.out = length)
  if(is.null(S0)) S0 <- matrix(0,n,n)
  if(is.null(L0)) L0 <- matrix(0,n,n)
  S <- S0
  L <- L0
  iter <- 1
  obs.idx <- which(omega)
  na.func <- function(x, prob){
    xp <- x
    xp[sample(obs.idx, round(prob * sum(omega)))] <- NA
    xp
  }
  A.list <- lapply(1:nb.boot, function(i) na.func(A, prob))
  pred.list <- lapply(A.list, function(xx) (is.na(xx) * (!is.na(A))))
  res.cv <- list()
  for(k in 1:nb.boot){
    print(k)
    res.cv[[k]] <- lapply(1:length(lambda1.grid.log), function(i) lapply(1:length(lambda2.grid.log),
                                                                         function(j) gsbm_mcgd(A.list[[k]], lambda1 = exp(lambda1.grid.log[i]),
                                                                                               lambda2 = exp(lambda2.grid.log[j]), epsilon=epsilon,
                                                                                               U = NULL, maxit = maxit, thresh = thresh,
                                                                                               S0 = S0, L0 = L0, R0 = NULL, trace.it = trace.it)))
    om <- pred.list[[k]]
    res.cv[[k]] <- lapply(1:length(lambda1.grid.log), function(i) lapply(1:length(lambda2.grid.log),
                                                                         function(j) sum(((res.cv[[k]][[i]][[j]]$S*om+t(res.cv[[k]][[i]][[j]]$S*om)+res.cv[[k]][[i]][[j]]$L*om) - A*om)^2, na.rm = T)))

  }
  for(k in 1:nb.boot){
    res.cv[[k]] <- sapply(1:length(lambda1.grid.log), function(i) unlist(res.cv[[k]][[i]]))
    colnames(res.cv[[k]]) <- exp(lambda2.grid.log)
    rownames(res.cv[[k]]) <- exp(lambda1.grid.log)
  }
  error <- Reduce('+', res.cv)/nb.boot
  idx <- which(error == min(error, na.rm =T), arr.ind = TRUE)
  lambda1.cv <- exp(lambda1.grid.log)[idx[1]]
  lambda2.cv <- exp(lambda2.grid.log)[idx[2]]
  estim.cv <- gsbm_mcgd(A, lambda1 = lambda1.cv, lambda2 = lambda2.cv,
                        epsilon=epsilon, U = NULL, maxit = maxit, thresh = thresh,
                        S0 = S0, L0 = L0, R0 = NULL, trace.it = trace.it)
  return(list(lambda1 = lambda1.cv, lambda2 = lambda2.cv, estim.cv = estim.cv, error = error,
              lambda1.grid = exp(lambda1.grid.log), lambda2.grid = exp(lambda2.grid.log)))

}
