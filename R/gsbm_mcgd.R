#' Fit sparse main effects plus low-rank interactions model
#' with least-squares data fitting term (numeric variables only)
#'
#' @param A nxn adjacency matrix
#' @param lambda1 regularization parameter for nuclear norm penalty (positive number)
#' @param lambda2 regularization parameter for ell_2,1 penalty (positive number)
#' @param U lower bound on nuclear norm (positive number, if NULL, default method is applied)
#' @param maxit maximum number of iterations (positive integer)
#' @param thresh convergence tolerance (positive number)
#' @param S0 initial value for the sparse component
#' @param L0 initial value for the low-rank component
#' @param R0 lower bound on nuclear norm of L0 (positive number, if NULL, default method is applied)
#' @param trace.it whether messages about convergence should be printed (boolean)
#'
#' @return
#' @importFrom stats aggregate
#' @export
#'
#' @examples
gsbm_mcgd <- function(A, lambda1, lambda2, epsilon=0.1, U = NULL, maxit = 100, thresh = 1e-6,
                      S0 = NULL, L0 = NULL, R0 = NULL, trace.it = T){
  n <- nrow(A)
  Omega <- !is.na(A)
  if(is.null(S0)) {
    S0 <- matrix(0,n,n)
  }
  if(is.null(L0)) L0 <- matrix(0,n,n)
  if(is.null(R0)) R0 <- svd(L0, nu = 1, nv = 1)$d[1]
  if(is.null(U)) {
    U <- (1/2)*sum((A - S0)^2, na.rm = T)/lambda1
  }
  S <- S0
  L <- L0
  R <- R0
  objective <- 1/2*sum(A^2, na.rm=T)
  error <- 1
  iter <- 0
  A0 <- A
  A[is.na(A)] <- 0
  while((error > thresh) && (iter < maxit)){
    iter <- iter + 1
    S.tmp <- S
    L.tmp <- L
    R.tmp <- R
    G_S <- -2*Omega*(A-L-S-t(S))+epsilon*S # gradient wrt S
    obj0 <- (1/2)*sum((A0-S-t(S)-L)^2, na.rm = T)+lambda1*R+lambda2*sum(sqrt(colSums(S^2)))+epsilon*(norm(L, type="F")^2+norm(S, type="F")^2)
    step <- 1
    flag <- TRUE
    while(flag){
      step <- 0.5*step
      mat <- sapply(1:n, function(j){
        max(1-step*lambda2/sqrt(sum((S[,j]-step*G_S[,j])^2, na.rm=T)),0)*(S[,j]-step*G_S[,j])
      })
      mat <- pmax(mat,0)
      obj <- (1/2)*sum((A0-mat-t(mat)-L)^2, na.rm = T)+lambda1*R+lambda2*sum(sqrt(colSums(mat^2)))+epsilon*(norm(L, type="F")^2+norm(mat, type="F")^2)
      flag <- obj > obj0
    }
    S <- mat
    G_L <- -Omega*(A - S - t(S) - L) + epsilon*L
    obj0 <- (1/2)*sum((A0-S-t(S)-L)^2, na.rm = T)+lambda1*R+lambda2*sum(sqrt(colSums(S^2)))+epsilon*(norm(L, type="F")^2+norm(S, type="F")^2)
    step <- 2
    flag <- TRUE
    svdL <- svd(G_L, nu = 1, nv = 1)
    D_t <- - svdL$u%*%t(svdL$v)
    while(flag){
      step <- 0.5*step
      if(lambda1 >= - sum(D_t*G_L)){
        R_tilde <- 0
        L_tilde <- matrix(0,n,n)
        L <- L.tmp + step*(L_tilde - L.tmp)
        L <- (L+t(L))/2 # to avoid propagation of numerical errors
        R <- R.tmp + step*(R_tilde - R.tmp)
      } else{
        R_tilde <- U
        L_tilde <- U*D_t
        #beta <- min(1, (sum((L - L_tilde)*L_tilde)+(R - R_tilde)*lambda1)/norm(L_tilde - L, type = "F")^2)
        L <- L.tmp + step*(L_tilde - L.tmp)
        L <- (L+t(L))/2 # to avoid propagation of numerical errors
        R <- R.tmp + step*(R_tilde - R.tmp)
      }
      obj <- (1/2)*sum((A0-S-t(S)-L)^2, na.rm = T)+lambda1*R+lambda2*sum(sqrt(colSums(S^2)))+epsilon*(norm(L, type="F")^2+norm(S, type="F")^2)
      flag <- obj > obj0
    }
    obj <- (1/2)*sum((A0-S-t(S)-L)^2, na.rm = T)+lambda1*R+lambda2*sum(sqrt(colSums(S^2)))+epsilon*(norm(L, type="F")^2+norm(S, type="F")^2)
    objective <- c(objective, obj)
    U <- obj/lambda1
    if(iter == 1) {
      error <- 1
    } else{
      error <- abs(objective[iter]-objective[iter - 1]) /abs(objective[iter])
    }
    if(trace.it && (iter%%10 == 0) ){
      print(paste("iter ", iter, ": error ", error, " - objective: ", objective[iter]))
    }
  }
  return(list(A = A0, L = L, S = S, objective = objective, R = R,
              iter = iter))
}
