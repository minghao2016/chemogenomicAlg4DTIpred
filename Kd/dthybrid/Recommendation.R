.doCall <- function (fun, args) {
  if ((is.character(fun) && length(fun) == 1) || is.name(fun)) {
    fun <- get(as.character(fun), envir=.GlobalEnv, mode="function")
  }
  do.call("fun", lapply(args, enquote))
}

.splitRows <- function (x, ncl) {
  lapply(.splitIndices(nrow(x), ncl), 
         function(i) x[i, , drop=FALSE])
}

.splitIndices <- function (nx, ncl) {
  i <- 1:nx
  if (ncl == 1) {
    i
  } else {
    fuzz <- min((nx - 1)/1000, 0.4 * nx/ncl)
    breaks <- seq(1 - fuzz, nx + fuzz, length=ncl + 1)
    structure(split(i, cut(i, breaks)), names=NULL)
  }
}

.parMM <- function (cl, A, B) {
  if (!all(is.na(cl)) && is.object(cl)) {
    R <- .doCall(rbind, 
                clusterApply(cl=cl, x=.splitRows(A, length(cl)), 
                             get("%*%"), B))
  } else {
    R <- A %*% B
  }
  return (R)
}

graphWeights <- function (n, m, A, lambda=0.5, alpha=0.5, S=NA, S1=NA, cl=NA) {
  if (nrow(A) != n || ncol(A) != m) {
    stop("The matrix A should be an n by m matrix.")
  }
  
  has.similarity <- (!all(is.na(S)) && is.matrix(S) && !all(is.na(S1)) && is.matrix(S1))
  
  if (has.similarity) {
    if (nrow(S1) != m || ncol(S1) != m) {
      stop("The matrix S1 should be an m by m matrix.")
    }
    if (nrow(S) != n || ncol(S) != n) {
      stop("The matrix S should be an n by n matrix.")
    }
  }
  
  Ky <- diag(1/colSums(A))
  Ky[is.infinite(Ky) | is.na(Ky)] <- 0 #BugFix: 1/0=Infinite replaced with 0
  
  kx <- rowSums(A)
  Nx <- 1/(matrix(kx, nrow=n, ncol=n, byrow=TRUE)^(lambda) * 
             matrix(kx, nrow=n, ncol=n, byrow=FALSE)^(1-lambda))
  Nx[is.infinite(Nx) | is.na(Nx)] <- 0 #BugFix: 1/0=Infinite replaced with 0
  kx[is.infinite(kx) | is.na(kx)] <- 0 #BugFix: 1/0=Infinite replaced with 0
  
  W <- t(.parMM(cl, A, Ky))
  W <- .parMM(cl, A, W)  
  W <- Nx * W
  rownames(W) <- rownames(A)
  colnames(W) <- rownames(A)
  
  if (has.similarity) {
    X5 <- .parMM(cl, A, S1)
    X6 <- .parMM(cl, X5, t(A))
    X7 <- .parMM(cl, A, matrix(1, nrow=m, ncol=m))
    X8 <- .parMM(cl, X7, t(A))
    S2 <- X6 / X8
    W  <- W * (1 + (alpha * S) + ((1-alpha) * S2))
  }
  
  W[is.nan(W)] <- 0 #This should never happen
  return (W)
}

computeRecommendation <- function (A, lambda=0.5, alpha=0.5, S=NA, S1=NA, cl=NA) {
  n <- nrow(A)
  m <- ncol(A)
  W <- graphWeights(n=n, m=m, A=A, lambda=lambda, alpha=alpha, S=S, S1=S1, cl=cl)
  R <- .parMM(cl, W, A)
  return (R)
}