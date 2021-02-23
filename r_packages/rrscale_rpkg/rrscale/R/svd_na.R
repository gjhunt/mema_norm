#' The completed SVD 
#' 
#' This calculates right and left singular vectors of a data matrix possibly containing missing values.
#'
#' @param X the data matrix of which to calcluate the completed SVD. 
#' @param nu the number of left singular vectors to calculate
#' @param nv the nubmer of right singular vectors to calculate
#' @examples
#' Y <- rnorm(10)%*%t(rnorm(10))
#' Y[1,1] <- NA
#' svdc.out <- svdc(Y)
#' @export
svdc <- function(X, nu = NULL, nv = NULL) {
    X <- as.matrix(X)
    eU <- eigen(gram(X, chir = "left", na = 0), symmetric = TRUE)
    eV <- eigen(gram(X, chir = "right", na = 0), symmetric = TRUE)
    eU$vectors <- eU$vectors[, c(sort(which(eU$values > 0)), sort(which(eU$values < 
        0)), sort(which(eU$values == 0)))]
    eU$values <- sapply(eU$values, max, 0)
    eV$vectors <- eV$vectors[, c(sort(which(eV$values > 0)), sort(which(eV$values < 
        0)), sort(which(eV$values == 0)))]
    eV$values <- sapply(eV$values, max, 0)
    if (is.null(c(nu, nv))) 
        nu <- nv <- min(dim(X))
    if (is.null(nu)) 
        nu <- sum(eU$values > 0)
    if (is.null(nv)) 
        nv <- sum(eV$values > 0)
    return(list(u = eU$vectors[, 1:nu], du = eU$values[1:nu], v = eV$vectors[, 1:nv], 
        dv = eV$values[1:nv], missing = which(is.na(X))))
}

gram <- function(A, chir = "left", na = NA, project = TRUE) {
    if (chir == "right") 
        A <- t(A)
    gram <- A %:% t(A)
    gram_zero <- gram
    gram_zero[is.na(gram_zero)] <- 0
    if (project) 
        gram_zero <- spsd_project(gram_zero)
    gram_zero[is.na(gram)] <- na
    return(gram_zero)
}

spsd_project <- function(A) {
    e <- eigen(A, symmetric = TRUE)
    dp <- sapply(e$values, max, 0)
    A_spsd <- e$vectors %*% diag(dp) %*% t(e$vectors)
    return(A_spsd)
}

`%.%` <- function(A, B) {
    C <- array(NA, c(nrow(A), ncol(B)))
    for (i in 1:nrow(A)) {
        for (j in 1:ncol(B)) {
            prds <- A[i, , drop = TRUE] * B[, j, drop = TRUE]
            if (!all(is.na(prds))) {
                C[i, j] <- mean(prds, na.rm = TRUE)
            }
        }
    }
    return(ncol(A) * C)
}

`%:%` <- function(A, B) {
    A <- as.matrix(A)
    B <- as.matrix(B)
    A_na <- array(0, dim(A))
    A_na[!is.na(A)] <- 1
    B_na <- array(0, dim(B))
    B_na[!is.na(B)] <- 1
    M <- A_na %*% B_na
    A[is.na(A)] <- 0
    B[is.na(B)] <- 0
    return((A %*% B)/M * ncol(A))
}

mean_mtx_na <- function(mtx_list) {
    mtx_list <- lapply(mtx_list, as.matrix)
    mtx_na <- lapply(mtx_list, function(x) {
        M <- array(0, dim(x))
        M[!is.na(x)] <- 1
        return(M)
    })
    M <- Reduce("+", mtx_na)
    M[M == 0] <- 1
    mtx_list <- lapply(mtx_list, function(x) {
        x[is.na(x)] <- 0
        return(x)
    })
    return(Reduce("+", mtx_list)/M)
}
