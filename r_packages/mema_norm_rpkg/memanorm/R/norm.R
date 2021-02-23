#' Normalize feature matrices.
#'
#' @param f_mats feature matrices
#' @param adj_idx indices of which features to use to learn spatial effects
#' @param row_rep row replicate factor
#' @param col_rep column replicate factor
#' @param d dimension of spatial effects to remove
#' @param row_impute impute missing values in f_mats row-wise as col_rep grouped means
#' @param row_center subtract mean from f_mats rows
#' @param verbose print intermediate results
#' @param full.results return intermediate results. If FALSE just returns adjusted matrices.
#' @export
remove_spatial <- function(f_mats, adj_idx, row_rep, col_rep, d, row_impute = TRUE, 
    row_center = TRUE, verbose = TRUE, full.results = FALSE) {
    # pre-proscessing
    if (any(sapply(f_mats, class) != "matrix")) 
        f_mats <- lapply(f_mats, as.matrix)
    if (row_impute) 
        f_mats <- lapplyp(f_mats, row_impute_fn, col_rep, verbose = verbose, .name = "Imputing")
    if (row_center) 
        f_mats <- lapply(f_mats, function(mat) t(scale(t(mat), scale = FALSE, center = TRUE)))
    # replicate design matrices
    L <- stats::model.matrix(~-1 + row_rep)
    Et <- stats::model.matrix(~-1 + col_rep)
    # estimate S
    PL <- proj_mat(L)
    f_mats_L <- lapply(f_mats[adj_idx], function(mat) residop(mat, P = PL))
    if (d > 0) {
        asvd <- average_svd_na(f_mats_L, nu = 0, nv = d)
        s_hat <- t(asvd$v)
    } else {
        s_hat <- array(0, c(1, ncol(f_mats_L[[1]])))
    }
    # estimate Wi
    s_hat_Et <- residop(t(s_hat), Et)
    w_maker <- s_hat_Et %*% MASS::ginv(t(s_hat_Et) %*% s_hat_Et)
    w_hats <- lapply(f_mats, function(mat) mat %*% w_maker)
    # remove spatial effects
    rmv <- lapply(w_hats, function(w) w %*% s_hat)
    f_adj <- lapply(seq_along(f_mats), function(i) f_mats[[i]] - rmv[[i]])
    names(f_adj) <- names(f_mats)
    if (full.results) 
        return(list(f_adj = f_adj, w_hats = w_hats, s_hat = s_hat))
    return(list(f_adj = f_adj))
}

# R(B)%*%A; P = proj_mat(B)
residop <- function(A, B = NULL, P = NULL) {
    if (is.null(P)) 
        P <- proj_mat(B)
    return(A - P %*% A)
}

proj_mat <- function(B) {
    return(B %*% MASS::ginv(t(B) %*% B) %*% t(B))
}

row_impute_fn <- function(Y, col_rep) {
    for (rep in col_rep) {
        ii <- which(col_rep == rep)
        I <- t(apply(Y[, ii], 1, impute_fn))
        Y[, ii] <- I
    }
    return(Y)
}

impute_fn <- function(x) {
    x[is.na(x)] <- mean(x, na.rm = TRUE)
    return(x)
}
