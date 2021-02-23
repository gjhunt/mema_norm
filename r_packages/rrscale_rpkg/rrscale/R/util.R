#' Centers the data column-wise
#'
#' @param x the data. 
#' @export
center <- function(x) scale(x, scale = FALSE)

#' Winsorizes the data
#'
#' @param x the data.
#' @param fraction the top and bottom quantiles to cap.
#' @examples
#' Y <- rlnorm(10)%*%t(rlnorm(10))
#' Yw <- winsor(Y,1E-2)
#' @export
winsor <- function(x, fraction = 0.01) {
    if (length(fraction) != 1 || fraction < 0 || fraction > 0.5) {
        stop("Fraction must be in the interval (.5,1)")
    }
    lim <- stats::quantile(x, probs = c(fraction, 1 - fraction), na.rm = TRUE)
    x[x < lim[1]] <- lim[1]
    x[x > lim[2]] <- lim[2]
    x
}

#' Calculate the geometric mean
#'
#' @examples
#' Y <- rlnorm(10)
#' gm <- gm_mean(Y)
#' @param x the data.
#' @export
gm_mean <- function(x) {
    exp(mean(log(x[x > 0])))
}
