#' Traditional box-cox power transformation. Accepts one real parameter
#'
#' \itemize{
#' \item{T} the transformation with arguments Y, the data, lambda the parameter, and boolean inverse to calculate inverse transformation
#' \item{T_deriv} the transformation with arguments Y, the data, lambda the parameter.
#' }
#' @export
box_cox <- list(T = function(Y, lambda, inverse = FALSE) {
    if (!inverse) {
        if (lambda == 0) {
            Y <- log(Y)
        } else {
            Y <- (Y^lambda - 1)/lambda
        }
    } else {
        if (lambda == 0) {
            Y <- exp(Y)
        } else {
            Y <- (lambda * Y + 1)^(1/lambda)
        }
    }
    
    return(Y)
}, T_deriv = function(Y, lambda) {
    if (lambda == 0) {
        Y <- 1/Y
    } else {
        Y <- Y^(lambda - 1)
    }
    
    return(Y)
})

#' Simple power transformation
#'
#' \itemize{
#' \item{T} the transformation with arguments Y, the data, lambda the parameter, and boolean inverse to calculate inverse transformation.
#' \item{T_deriv} the transformation with arguments Y, the data, lambda the parameter.
#' }
#' @export
power <- list(T = function(Y, lambda, inverse = FALSE) {
    if (!inverse) {
        Y <- Y^lambda
    } else {
        Y <- Y^(1/lambda)
    }
    return(Y)
}, T_deriv = function(Y, lambda) {
    Y <- lambda * Y^(lambda - 1)
    return(Y)
})

#' Box-cox transformation of shifted variable
#'
#' \itemize{
#' \item{T} the transformation with arguments Y, the data, lambda the parameter, and boolean inverse to calculate inverse transformation. The parameter lambda has two real elements (1) the power and (2) the additive shift to the data.  
#' \item{T_deriv} the transformation with arguments Y, the data, lambda the parameter.
#' }
#' @export
box_cox_shift <- list(T = function(Y, lambda, inverse = FALSE) {
    if (!inverse) {
        if (lambda[1] == 0) {
            Y <- log(lambda[2] + Y)
        } else {
            Y <- ((lambda[2] + Y)^lambda[1] - 1)/lambda[1]
        }
    } else {
        if (lambd[1] == 0) {
            Y <- exp(Y) - lambda[2]
        } else {
            Y <- (lambda[1] * Y + 1)^(1/lambda) - lambda[2]
        }
    }
    
    return(Y)
}, T_deriv = function(Y, lambda) {
    if (lambda[1] == 0) {
        Y <- 1/(lambda[2] + Y)
    } else {
        Y <- (lambda[2] + Y)^(lambda[1] - 1)
    }
    
    return(Y)
})

#' Arc-hyperbolic-sine transformation
#'
#' \itemize{
#' \item{T} the transformation with arguments Y, the data, lambda the parameter, and boolean inverse to calculate inverse transformation.
#' \item{T_deriv} the transformation with arguments Y, the data, lambda the parameter.
#' }
#' @export
asinh <- list(T = function(Y, lambda, inverse = FALSE) {
    if (!inverse) {
        if (lambda == 0) {
            Y <- Y
        } else {
            Y <- asinh(lambda * Y)/lambda
        }
    } else {
        if (lambda == 0) {
            Y <- Y
        } else {
            Y <- sinh(lambda * Y)/lambda
        }
    }
    
    return(Y)
}, T_deriv = function(Y, lambda) {
    if (lambda == 0) {
        if (!is.null(dim(Y))) Y <- array(1, dim = dim(Y)) else Y <- rep(1, length(Y))
    } else {
        Y <- 1/sqrt(1 + lambda^2 * Y^2)
    }
    
    return(Y)
})

#' Log of the traditional box-cox transformation
#'
#' \itemize{
#' \item{T} the transformation with arguments Y, the data, lambda the parameter, and boolean inverse to calculate inverse transformation.
#' \item{T_deriv} the transformation with arguments Y, the data, lambda the parameter.
#' }
#' @export
log_box_cox <- list(T = function(Y, lambda, inverse = FALSE) {
    if (!inverse) {
        if (lambda == 0) {
            Y <- log(log(Y))
        } else {
            Y <- log((Y^lambda - 1)/lambda)
        }
    } else {
        if (lambda == 0) {
            Y <- exp(exp(Y))
        } else {
            Y <- (lambda * exp(Y) + 1)^(1/lambda)
        }
    }
    
    return(Y)
}, T_deriv = function(Y, lambda) {
    if (lambda == 0) {
        Y <- 1/(Y * log(Y))
    } else {
        Y <- (lambda * Y^(lambda - 1))/(Y^lambda - 1)
    }
    
    return(Y)
})

#' Box-cox transformation with a shift of 1 added to the data
#'
#' \itemize{
#' \item{T} the transformation with arguments Y, the data, lambda the parameter, and boolean inverse to calculate inverse transformation.
#' \item{T_deriv} the transformation with arguments Y, the data, lambda the parameter.
#' }
#' @export
box_cox_plus1 <- list(T = function(Y, lambda, inverse = FALSE) {
    if (!inverse) {
        if (lambda == 0) {
            Y <- log(1 + Y)
        } else {
            Y <- ((1 + Y)^lambda - 1)/lambda
        }
    } else {
        if (lambda == 0) {
            Y <- exp(Y) - 1
        } else {
            Y <- (lambda * Y + 1)^(1/lambda) - 1
        }
    }
    
    return(Y)
}, T_deriv = function(Y, lambda) {
    if (lambda == 0) {
        Y <- 1/(1 + Y)
    } else {
        Y <- (1 + Y)^(lambda - 1)
    }
    
    return(Y)
})

#' A generalized box-cox transformation that can handle negative data
#'
#' \itemize{
#' \item{T} the transformation with arguments Y, the data, lambda the parameter, and boolean inverse to calculate inverse transformation.
#' \item{T_deriv} the transformation with arguments Y, the data, lambda the parameter.
#' }
#' @export
box_cox_negative <- list(T = function(Y, lambda, inverse = FALSE) {
    if (!inverse) {
        Y <- (sign(Y) * abs(Y)^lambda - 1)/lambda
    } else {
        Y <- sign(lambda * x + 1) * abs(lambda * x + 1)
    }
    
    return(Y)
}, T_deriv = function(Y, lambda) {
    if (lambda == 0) {
        Y <- Inf
    } else {
        Y <- abs(Y)^(lambda - 1)
    }
    
    return(Y)
})

#' Box-cox transformation with the data shifted so that it is positive
#'
#' \itemize{
#' \item{T} the transformation with arguments Y, the data, lambda the parameter, and boolean inverse to calculate inverse transformation.
#' \item{T_deriv} the transformation with arguments Y, the data, lambda the parameter.
#' }
#' @export
box_cox_plusmin <- list(T = function(Y, lambda, inverse = FALSE) {
    eps <- 1e-05
    if (min(Y) > 0) {
        l2 <- 0
    } else if (min(Y) < 0) {
        l2 <- -min(Y)
    } else if (min(Y) == 0) {
        l2 <- eps
    }
    lambda <- c(lambda, l2)
    if (!inverse) {
        if (lambda[1] == 0) {
            Y <- log(lambda[2] + Y)
        } else {
            Y <- ((lambda[2] + Y)^lambda[1] - 1)/lambda[1]
        }
    } else {
        if (lambd[1] == 0) {
            Y <- exp(Y) - lambda[2]
        } else {
            Y <- (lambda[1] * Y + 1)^(1/lambda) - lambda[2]
        }
    }
    
    return(Y)
}, T_deriv = function(Y, lambda) {
    eps <- 1e-05
    if (min(Y) > 0) {
        l2 <- 0
    } else if (min(Y) < 0) {
        l2 <- -min(Y)
    } else if (min(Y) == 0) {
        l2 <- eps
    }
    lambda <- c(lambda, l2)
    if (lambda[1] == 0) {
        Y <- 1/(lambda[2] + Y)
    } else {
        Y <- (lambda[2] + Y)^(lambda[1] - 1)
    }
    
    return(Y)
})

#' Exponential of the tranditional box-cox transformation
#'
#' \itemize{
#' \item{T} the transformation with arguments Y, the data, lambda the parameter, and boolean inverse to calculate inverse transformation.
#' \item{T_deriv} the transformation with arguments Y, the data, lambda the parameter.
#' }
#' @export
box_cox_exp <- list(T = function(Y, lambda, inverse = FALSE) {
    if (!inverse) {
        if (lambda == 0) {
            Y <- Y
        } else {
            Y <- (exp(lambda * Y) - 1)/lambda
        }
    } else {
        if (lambda == 0) {
            Y <- Y
        } else {
            Y <- log(lambda * Y + 1)/lambda
        }
    }
    
    return(Y)
}, T_deriv = function(Y, lambda) {
    if (lambda == 0) {
        Y <- Y
    } else {
        Y <- exp(lambda * Y)
    }
    
    return(Y)
})


# yj helper functions
yjpower <- function(y, lambda) {
    if ((lambda != 0) & (y >= 0)) {
        return(((y + 1)^lambda - 1)/lambda)
    }
    if ((lambda == 0) & (y >= 0)) {
        return(log(y + 1))
    }
    if ((lambda != 2) & (y < 0)) {
        return(-((-y + 1)^(2 - lambda) - 1)/(2 - lambda))
    }
    if ((lambda == 2) & (y < 0)) {
        return(-log(-y + 1))
    }
    print("here")
}

yjpower_inv <- function(y, lambda) {
    if ((lambda != 0) & (y >= 0)) {
        return((y * lambda + 1)^(1/lambda) - 1)
    }
    if ((lambda == 0) & (y >= 0)) {
        return(exp(y) - 1)
    }
    if ((lambda != 2) & (y < 0)) {
        return(-((-y * (2 - lambda) + 1)^(1/(2 - lambda)) + 1))
    }
    if ((lambda == 2) & (y < 0)) {
        return(exp(-y) - 1)
    }
    print("here")
}

yjpower_deriv <- function(y, lambda) {
    if ((lambda != 0) & (y >= 0)) {
        return((y + 1)^(lambda - 1))
    }
    if ((lambda == 0) & (y >= 0)) {
        return(1/(1 + y))
    }
    if ((lambda != 2) & (y < 0)) {
        return((-y + 1)^(1 - lambda))
    }
    if ((lambda == 2) & (y < 0)) {
        return(1/(1 - y))
    }
    print("here")
}


#' Yeo-Johnson transformation.
#'
#' \itemize{
#' \item{T} the transformation with arguments Y, the data, lambda the parameter, and boolean inverse to calculate inverse transformation
#' \item{T_deriv} the transformation with arguments Y, the data, lambda the parameter.
#' }
#' @export
yeo_johnson <- list(T = function(Y, lambda, inverse = FALSE) {
    if (!inverse) {
        if (lambda != 0) Y[which(Y >= 0)] <- ((Y[which(Y >= 0)] + 1)^lambda - 1)/lambda
        if (lambda == 0) Y[which(Y >= 0)] <- log(1 + Y[which(Y >= 0)])
        if (lambda != 2) Y[which(Y < 0)] <- -((-Y[which(Y < 0)] + 1)^(2 - lambda) - 
            1)/(2 - lambda)
        if (lambda == 2) Y[which(Y < 0)] <- -log(-Y[which(Y < 0)] + 1)
        
        return(Y)
    } else {
        if (is.null(dim(Y))) return(sapply(Y, yjpower_inv, lambda))
        return(apply(Y, c(1, 2), yjpower_inv, lambda))
    }
}, T_deriv = function(Y, lambda) {
    if (lambda != 0) Y[which(Y >= 0)] <- (Y[which(Y >= 0)] + 1)^(lambda - 1)
    if (lambda == 0) Y[which(Y >= 0)] <- 1/(1 + Y[which(Y >= 0)])
    if (lambda != 2) Y[which(Y < 0)] <- (-Y[which(Y < 0)] + 1)^(1 - lambda)
    if (lambda == 2) Y[which(Y < 0)] <- 1/(1 - Y[which(Y < 0)])
    return(Y)
})

old_yeo_johnson <- list(T = function(Y, lambda, inverse = FALSE) {
    if (!inverse) {
        if (is.null(dim(Y))) return(sapply(Y, yjpower, lambda))
        return(apply(Y, c(1, 2), yjpower, lambda))
    } else {
        if (is.null(dim(Y))) return(sapply(Y, yjpower_inv, lambda))
        return(apply(Y, c(1, 2), yjpower_inv, lambda))
    }
}, T_deriv = function(Y, lambda) {
    if (is.null(dim(Y))) return(sapply(Y, yjpower_deriv, lambda))
    return(apply(Y, c(1, 2), yjpower_deriv, lambda))
})


#' List possible transformations
#' 
#' Returns list of transformations. Each transformation is a transformation function (``T'') accepting a parameter and the derivative of this transformation function (``T_deriv'').
#' @export
list_transformations <- function() {
    return(list(box_cox = box_cox, power = power, box_cox_shift = box_cox_shift, 
        asinh = asinh, log_box_cox = log_box_cox, box_cox_plus1 = box_cox_plus1, 
        box_cox_negative = box_cox_negative, box_cox_plusmin = box_cox_plusmin, box_cox_exp = box_cox_exp))
}
