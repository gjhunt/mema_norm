context("always")

# deoptim
set.seed(3251991)
A <- rlnorm(10)
B <- rlnorm(10)
Y <- t(t(B)) %*% t(A)
out <- rrscale(Y, trans_list = list(box_cox_negative = box_cox_negative, asinh = asinh), 
    lims_list = list(box_cox_negative = c(-100, 100), asinh = list(0, 100)), seed = 3251991, 
    verbose = FALSE)
out2 <- out[c("pars", "par_hat", "NT", "RR", "G", "Z", "O", "T_name")]
out2$NT <- as.matrix(out2$NT)
out2$G <- as.matrix(out2$G)
out2$Z <- as.matrix(out2$Z)
out2$O <- as.matrix(out2$O)
out2$RR <- as.matrix(out2$RR)
expect_equal_to_reference(out2, "basic_ref10.rds", tolerance = 1e-05, check.attributes = FALSE)

# nloptr
set.seed(3251991)
A <- rlnorm(10)
B <- rlnorm(10)
Y <- t(t(B)) %*% t(A)
out <- rrscale(Y, trans_list = list(box_cox_negative = box_cox_negative, asinh = asinh), 
    lims_list = list(box_cox_negative = c(-100, 100), asinh = list(0, 100)), seed = 3251991, 
    verbose = FALSE, opt_method = "nloptr", opts = TRUE)
out2 <- out[c("pars", "par_hat", "NT", "RR", "G", "Z", "O", "T_name")]
out2$NT <- as.matrix(out2$NT)
out2$G <- as.matrix(out2$G)
out2$Z <- as.matrix(out2$Z)
out2$O <- as.matrix(out2$O)
out2$RR <- as.matrix(out2$RR)
expect_equal_to_reference(out2, "basic_ref10.rds", tolerance = 1e-05, check.attributes = FALSE)
