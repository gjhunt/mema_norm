context("test-basic")

.runThisTest <- Sys.getenv("RunAllRRTests") == "yes"

if (.runThisTest) {
    
    # RR transformations Z, G, O and RR
    set.seed(3251991)
    A <- rlnorm(100)
    B <- rlnorm(100)
    Y <- t(t(B)) %*% t(A)
    out <- rrscale(Y, trans_list = list(box_cox_negative = box_cox_negative, asinh = asinh), 
        lims_list = list(box_cox_negative = c(-100, 100), asinh = list(0, 100)), 
        seed = 3251991, verbose = FALSE)
    out2 <- out[c("pars", "par_hat", "NT", "RR", "G", "Z", "O", "T_name")]
    out2$NT <- as.matrix(out2$NT)
    expect_equal_to_reference(out2, "basic_ref.rds", tolerance = 1e-05, check.attributes = FALSE)
    
    Y_tmp <- data.frame(Y)
    expect_equal(Y_tmp, out$rr(Y, G = FALSE, Z = FALSE, O = FALSE), check.attributes = FALSE)
    Yw <- as.matrix(winsor(Y_tmp, 0.001))
    mu <- mean(as.matrix(Yw), na.rm = TRUE)
    Ywc <- Yw - mu
    norm <- sqrt(mean(Ywc^2, na.rm = TRUE))
    Y_tmp <- (Y_tmp - mu)/norm
    Z <- Y_tmp
    expect_equal(Z, out$rr(Y, G = FALSE, Z = TRUE, O = FALSE), check.attributes = FALSE)
    O <- data.frame(Y)
    O[abs(Z) > 4] <- NA
    expect_equal(O, out$rr(Y, G = FALSE, Z = FALSE, O = TRUE), check.attributes = FALSE)
    
    ## trans_list test
    out2 <- rrscale(Y, trans_list = list(asinh = asinh), lims_list = list(box_cox_negative = c(-100, 
        100), asinh = list(0, 100)), seed = 3251991, verbose = FALSE)
    expect_equal(out2$T_name, "asinh")
    out3 <- rrscale(Y, trans_list = list(box_cox_negative = box_cox_negative), lims_list = list(box_cox_negative = c(-100, 
        100), asinh = list(0, 100)), seed = 3251991, verbose = FALSE)
    expect_equal(out3$T_name, "box_cox_negative")
    
    ## lims_list test
    out4 <- rrscale(Y, trans_list = list(box_cox_negative = box_cox_negative), lims_list = list(box_cox_negative = c(0, 
        100)), seed = 3251991, verbose = FALSE)
    expect_equal(out4$par_hat, 0, tolerance = 1e-15)
    
    ## z cutoff test
    out5 <- rrscale(Y, trans_list = list(box_cox_negative = box_cox_negative), lims_list = list(box_cox_negative = c(-100, 
        100)), z = 2, seed = 3251991, verbose = FALSE)
    expect_equal(out5$alg_control$z, 2)
    expect_gte(2, max(abs(out5$RR), na.rm = TRUE))
    expect_true(all(is.na(as.matrix(out5$O)[which(out5$Z > 2)])))
    out6 <- rrscale(Y, trans_list = list(box_cox_negative = box_cox_negative), lims_list = list(box_cox_negative = c(-100, 
        100)), z = 5, seed = 3251991, verbose = FALSE)
    expect_equal(out6$alg_control$z, 5)
    expect_true(all(is.na(as.matrix(out6$O)[which(out6$Z > 5)])))
    
    ## q cutoff
    out7 <- rrscale(Y, trans_list = list(box_cox_negative = box_cox_negative), lims_list = list(box_cox_negative = c(-100, 
        100)), q = 0.01, seed = 3251991, verbose = FALSE)
    Y_tmp <- as.matrix(Y)
    Yw <- winsor(Y_tmp, 0.01)
    mu <- mean(Yw, na.rm = TRUE)
    Ywc <- Yw - mu
    norm <- sqrt(mean(Ywc^2, na.rm = TRUE))
    Y_tmp <- (Y_tmp - mu)/norm
    Z <- Y_tmp
    expect_equal(as.matrix(out7$Z), Z, check.attributes = FALSE)
    
    out8 <- rrscale(Y, trans_list = list(box_cox_negative = box_cox_negative), lims_list = list(box_cox_negative = c(-100, 
        100)), q = 0, seed = 3251991, verbose = FALSE)
    Y_tmp <- as.matrix(Y)
    Yw <- winsor(Y_tmp, 0)
    mu <- mean(Yw, na.rm = TRUE)
    Ywc <- Yw - mu
    norm <- sqrt(mean(Ywc^2, na.rm = TRUE))
    Y_tmp <- (Y_tmp - mu)/norm
    Z <- Y_tmp
    expect_equal(as.matrix(out8$Z), Z, check.attributes = FALSE)
    
    ## zeros test
    Y2 <- Y
    Y2[1, 1] <- 0
    out10 <- rrscale(Y2, trans_list = list(box_cox_negative = box_cox_negative), 
        lims_list = list(box_cox_negative = c(-100, 100)), seed = 3251991, verbose = FALSE, 
        zeros = 1)
    expect_true(out10$NT[1, 1] == 1)
    expect_true(!is.na(out10$G[1, 1]))
    
    out11 <- rrscale(Y2, trans_list = list(box_cox_negative = box_cox_negative), 
        lims_list = list(box_cox_negative = c(-100, 100)), seed = 3251991, verbose = FALSE, 
        zeros = NA)
    expect_true(is.na(out11$NT[1, 1]))
    expect_true(is.na(out11$G[1, 1]))
    
    ## opt_control test
    opt_control <- DEoptim.control(trace = 2, reltol = 1e-10, itermax = 1000, steptol = 100)
    out12 <- rrscale(Y2, trans_list = list(box_cox_negative = box_cox_negative), 
        lims_list = list(box_cox_negative = c(-100, 100)), seed = 3251991, verbose = FALSE, 
        opt_control = opt_control, zeros = NA)
    expect_equal(out12$alg_control$opt_control, opt_control)
}
