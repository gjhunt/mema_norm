
pp <- function(pl, fn, ...) {
    dir.create(dirname(fn), showWarnings = FALSE, recursive = TRUE)
    if ("ggplot" %in% class(pl))
        pl <- list(pl)
    pb <- progress::progress_bar$new(total = length(pl), width = 120, format = paste0("Saving ",
        fn, " [:bar] :current/:total (:percent) elapsed: :elapsed eta: :eta"))
    grDevices::pdf(file = fn, ...)
    for (p in pl) {
        print(p)
        pb$tick()
    }
    grDevices::dev.off()
}


`%+%` <- function(a, b) {
    paste0(a, b)
}

lapplyp <- function(L, FUN, ..., verbose = TRUE, .name = "Applying") {
    if (!verbose)
        return(lapply(L, FUN, ...))
    ii <- 1:length(L)
    pb <- progress::progress_bar$new(total = length(L), width = 120, format = paste0(.name,
        " [:bar] :current/:total (:percent) elapsed: :elapsed eta: :eta"))
    ret <- lapply(L, function(l) {
        pb$tick()
        FUN(l, ...)
    })
    return(ret)
}
