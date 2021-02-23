subClosure <- function(closure, sub) {
    e <- environment(closure)
    restr <- new.env(parent = parent.env(e))
    for (ss in sub) {
        assign(ss, get(ss, envir = e), envir = restr)
    }
    environment(closure) <- restr
    return(closure)
}
