.centr <-
function(X,m) t (t(X) - m)

 .First.lib <-
 function(lib,pkg)
 {
    library.dynam("pcaPP",pkg,lib)
    library(mvtnorm)
    cat("pcaPP 0.1-1 loaded\n")
 }

