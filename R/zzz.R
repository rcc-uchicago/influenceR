# Load this package's dynamic library.
.onLoad <- function(libname, pkgname){
      library.dynam("influenceR", package=pkgname, lib.loc=libname)
}
