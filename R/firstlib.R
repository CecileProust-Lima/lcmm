############ First.lib ###############

.onLoad <- function(lib, pkg){
   library.dynam("lcmm", pkg, lib)
}

.onUnload <- function(libpath)
    library.dynam.unload("lcmm", libpath)

############ End of .First.lib ###############




