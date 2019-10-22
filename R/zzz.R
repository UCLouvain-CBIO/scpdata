# Put onload and onAttach functionality

.onAttach <- function(libname, pkgname) {
  msg <- paste0("\nThis is scpdata version ",
                packageVersion("scpdata"), ".\n",
                "Use 'scpdata()' to list available data sets.")
  packageStartupMessage(msg)  
}