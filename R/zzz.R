.onAttach <- function(libname, pkgname) {
  version <- utils::packageVersion(pkgname)
  packageStartupMessage(
    pkgname,
    " ",
    version,
    " loaded\n",
    "Reusable bioinformatics analysis and visualization functions."
  )
}
