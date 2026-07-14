.onAttach <- function(libname, pkgname) {
  version <- utils::packageVersion(pkgname)
  banner <- c(
    "",
    "    _    _          _ ____",
    "   / \\  | |__   ___| |  _ \\",
    "  / _ \\ | '_ \\ / _ \\ | |_) |",
    " / ___ \\| |_) |  __/ |  _ <",
    "/_/   \\_\\_.__/ \\___|_|_| \\_\\",
    "",
    paste0("  AbelR ", version, " loaded"),
    "  Bioinformatics analysis, modelling, and visualization",
    "  (\u2022\u0300\u1d17\u2022\u0301)\u0648  Ready to explore your data!",
    ""
  )
  packageStartupMessage(paste(banner, collapse = "\n"))
}
