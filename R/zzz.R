#' @importFrom utils packageDescription
#' @importFrom glue glue
#'
.onAttach <- function (
    libname,
    pkgname
) {
  pkgVersion <- packageDescription(pkgname, fields = "Version")
  msg <- glue("{pkgname} v{pkgVersion}\n\n")

  citation <- glue("\n\nAuthors:",
                     "\nLudvig Larsson and Lovisa Franzen")
  packageStartupMessage(paste0(strwrap(paste0(msg, citation)), collapse = "\n"))
}
