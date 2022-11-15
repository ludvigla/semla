#' @importFrom utils packageDescription
#' @importFrom glue glue
#'
.onAttach <- function (
    libname,
    pkgname
) {
  pkgVersion <- packageDescription(pkgname, fields = "Version")
  msg <- glue("{pkgname} v{pkgVersion}\n\n")

  citation <- glue("\n\nIf you use {pkgname} in published research, please cite the following paper:",
                     "\n\nL. Larsson, L. Franzen, 'Placeholder for title'")
  packageStartupMessage(paste0(strwrap(paste0(msg, citation)), collapse = "\n"))
}
