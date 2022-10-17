#' @importFrom utils packageDescription
#' @importFrom glue glue
#' @importFrom pillar pillar
.onAttach <- function (
    libname,
    pkgname
) {
  pkgVersion <- packageDescription(pkgname, fields = "Version")
  msg <- glue("{pkgname} v{pkgVersion}\n\n")

  citation <- glue("\n\nIf you use {pkgname} in published research, please cite the following paper:",
                     "\n\nL. Larsson, L. Franzén, 'Placeholder for title'")
  packageStartupMessage(paste0(strwrap(pillar::style_subtle(paste0(msg, citation))), collapse = "\n"))
}
