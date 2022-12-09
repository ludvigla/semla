#' Host a file server
#'
#' Hosts a file server from a specified directory. See
#' [beakr](https://github.com/MazamaScience/beakr) for more information.
#'
#' @param host A string with a valid IPv4 or IPv6 address to listen on,
#' Defaults to localhost "127.0.0.1"
#' @param port An integer that indicates the port to listen on.
#'
#' @import rlang
#' @import glue
#'
#' @return A beakr server with an active instance
#'
#' @examples
#' \donttest{
#'
#' # Host files in current directory
#' fs <- file_server(hostDir = "./")
#'
#' # Stop server
#' beakr::stopServer(fs)
#'
#' # Or stop all servers
#' beakr::stopAllServers()
#' }
#'
#' @export
file_server <- function (
  hostDir,
  host = "127.0.0.1",
  port = 8080
) {

  # Normalize path
  hostDir <- normalizePath(hostDir)

  # Check if directory exists
  if (!dir.exists(hostDir)) abort(glue("Directory '{hostDir}' does not exist"))

  # Check if beakr is installed
  if (!requireNamespace("beakr"))
    install.packages("beakr")

  # Create a new beakr instance
  beakr <- newBeakr()

  # Start a file server with CORS allowed
  fs <- beakr |>
    beakr::cors(origin = "*", headers = c("Origin", "X-Requested-With", "Content-Type", "Accept", "Range")) |>
    # Host the directory of static files
    beakr::serveStaticFiles(urlPath = "/", rootPath = hostDir, verbose = TRUE) |>
    # Start the server on port 25118
    beakr::listen(host = host, port = port, daemon = TRUE)

  return(fs)
}
