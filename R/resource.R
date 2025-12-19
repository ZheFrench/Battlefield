####################################################
###     Create / Get Ressources                  ###
####################################################

#' Create all resources.
#'
#' Create cache for all resources (pathways, or PWC network)
#' downloaded from the web when library is first loaded.
#' This part is handled with BiocFileCache.
#' Otherwise datatabase, is handled by another process
#' not relying on BiocFileCache instance.
#'
#' @param onRequest logical True if you force
#' download again. This will overwrite
#' pre-existing database. Default is True.
#' @param verbose Default is FALSE
#' @return Returns `NULL`, invisibly. 
#' @importFrom curl has_internet
#' @importFrom cli cli_alert_danger cli_alert
#' @export
#' @examples
#' createResources(onRequest=FALSE)
createResources <- function(onRequest = TRUE, verbose = FALSE) {
    cacheDir <- .SignalR$BulkSignalR_CACHEDIR
    resourcesCacheDir <- paste(cacheDir, "resources", sep = "/")

    hasInternet <- tryCatch(expr={curl::has_internet()}, 
        error = FALSE)
    
    if (!hasInternet & 
    !file.exists(resourcesCacheDir)) {
        cli::cli_alert_danger("Your internet connection is off :")
        stop(
        "- Remote resources can't be downloaded.\n"
        )   
    }

    if (!hasInternet & 
        onRequest) {
        cli::cli_alert_danger("Your internet connection is off :")
        stop(
        "- Remote resources can't be downloaded.\n"
        )   
    }

    # Do it once, onLoad
    if (!dir.exists(resourcesCacheDir) | onRequest) {
        .cacheAdd(fpath = .SignalR$BulkSignalR_GO_URL,
            cacheDir = resourcesCacheDir,
            resourceName = "GO-BP", 
            verbose = verbose, download = TRUE)
        .cacheAdd(fpath = .SignalR$BulkSignalR_Reactome_URL,
            cacheDir = resourcesCacheDir, 
            resourceName = "Reactome",
            verbose = verbose, download = TRUE)
        .cacheAdd(fpath = .SignalR$BulkSignalR_Network_URL,
            cacheDir = resourcesCacheDir, 
            resourceName = "Network",
            verbose = verbose, download = TRUE)

    }
    cacheVersion(dir="resources")
    
    return(invisible(NULL))
}


