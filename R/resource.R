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
#' download again. This will overwrite
#' pre-existing database. Default is True.
#' @return Returns `NULL`, invisibly. 
#' @export
#' @examples
#' createResources()
createResources <- function() {

    return(invisible(NULL))
}


