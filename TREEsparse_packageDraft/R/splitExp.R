#' Group the expression values of TFs and target gene according to the single cell cluster tree.
#'
#' @param y Expression vector of the target gene in all cells.
#' @param expr_tf_matrix Expression matrix of the TFs in all cells with rows correspond cells and columns corresponds to TFs.
#' @param cluster_groups A list of cell clusters.
#'
#' @return A list of expression values for the target gene and TFs in each cell cluster.
#'
#' @export
#'
#' @examples
#'
#' # Specific the samples for each single cell cluster.
#'
#' y  <-  runif(50)
#'
#' expr_tf_matrix <- matrix(runif(50*6), ncol = 6, byrow = F)
#'
#' cluster_groups <- list(sample(50, 20), sample(50, 20), sample(50, 20))
#'
#' tnl <- splitExp(y, expr_tf_matrix, cluster_groups)
#'
#' tam = matrix(0, nrow = 3, ncol = 3) # tree adjacency matrix with tam(i,j) = 1 if i is a parent of j.
#'
#' tam[1,2:3]=c(1,1)
#'
#' grns <- grnTree(tam, tnl, lambda1, lambda2)
#'
#' plot(grns)
#'
splitExp <- function(y, expr_tf_matrix, cluster_groups){

  if ( length(y) != nrow(expr_tf_matrix))
    stop("The number of samples for target gene in y and the TFs in expr_tf_matrix are not the same!")

  if (! is.list(cluster_groups))
    stop("The groups of cells in cluster_groups are not a list!")

  if( ! is.integer(unlist(cluster_groups)))
    stop("Wrong index in cluster_groups!")

  res <- vector(mode = "list", length = length(cluster_groups))

  for(ii in 1:length(cluster_groups)){

    res[[ii]] <- list(y = y[cluster_groups[[ii]]]
                      , X = expr_tf_matrix[cluster_groups[[ii]],])

  }

  return(res)

}
