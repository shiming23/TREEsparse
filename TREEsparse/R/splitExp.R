splitExp <- function(y, expr_tf_matrix, cluster_groups){

  if ( length(y) != nrow(expr_tf_matrix))
    stop("The number of samples for target gene in y and the TFs in expr_tf_matrix are not the same!")

  if (! is.list(cluster_groups))
    stop("The groups of cells in cluster_groups are not a list!")

  if( ! is.numeric(unlist(cluster_groups)))
    stop("Wrong index in cluster_groups!")

  res <- vector(mode = "list", length = length(cluster_groups))

  for(ii in 1:length(cluster_groups)){

    res[[ii]] <- list(y = y[cluster_groups[[ii]]]
                      , X = expr_tf_matrix[cluster_groups[[ii]],])

  }

  return(res)

}
