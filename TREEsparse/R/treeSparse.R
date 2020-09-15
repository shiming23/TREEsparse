#' Fast computation for tree-structured sparsity
#'
#' @param tam tree matrix
#' @param tnl tree node list
#' @param lambda1 turning parameter alpha 1
#' @param lambda2 turning parameter alpha 2
#'
#' @return coefficient for the tree-structured sparsity
#' @export
#'
#' @examples
#'
treeSparse <- function(tam, tnl, lambda1 = 1.0 , lambda2 = 0.0){

  if(!nrow(tam)==ncol(tam))
    stop("Error: The tree adjancy matrix must be a square 0-1 matrix!")

  if(sum(abs(diag(tam))))
    stop("Error: The tree must NOT contain self-loops!")

  if(!length(tnl)==nrow(tam))
    stop("Error: The number of nodes in the tree-nodes-lists and the tree-adjacency-matrix are not equal!")

  if(!(prod(sapply(tnl, function(node) is.vector(node$y))) &
       prod(sapply(tnl, function(node) is.matrix(node$X)))))
    stop("Error: The predictors and responses must be matrix and vector, respectively!")

  if(! length( unique(sapply(tnl, function(node) ncol(node$X))) )==1)
    stop("Error: The predictor variables in each node are not equal!")

  if(! ((lambda1>0)&(lambda2>=0)) )
    stop("The value of lambda1 should be positive and lambda2 should be non-negative!")

  sample_number <- sapply(tnl, simplify = T, function(node) length(node$y))

  predictor_number <- sapply(tnl, simplify = T, function(node) ncol(node$X))

  y_stack <- as.vector(unlist( sapply(tnl, simplify = T, function(node) node$y) ))

  X_stack <- matrix(0.0, ncol = sum(predictor_number), nrow = sum(sample_number))

  cum_sample_number <- c( 0, cumsum(sample_number)) + 1

  cum_predictor_number <- c(0, cumsum(predictor_number))+ 1

  for( ii in 1:length(tnl)){

    X_stack[ cum_sample_number[ii]: (cum_sample_number[ii+1]-1)
              , cum_predictor_number[ii]: (cum_predictor_number[ii+1]-1) ] = tnl[[ii]]$X

  }

  D1 <- diag(x = 1, nrow = length(predictor_number))

  edges_list <- which(tam!=0, arr.ind = T)

  D2 <- matrix(0.0, nrow = nrow(edges_list), ncol = length(predictor_number))

  for(ii in 1:nrow(edges_list)){

    D2[ii, edges_list[ii,1]] = 1

    D2[ii, edges_list[ii,2]] = -1

  }

  D <- rbind(D1, D2)

  D_kronecker <- kronecker( D , diag(x = 1, nrow = predictor_number[1]) )

  D_kroInverse <- kronecker( solve(t(D)%*%D)%*%t(D), diag(x = 1, nrow = predictor_number[1]) )

  library(lars)

  res <- lars(x = X_stack%*%D_kroInverse, y = y_stack, max.steps = length(tnl)
              , type = "lar"
              , normalize = F
              , intercept = F
              , use.Gram = F)

#  return(D_krooInverse %*% coef(res)[length(tnl)+1,])

  return(D_kroInverse %*% t(coef(res)))

#  return(res)

}
