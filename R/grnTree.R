#' Build GRNs for a tree of single cell clusters.
#'
#' @param tam Adjacency matrix for the tree of cell clusters with tam(i,j) = 1 if and only i is the parent of j.
#' @param tnl A list that contains the expression of target gene (y) and TFs (X) for each cell cluster.
#' @param lambda1 The first regularization parameter that imposed on sum_i{|b_i|_1}
#' @param lambda2 The second regularization parameter that imposed on sum_(i,j){|b_i-b_j|_1}, i is the parent of j
#'
#' @return The solution path for the all the coefficients |b_i| with i=1,2,...,K, where K is the number of single cell clusters.
#' @export
#'
#' @examples
#' # Input:
#'
#' #  tam: tree adjency matrix (0-1)^(K\timesK)
#' #           tam(i,j) = 1 if and only if j is a child of i
#' #  tnl: tree nodes list
#' #          node_i = list(y_i, X_i), i = 1,2,...,K
#'
#' tam = matrix(0, nrow = 3, ncol = 3)
#' tnl = vector(mode = "list", length = nrow(tam))
#' lambda1 = 1.0
#' lambda2 = 1.0
#'
#' # Specify the structure of the cell cluster tree
#' tam[1,2:3]=c(1,1)
#'
#' # Specify the background network for each single cell cluster.
#' B = matrix(0.0, ncol = 3, nrow = 6)
#' B[1,] = 1
#' B[2, 2] = 0.1
#' B[3, 3] = 0.1
#' print(B)
#'
#' for(ii in 1: length(tnl)){
#'
#'   X = matrix (runif(120), ncol = 6)
#'   b = B[,ii]
#'   y = X%*%b
#'   tnl[[ii]] = list(y = y, X= X)
#' }
#'
#' lambda2 = 0.0
#'
#' res <- grnTree(tam, tnl, lambda1, lambda2)
#'
#' plot(res, main = parse(text =  paste0('"grnTree"','~lambda[1] ==', lambda1,'*","~lambda[2] ==', lambda2)) )
#'
#' lambda2 = 1.0*0.5
#'
#' res <- grnTree(tam, tnl, lambda1, lambda2)
#'
#' plot(res, main = parse(text =  paste0('"grnTree"','~lambda[1] ==', lambda1,'*","~lambda[2] ==', lambda2)) )
#'
#'
grnTree <- function(tam, tnl, lambda1 = 1.0, lambda2 = 1.0){

  pt <- sapply(tnl, function(x) ncol(x$X))

  #Detect the number of candidate TFs in each node

  if(  length(as.vector(table(pt)))!=1 )
    stop("Number of TFs in the nodes are not the same!")

  nt <- sapply(tnl, function(x) nrow(x$X))

  Xt <- matrix(0.0, ncol = sum(pt), nrow = sum(nt))

  yt <- vector(mode = "numeric", length = sum(nt))


  Xt_colbreaks <- matrix(c(pt, cumsum(pt), c(1, cumsum(pt) + 1)[-length(pt)-1] )
                         , nrow = 3
                         , byrow = T)

  Xt_rowbreaks <- matrix(c(nt, cumsum(nt), c(1, cumsum(nt) + 1)[-length(nt)-1] )
                         , nrow = 3
                         , byrow = T)

  for(ii in 1: ncol(Xt_colbreaks)){

    Xt[Xt_rowbreaks[3, ii]: Xt_rowbreaks[2, ii], Xt_colbreaks[3, ii]: Xt_colbreaks[2, ii]] = tnl[[ii]]$X

    yt[Xt_rowbreaks[3, ii]: Xt_rowbreaks[2, ii]] = tnl[[ii]]$y

  }

  # Construct the multiplier matrix Dt
  if( lambda1 > 0 & lambda2 >0 )

    Dt <- kronecker( rbind( lambda1 * diag(nrow(tam)), lambda2 * t( (-tam+ diag(nrow(tam)))[,apply(tam,2, sum)>0]) )
                     , diag(pt[1]))

  else if(lambda1 > 0 & lambda2 ==0 )

    Dt <- kronecker(lambda1 * diag(nrow(tam)), diag(pt[1]))

  else if( lambda1 ==0 & lambda2 > 0)

    Dt <- kronecker(lambda2 * t( (-tam+ diag(nrow(tam)))[,apply(tam,2, sum)>0])
                    , diag(pt[1]))

  else if(lambda1 == 0 & lambda2 ==0)

    Dt = matrix(0.0, ncol = sum(pt), nrow = 1)

  library(genlasso)

  res <- genlasso(yt, Xt, Dt)

  return(res)

}
