#' Calculate overall expression (OE)
#'
#' \code{get_OE_bulk} obtained from literature to calculate Immune resistance program  (Jerby-Arnon et al., 2018)
#'
#' @importFrom arules discretize
#'
#' @param r list
#' @param gene.sign string
#' @param num.rounds integer
#' @param full.flag boolean
#'
#' @return Random score
#'
#-------------------------------------------------------------------------------------------------------------
# function: calculate overall expression (OE)
get_OE_bulk <- function(r,gene.sign = NULL,num.rounds = 1000,full.flag = FALSE){
  set.seed(1234)
  r$genes.mean <- rowMeans(r$tpm)
  r$zscores <- sweep(r$tpm,1,r$genes.mean,FUN = '-')
  r$genes.dist <- r$genes.mean
  r$genes.dist.q <- arules::discretize(r$genes.dist,n.cat = 50)
  r$sig.scores <- matrix(data = 0,nrow = ncol(r$tpm),ncol = length(gene.sign))
  sig.names <- names(gene.sign)
  colnames(r$sig.scores) <- sig.names
  r$sig.scores.raw <- r$sig.scores
  rand.flag <- is.null(r$rand.scores)|!all(is.element(names(gene.sign),colnames(r$rand.scores)))
  if(rand.flag){
    print("Computing also random scores.")
    r$rand.scores <- r$sig.scores
  }
  for (i in sig.names){
    b.sign <- is.element(r$genes,gene.sign[[i]])
    if(sum(b.sign) < 2){next()}
    if(rand.flag){
      rand.scores <- get_semi_random_OE(r,r$genes.dist.q,b.sign,num.rounds = num.rounds)
    }else{
      rand.scores <- r$rand.scores[,i]
    }
    raw.scores <- colMeans(r$zscores[b.sign,])
    final.scores <- raw.scores-rand.scores
    r$sig.scores[,i] <- final.scores
    r$sig.scores.raw[,i] <- raw.scores
    r$rand.scores[,i] <- rand.scores
  }
  if(full.flag){return(r)}
  sig.scores <- r$sig.scores
  return(sig.scores)
}
