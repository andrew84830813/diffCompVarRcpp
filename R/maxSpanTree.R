#' Computes the maximum spanning tree (MST) of a set of log-ratios
#'
#'
#' Uses scoring from dcvScores() to construct and prune a log-ratio graph witht the maximum spanning tree method
#'
#' @importFrom magrittr %>%
#' @param featMatrix a n-sample by p-logratio matrix
#' @param dcvRanking  a p-row log-ratio scoring matrix usually from dcvScore()$lrs
#'
#' @examples
#' mstAll()
#'

#' @references
#' Hinton, A.L., Mucha, P.J., (2021). Simultaneous variable selection and group association testing in sparse high dimensional omics data. XXXX.
#'
#' @return A data.frame with the MST selected log-ratio features
#'
#' @export
mstAll <-
  function(featMatrix,dcvRanking){
    dcvRanking = dcvRanking[Ratio %in% colnames(featMatrix)]
    keyRats = tidyr::separate(dcvRanking,1,into = c("Num","Denom"),sep = "___",remove = F)
    el_= data.frame(keyRats$Num,keyRats$Denom,keyRats$Ratio)
    g = igraph::graph_from_edgelist(as.matrix(el_[,1:2]),directed = T)
    igraph::E(g)$weight = dcvRanking$rowmean
    g <-  igraph::minimum.spanning.tree(g,weights = -igraph::E(g)$weight)
    Ratios = data.frame(igraph::get.edgelist(g,names = T))
    Ratios_na = data.frame(paste(Ratios[,1],Ratios[,2],sep = "___"))
    el_ = cbind(Ratios,Ratios_na)
    colnames(el_)[1:3] = c("Num","Denom","Ratio")
    
    return(ratios = subset(featMatrix,select = el_$Ratio))
  }