#' compute_CC_pairs_grouped
#'
#' \code{compute_CC_pairs_grouped} computes CC pairs (considering cell groups instead of individual cell types) from tpm,
#'  using the null model for CC interaction computed on TCGA data.
#'
#' @export
#'
#' @param lrpairs Ligand-leceptor pairs weights matrix
#' @param cancertype string character
#'
#' @return A list with the following elements:
#'         \describe{
#'               \item{score}{cell-cell interaction scores matrix with rows=samples and columns=cells}
#'               \item{pval}{corresponding p-values for each cell-cell pair, comparing the score with
#'               the null hypothesis}
#'         }
#'
#--------------------------------------------------------------------
# Compute CC pairs grouped from transcriptomics data.
#--------------------------------------------------------------------

compute_CC_pairs_grouped <- function(lrpairs, cancertype){

  # remove ligand receptor pairs that are always NA
  na_lrpairs <- apply(lrpairs, 2, function(x){all(is.na(x))})
  lrpairs <- lrpairs[, na_lrpairs == FALSE]

  # binarize the data: set a threshold to 10 TPM, only pairs where both ligand and receptor have
  # TPM > 10 are kept
  lrpairs_binary <- ifelse(lrpairs > log2(10+1), 1, 0)
  # sum(LR.pairs.binary)/(ncol(LR.pairs.binary)*nrow(LR.pairs.binary))*100 #percentage of

  # keep only the LR.pairs for which I have (non-zero) frequencies in the TCGA
  lrpairs_binary <- lrpairs_binary[, colnames(lrpairs_binary) %in% names(lr_frequency)]

  # cancer type specific network
  intercell_network <- intercell_network_cancer_spec[[cancertype]]

  # compute the CC score for each patient
  celltypes <- unique(c(as.character(intercell_network$cell1), as.character(intercell_network$cell2)))

  CC.pairs.score <- do.call(cbind, lapply(celltypes, function(celltype1){
    do.call(cbind, lapply(celltypes, function(celltype2){
      compute_CCpair_score(celltype1, celltype2, intercell_network,
                           lrpairs_binary, lr_frequency, compute_log=TRUE)
    }))
  }))

  metadata.CC.pairs <- do.call(rbind, lapply(celltypes, function(celltype1){
    do.call(rbind, lapply(celltypes, function(celltype2){
      data.frame(CCpair = gsub(" ", "", paste(celltype1, celltype2, sep="_")),
                 celltype1 = celltype1, celltype2 = celltype2)
    }))
  }))

  colnames(CC.pairs.score) <- metadata.CC.pairs$CCpair

  # CC.pairs.pval <- do.call(cbind, lapply(colnames(CC.pairs.score), function(CC.pair){
  #   sapply(CC.pairs.score[, CC.pair], function(x){
  #     compute_CCpair_pval(x, CC.pairs.score.random_grouped[, CC.pair])
  #   })
  # }))
  # colnames(CC.pairs.pval) <- colnames(CC.pairs.score)

  CC.pair <- list(score = as.data.frame(CC.pairs.score))

  return(CC.pair)
}
