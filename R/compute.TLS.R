#' Compute tertiary lymphoid structures signature
#'
#' \code{compute_TLS} computes TLS signature as the geometric-mean of TLS signature genes
#' (Cabrita et al., 2020).
#'
#' @importFrom stats na.omit
#'
#' @param RNA.tpm numeric matrix with rows=genes and columns=samples
#'
#' @return numeric matrix with rows=samples and columns=TLS signature
#'
#' @export
#-------------------------------------------------------------------------------------------------------------

compute.TLS <- function(RNA.tpm){

  # Literature genes
  TLS.read <- c("CD79B", "CD1D", "CCR6", "LAT", "SKAP1", "CETP", "EIF1AY", "RBP5", "PTGDS")
  match_TLS.read <- match(TLS.read, rownames(RNA.tpm))

  if (anyNA(match_TLS.read)){
    warning(c("differenty named or missing signature genes : \n", paste(TLS.read[!TLS.read %in% rownames(RNA.tpm)], collapse = "\n")))
    match_TLS.read <- stats::na.omit(match_TLS.read)
  }

  # Subset RNA.tpm
  sub_gene.tpm <- RNA.tpm[match_TLS.read, ]

  # Calculation: geometric mean (so-called log-average) [TPM, 1 offset]
  geom_mean <- apply(sub_gene.tpm, 2, function(X) exp(mean(log2(X + 1))))

  message("TLS score computed")
  return(data.frame(TLS = geom_mean, check.names = FALSE))
}
