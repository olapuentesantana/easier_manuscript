### EstimAte Systems ImmunE Response (EaSIeR)

---
 
### How to install EaSIeR

´´´R
library(devtools)
install_github("olapuentesantana/easier_manuscript")
´´´

### Basic info on EaSIeR
For the essentials on how EaSIeR was developed, we recommend reading the following preprint [article](https://www.biorxiv.org/content/10.1101/2021.02.05.429977v1).

### Features of EaSIeR

# Computation of system-based signatures of the tumor microenvironment

´´´R
# Computation of cell fractions
cell_fractions <- compute_cell_fractions(RNA.tpm=tpm)
# Computation of pathway activity
pathways_activity <- compute_pathways_scores(RNA.countss=counts, remove.genes.ICB_proxies=TRUE)
# Computation of TF activity
tf_activity <- compute_TF_activity(RNA.tpm=tpm, remove.genes.ICB_proxies=FALSE)
# Computation of LR pairs weights
lrpairs_weights <- compute_LR_pairs(RNA.tpm=tpm, remove.genes.ICB_proxies=FALSE, compute.cytokines.pairs=FALSE, cancertype="pancan")
# Computation of Cell-Cell scores
ccpairs_scores <- compute_CC_pairs_grouped(lrpairs=lrpairs_weights$LRpairs, cancertype="pancan")
´´´

# Computation of different hallmarks of the immune response

´´´R
tasks <- c("CYT", "Roh_IS", "chemokines", "Davoli_IS", "IFNy", "Ayers_expIS", "Tcell_inflamed", "RIR", "TLS")
immune_response <- compute_gold_standards(RNA.tpm=tpm, list_gold_standards=tasks, cancertype=cancer_type, output_file_path=tmp_file_path)
´´´

# Predictions of patients' immune response
´´´R
predictions_immune_response <- predict_immune_response(pathways = pathways_activity,
                                                       immunecells = cell_fractions,
                                                       lrpairs = lrpairs_weights,
                                                       tfs = tf_activity,
                                                       ccpairs = ccpairs_scores,
                                                       cancertype = cancertype)
´´´