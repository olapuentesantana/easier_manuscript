# EstimAte Systems ImmunE Response (EaSIeR)

---
 
### How to install EaSIeR

```R
library(devtools)
```

First you will need to clone the repositroy and install it using `devtools::install()`. Then, load the namespace of `easier` package:
```R
library(easier)
```

### Basic info on EaSIeR
For the essentials on how EaSIeR was developed, we recommend reading the following preprint [article](https://www.biorxiv.org/content/10.1101/2021.02.05.429977v1).

### Features of EaSIeR

#### Computation of quantitative descriptors of the tumor microenvironment

| Quantitative Descriptors  | Descriptor conception  | Prior knowledge |
|------------------------------------------------------------------------ | ------------------ | ------------------------------ |
| Pathway activity  | Holland et al., BBAGRM, 2019; Schubert et al., Nat Commun, 2018 | Holland et al., BBAGRM, 2019; Schubert et al., Nat Commun, 2018 |
| Immune cell quantification  | Finotello et al., Genome Med, 2019 | Finotello et al., Genome Med, 2019|
| Transcription factor activity | Garcia-Alonso et al., Genome Res, 2019 | Garcia-Alonso et al., Genome Res, 2019 |
| Ligand-Receptor pairs | Lapuente-Santana et al., bioRxiv, 2021 | Ramilowski et al., Nat Commun, 2015 |
| Cell-cell interaction | Lapuente-Santana et al., bioRxiv, 2021 | Ramilowski et al.,  Nat Commun, 2015 |

```R
# Computation of cell fractions
cell_fractions <- compute_cell_fractions(RNA.tpm=tpm)
# Computation of pathway activity
pathways_activity <- compute_pathways_scores(RNA.counts=counts, remove.genes.ICB_proxies=TRUE)
# Computation of TF activity
tf_activity <- compute_TF_activity(RNA.tpm=tpm, remove.genes.ICB_proxies=FALSE)
# Computation of LR pairs weights
lrpairs_weights <- compute_LR_pairs(RNA.tpm=tpm, remove.genes.ICB_proxies=FALSE, compute.cytokines.pairs=FALSE, cancertype="pancan")
# Computation of Cell-Cell scores
ccpairs_scores <- compute_CC_pairs_grouped(lrpairs=lrpairs_weights$LRpairs, cancertype="pancan")
```

#### Computation of different hallmarks of the immune response

| Hallmark of the immune response | Original study |
|------------------------------------------- | ------------------ |
| Cytolytic activity (CYT) | Rooney et al, Cell, 2015 |
| Roh immune score (Roh_IS) | Roh et al., Sci. Transl. Med., 2017 |
| Chemokine signature (chemokines) | Messina et al., Nat. Sci. Rep., 2012 |
| Davoli immune signature (Davoli_IS) | Davoli et al., Science 2017 |
| IFNy signature (IFNy) | Ayers et al., JCI, 2017 |
| Expanded immune signature (Ayers_expIS) | Ayers et al., JCI, 2017 |
| T-cell inflamed signature (Tcell_inflamed) | Ayers et al., JCI, 2017 |
| Repressed immune resistance (RIR) | Jerby-Arnon et al., Cell, 2018 |
| Tertiary lymphoid structures signature (TLS) | Cabrita et al., Nature, 2020 |

```R
tasks <- c("CYT", "Roh_IS", "chemokines", "Davoli_IS", "IFNy", "Ayers_expIS", "Tcell_inflamed", "RIR", "TLS")
immune_response <- compute_gold_standards(RNA.tpm=tpm, list_gold_standards=tasks, cancertype=cancer_type, output_file_path=tmp_file_path)
```

#### Predictions of patients' immune response
```R
predictions_immune_response <- predict_immune_response(pathways = pathways_activity,
                                                       immunecells = cell_fractions,
                                                       lrpairs = lrpairs_weights,
                                                       tfs = tf_activity,
                                                       ccpairs = ccpairs_scores,
                                                       cancertype = cancertype)
```
