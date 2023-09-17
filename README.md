# identify critical AA using ESM-variant

A variant can meet the PM1 criteria if it is located within a mutational hotspot and/or a critical, well-established functional domain (e.g., the active site of an enzyme) with no benign variations.

Not all functional domains in proteins have been annotated, and our knowledge of protein hotspots is not complete. Therefore, a reliable computational approach can provide additional assistance when classifying a variant.

In this post, I utilized the ESM1-variant tool ([Nature Genetics, 2023](https://www.nature.com/articles/s41588-023-01465-0)) to identify critical residues in the *RPE65* gene. The EMS-variant workflow employs ESM1b, a 650-million-parameter protein language model, to predict the potential effects of approximately 450 million missense variants in the human genome. ESM1b demonstrated superior performance compared to existing methods in classifying approximately 150,000 ClinVar/HGMD missense variants as pathogenic or benign, as well as predicting outcomes across 28 deep mutational scan datasets. Further information can be found [here](https://www.nature.com/articles/s41588-023-01465-0).

ESM-variant's predictions are presented as log-likelihood ratios (LLR), with a cutoff of -7.5. Employing this LLR threshold of -7.5 to distinguish between pathogenic and benign variants resulted in an 81% true-positive rate and an 82% true-negative rate in both ClinVar and HGMD datasets.

A critically important residue can be thought of as a residue where the substitution of the wild-type amino acid with almost any other amino acid, leads to damaging effects. To identify such residue, I established the following criteria: (i) Given the 19 possible LLR scores at each position, a residue must not possess an LLR score greater than -7.5, and (ii) the mean LLR for all possible missense alterations must exceed -7.5.


The resulting data contains four columns:

`pos_variant`: Indicates the amino acid and its position.

`mean_score`: Represents the average Log-Likelihood ratio (LLR) score for all potential missense variations at a given position.

`Is_critical`: If the mean score is > -6.5 , the residue is believed to be non-critical. If the mean score is <-8.5, it is believed to be critical for protein function. Residues with a mean score  between -8.5 and -6.5 are labeled as 'in_gray_zone'.

`Is_in_guideline`: Indicates whether the position is identified as critical by the ClinGen LCA / eoRD VCEP guidelines. As can be seen in the result, the analysis could correctly identify critical residues indicated in the guideline as critical. 

The resulting data looks like the following table:

| pos_variant | mean_score | Is_critical  | Is_in_guideline |
|-------------|------------|--------------|-----------------|
| M_1         | -9.13      | yes          | no              |
| S_2         | -5.79      | no           | no              |
| I_3         | -3.11      | no           | no              |
| Q_4         | -3.45      | no           | no              |
| V_5         | -3.32      | no           | no              |
| E_6         | -6.04      | no           | no              |
| H_7         | -6.44      | no           | no              |
| P_8         | -4.57      | no           | no              |
| A_9         | -4.5       | no           | no              |
| G_10        | -5.56      | no           | no              |
| G_11        | -8.31      | in_gray_zone | no              |
| Y_12        | -7.15      | in_gray_zone | no              |


### Analysis result is available for the following genes:

#### *RPE65* : [link](https://github.com/hamidghaedi/identify_critical_AA_using_ESM-variant/blob/main/RPE65_critical_residues.csv)
#### *CE290* : [link](https://github.com/hamidghaedi/identify_critical_AA_using_ESM-variant/blob/main/CE290_critical_residues.csv)
#### *GUCY2D* : [link](https://github.com/hamidghaedi/identify_critical_AA_using_ESM-variant/blob/main/GUCY2D_critical_residues.csv)






### Steps in this guide

1- Obtaining the ESM-variant prediction from the model deployment in [huggingface](https://huggingface.co/spaces/ntranoslab/esm_variants)

2- Labeling the residues according to the abovementioned criteria.

``` r
library(dplyr)

# Reading the prediction score in R
gn<- read.csv("~/RPE65_(RPE65) _ Q16518.csv")
# dropping wt_residue_wt substitution
gn <- gn[gn$score!=0,]

# Filter and calculate the mean 
result_df <- gn %>%
  mutate(pos_variant = paste0(substr(variant, 1, 1), "_", pos)) %>%
  group_by(pos, pos_variant) %>%
  summarise(mean_score = round(mean(score),2))
  
# Adding more info
result_df$Is_critical <- cut(result_df$mean_score, 
                            breaks = c(-Inf,-8.5, -6.5, +Inf), 
                            labels = c("yes", "in_gray_zone", "no")) 
                            
# Currently there are some residues identified in the Variant Curation Expert Panel Specifications to the ACMG/AMP Variant Interpretation Guidelines 
# as critically important 
known_aa<- c("H_180", "H_182","His_241","His_313", "E_417", "H_527", "A_107", "G_125" )

# Adding more data
result_df$Is_in_guideline <- ifelse(result_df$pos_variant %in% known_aa, "yes",
"no")
RPE65_critical_residues <- result_df[, -1]
write.csv(RPE65_critical_residues, "RPE65_critical_residues.csv", row.names = F)
```
