# identify critical AA using ESM-variant

A variant can meet the PM1 criteria if it is located within a mutational
hotspot and/or a critical, well-established functional domain (e.g., the
active site of an enzyme) with no benign variations.

Not all functional domains in proteins have been annotated, and our
knowledge of protein hotspots is not complete. Therefore, a reliable
computational approach can provide additional assistance when
classifying a variant.

In this post, I utilized the ESM1-variant tool ([Nature Genetics,
2023](https://www.nature.com/articles/s41588-023-01465-0)) to identify
critical residues in the *RPE65* gene. 
The EMS-variant workflow employs ESM1b, a 650-million-parameter protein language model, to predict the potential effects of approximately 450 million missense variants in the
human genome. ESM1b demonstrated superior performance compared to
existing methods in classifying approximately 150,000 ClinVar/HGMD
missense variants as pathogenic or benign, as well as predicting
outcomes across 28 deep mutational scan datasets. Further information
can be found [here](https://www.nature.com/articles/s41588-023-01465-0).

ESM-variant's predictions are presented as log-likelihood ratios (LLR),
with a cutoff of -7.5. Employing this LLR threshold of -7.5 to
distinguish between pathogenic and benign variants resulted in an 81%
true-positive rate and an 82% true-negative rate in both ClinVar and
HGMD datasets.

To identify critically important residues, I established the following
criteria: (i) a residue must not possess an LLR score greater than -7.5,
and (ii) the mean LLR for all possible missense alterations must exceed
-7.5.

The provided text looks quite clear and informative. However, I've made a few minor adjustments for clarity and consistency:

The resulting data contains four columns:

`pos_variant`: Indicates the amino acid and its position.

`mean_score`: Represents the average Log-Likelihood ratio (LLR) score for all potential missense variations at a given position.

`Is_critical`: Based on established criteria, variants falling between -8.5 and -6.5 are labeled as 'in_gray_zone'.

`Is_in_guideline`: Indicates whether the position is identified as critical by the ClinGen LCA / eoRD VCEP guidelines.


Steps in this guide


### Steps in this guide

1- Obtaining the ESM-variant prediction from the model deployment in [huggingface](https://huggingface.co/spaces/ntranoslab/esm_variants)

2- Labeling the residues according to the abovementioned criteria. 



```R
library(dplyr)

# Reading the prediction score in R
gn<- read.csv("~/RPE65_(RPE65) _ Q16518.csv")

# Filter and calculate mean 
# Filter and calculate mean 
result_df <- gn %>%
  mutate(pos_variant = paste0(substr(variant, 1, 1), "_", pos)) %>%
  filter(!grepl("(.)\\d\\1", variant)) %>%  # Exclude rows with same letter before and after number
  group_by(pos, pos_variant) %>%
  summarise(mean_score = mean(score))
  
# Adding more info
result_df$is_critical <- cut(result_df$mean_score, 
                            breaks = c(-Inf,-8.5, -6.5, +Inf), 
                            labels = c("yes", "in_gray_zone", "no")) 
                            
# Currently there are some residues identified in the Variant Curation Expert Panel Specifications to the ACMG/AMP Variant Interpretation Guidelines 
# as critically important 
known_aa<- c("aa_180", "aa_182","aa_241","aa_313", "aa_417", "aa_527", "aa_107", "aa_125" )

# adding more data
result_df$Is_in_guideline <- ifelse(result_df$pos %in% known_aa, "yes",
"no")
RPE65_critical_residues <- result_df[, -1]
write.csv(RPE65_critical_residues, "RPE65_critical_residues.csv", row.names = F, col.names = T)

```
The resulting data is available in the file named "RPE65_critical_residues.csv".


