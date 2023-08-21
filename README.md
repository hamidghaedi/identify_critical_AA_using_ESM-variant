# identify critical AA using ESM-variant

A variant can meet the PM1 criteria if it is located within a mutational hotspot and/or a critical, well-established functional domain (e.g., the active site of an enzyme) with no benign variations.

Not all functional domains in proteins have been annotated, and our knowledge of protein hotspots is not complete. Therefore, a reliable computational approach can provide additional assistance when classifying a variant.

In this post, I utilized the ESM1-variant tool ([Nature Genetics, 2023](https://www.nature.com/articles/s41588-023-01465-0)) to identify critical residues in the *RPE65* gene. The EMS-variant workflow employs ESM1b, a 650-million-parameter protein language model, to predict the potential effects of approximately 450 million missense variants in the human genome. ESM1b demonstrated superior performance compared to existing methods in classifying approximately 150,000 ClinVar/HGMD missense variants as pathogenic or benign, as well as predicting outcomes across 28 deep mutational scan datasets. Further information can be found [here](https://www.nature.com/articles/s41588-023-01465-0).

ESM-variant's predictions are presented as log-likelihood ratios (LLR), with a cutoff of -7.5. Employing this LLR threshold of -7.5 to distinguish between pathogenic and benign variants resulted in an 81% true-positive rate and an 82% true-negative rate in both ClinVar and HGMD datasets.

To identify critically important residues, I established the following criteria: (i) a residue must not possess an LLR score greater than -7.5, and (ii) the mean LLR for all possible missense alterations must exceed -7.5.

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

``` r
library(dplyr)

# Reading the prediction score in R
gn<- read.csv("~/RPE65_(RPE65) _ Q16518.csv")

# Filter and calculate mean 
# Filter and calculate mean 
result_df <- gn %>%
  mutate(pos_variant = paste0(substr(variant, 1, 1), "_", pos)) %>%
  filter(!grepl("(.)\\d\\1", variant)) %>%  # Exclude rows with same letter before and after number
  group_by(pos, pos_variant) %>%
  summarise(mean_score = round(mean(score),2)
  
# Adding more info
result_df$Is_critical <- cut(result_df$mean_score, 
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

The resulting data is available belwo and also fo download as the file named "RPE65_critical_residues.csv" from this repo.

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
| G_10        | -5.28      | no           | no              |
| G_11        | -7.9       | in_gray_zone | no              |
| Y_12        | -6.8       | in_gray_zone | no              |
| K_13        | -5.87      | no           | no              |
| K_14        | -6.4       | no           | no              |
| L_15        | -11.83     | yes          | no              |
| F_16        | -10.14     | yes          | no              |
| E_17        | -7.28      | in_gray_zone | no              |
| T_18        | -8.38      | in_gray_zone | no              |
| V_19        | -7.75      | in_gray_zone | no              |
| E_20        | -7.81      | in_gray_zone | no              |
| E_21        | -11.9      | yes          | no              |
| L_22        | -6.75      | in_gray_zone | no              |
| S_23        | -5.02      | no           | no              |
| S_24        | -4.91      | no           | no              |
| P_25        | -11.13     | yes          | no              |
| L_26        | -8.81      | yes          | no              |
| T_27        | -4.94      | no           | no              |
| A_28        | -11.19     | yes          | no              |
| H_29        | -1.74      | no           | no              |
| V_30        | -10.17     | yes          | no              |
| T_31        | -5.79      | no           | no              |
| G_32        | -12.45     | yes          | no              |
| R_33        | -5.28      | no           | no              |
| I_34        | -11.38     | yes          | no              |
| P_35        | -13.64     | yes          | no              |
| L_36        | -2.96      | no           | no              |
| W_37        | -13.19     | yes          | no              |
| L_38        | -12.93     | yes          | no              |
| T_39        | -5.26      | no           | no              |
| G_40        | -14.52     | yes          | no              |
| S_41        | -7.65      | in_gray_zone | no              |
| L_42        | -11.02     | yes          | no              |
| L_43        | -11.72     | yes          | no              |
| R_44        | -13.53     | yes          | no              |
| C_45        | -7.97      | in_gray_zone | no              |
| G_46        | -13.94     | yes          | no              |
| P_47        | -12.65     | yes          | no              |
| G_48        | -12.76     | yes          | no              |
| L_49        | -9.49      | yes          | no              |
| F_50        | -10.6      | yes          | no              |
| E_51        | -10.44     | yes          | no              |
| V_52        | -10.67     | yes          | no              |
| G_53        | -11.7      | yes          | no              |
| S_54        | -5.83      | no           | no              |
| E_55        | -10.25     | yes          | no              |
| P_56        | -7.35      | in_gray_zone | no              |
| F_57        | -10.79     | yes          | no              |
| Y_58        | -5.72      | no           | no              |
| H_59        | -12.54     | yes          | no              |
| L_60        | -10.05     | yes          | no              |
| F_61        | -13.49     | yes          | no              |
| D_62        | -14.2      | yes          | no              |
| G_63        | -12.49     | yes          | no              |
| Q_64        | -8.52      | yes          | no              |
| A_65        | -12.57     | yes          | no              |
| L_66        | -12.18     | yes          | no              |
| L_67        | -13.17     | yes          | no              |
| H_68        | -9.87      | yes          | no              |
| K_69        | -10.82     | yes          | no              |
| F_70        | -14.75     | yes          | no              |
| D_71        | -6.95      | in_gray_zone | no              |
| F_72        | -9.43      | yes          | no              |
| K_73        | -8.22      | in_gray_zone | no              |
| E_74        | -8.33      | in_gray_zone | no              |
| G_75        | -12.01     | yes          | no              |
| H_76        | -8.3       | in_gray_zone | no              |
| V_77        | -12.84     | yes          | no              |
| T_78        | -9.59      | yes          | no              |
| Y_79        | -13.46     | yes          | no              |
| H_80        | -5.61      | no           | no              |
| R_81        | -7.27      | in_gray_zone | no              |
| R_82        | -10.75     | yes          | no              |
| F_83        | -11.88     | yes          | no              |
| I_84        | -7.98      | in_gray_zone | no              |
| R_85        | -6.86      | in_gray_zone | no              |
| T_86        | -8.99      | yes          | no              |
| D_87        | -9.84      | yes          | no              |
| A_88        | -8.92      | yes          | no              |
| Y_89        | -13.38     | yes          | no              |
| V_90        | -5.86      | no           | no              |
| R_91        | -8.34      | in_gray_zone | no              |
| A_92        | -7.52      | in_gray_zone | no              |
| M_93        | -8.12      | in_gray_zone | no              |
| T_94        | -4.99      | no           | no              |
| E_95        | -7.3       | in_gray_zone | no              |
| K_96        | -8         | in_gray_zone | no              |
| R_97        | -11.15     | yes          | no              |
| I_98        | -10.82     | yes          | no              |
| V_99        | -12.85     | yes          | no              |
| I_100       | -8.61      | yes          | no              |
| E_102       | -11.47     | yes          | no              |
| F_103       | -11.06     | yes          | no              |
| G_104       | -12.72     | yes          | no              |
| T_105       | -12.29     | yes          | no              |
| C_106       | -3.09      | no           | no              |
| A_107       | -10.68     | yes          | no              |
| F_108       | -8.51      | yes          | no              |
| P_109       | -12.13     | yes          | no              |
| D_110       | -13.92     | yes          | no              |
| C_112       | -14.44     | yes          | no              |
| K_113       | -10.22     | yes          | no              |
| N_114       | -11.88     | yes          | no              |
| I_115       | -13.19     | yes          | no              |
| F_116       | -14.75     | yes          | no              |
| S_117       | -9.29      | yes          | no              |
| R_118       | -10.58     | yes          | no              |
| F_119       | -13.78     | yes          | no              |
| F_120       | -10.24     | yes          | no              |
| Y_122       | -9.39      | yes          | no              |
| F_123       | -11.5      | yes          | no              |
| R_124       | -7.76      | in_gray_zone | no              |
| G_125       | -9.99      | yes          | no              |
| V_126       | -8.2       | in_gray_zone | no              |
| E_127       | -9.13      | yes          | no              |
| V_128       | -5.82      | no           | no              |
| T_129       | -11.68     | yes          | no              |
| D_130       | -15        | yes          | no              |
| A_132       | -9.51      | yes          | no              |
| L_133       | -9.15      | yes          | no              |
| V_134       | -11.79     | yes          | no              |
| N_135       | -11.09     | yes          | no              |
| V_136       | -10.38     | yes          | no              |
| Y_137       | -9.99      | yes          | no              |
| P_138       | -7.69      | in_gray_zone | no              |
| V_139       | -9.52      | yes          | no              |
| G_140       | -8.37      | in_gray_zone | no              |
| D_142       | -10.03     | yes          | no              |
| Y_143       | -10.72     | yes          | no              |
| Y_144       | -13.83     | yes          | no              |
| A_145       | -10.49     | yes          | no              |
| C_146       | -7.5       | in_gray_zone | no              |
| T_147       | -10.68     | yes          | no              |
| E_148       | -14.33     | yes          | no              |
| T_149       | -11.11     | yes          | no              |
| N_150       | -9.89      | yes          | no              |
| I_152       | -12.83     | yes          | no              |
| T_153       | -8.81      | yes          | no              |
| K_154       | -10.68     | yes          | no              |
| I_155       | -12.31     | yes          | no              |
| N_156       | -10.01     | yes          | no              |
| P_157       | -9.01      | yes          | no              |
| E_158       | -7.41      | in_gray_zone | no              |
| T_159       | -9.98      | yes          | no              |
| L_160       | -12.74     | yes          | no              |
| T_162       | -10.58     | yes          | no              |
| I_163       | -6.91      | in_gray_zone | no              |
| K_164       | -9.89      | yes          | no              |
| Q_165       | -7.92      | in_gray_zone | no              |
| V_166       | -10.31     | yes          | no              |
| D_167       | -8.69      | yes          | no              |
| L_168       | -7.47      | in_gray_zone | no              |
| C_169       | -5.74      | no           | no              |
| N_170       | -7.48      | in_gray_zone | no              |
| V_172       | -11.3      | yes          | no              |
| S_173       | -8.19      | in_gray_zone | no              |
| V_174       | -12.09     | yes          | no              |
| N_175       | -10.77     | yes          | no              |
| G_176       | -8.04      | in_gray_zone | no              |
| A_177       | -8.59      | yes          | no              |
| T_178       | -11.51     | yes          | no              |
| A_179       | -11.55     | yes          | no              |
| H_180       | -14.75     | yes          | no              |
| H_182       | -12.34     | yes          | no              |
| I_183       | -6.57      | in_gray_zone | no              |
| E_184       | -7.37      | in_gray_zone | no              |
| N_185       | -4         | no           | no              |
| D_186       | -12.53     | yes          | no              |
| G_187       | -11.42     | yes          | no              |
| T_188       | -9.45      | yes          | no              |
| V_189       | -8.39      | in_gray_zone | no              |
| Y_190       | -10.75     | yes          | no              |
| I_192       | -9.46      | yes          | no              |
| G_193       | -11.96     | yes          | no              |
| N_194       | -9.63      | yes          | no              |
| C_195       | -8.59      | yes          | no              |
| F_196       | -10.7      | yes          | no              |
| G_197       | -9.12      | yes          | no              |
| K_198       | -7.48      | in_gray_zone | no              |
| N_199       | -5.34      | no           | no              |
| F_200       | -4.88      | no           | no              |
| S_201       | -4.12      | no           | no              |
| A_203       | -3.16      | no           | no              |
| Y_204       | -12.61     | yes          | no              |
| N_205       | -7.99      | in_gray_zone | no              |
| I_206       | -11.16     | yes          | no              |
| V_207       | -10.55     | yes          | no              |
| K_208       | -8.41      | in_gray_zone | no              |
| I_209       | -9.37      | yes          | no              |
| P_210       | -9.64      | yes          | no              |
| P_211       | -9.57      | yes          | no              |
| Q_213       | -4.04      | no           | no              |
| A_214       | -4.49      | no           | no              |
| D_215       | -6.26      | no           | no              |
| K_216       | -4.22      | no           | no              |
| E_217       | -5.59      | no           | no              |
| D_218       | -6.84      | in_gray_zone | no              |
| P_219       | -6.76      | in_gray_zone | no              |
| I_220       | -4.83      | no           | no              |
| S_221       | -4.76      | no           | no              |
| S_223       | -4.86      | no           | no              |
| E_224       | -8.21      | in_gray_zone | no              |
| I_225       | -8.09      | in_gray_zone | no              |
| V_226       | -8.08      | in_gray_zone | no              |
| V_227       | -5.62      | no           | no              |
| Q_228       | -4.49      | no           | no              |
| F_229       | -7.26      | in_gray_zone | no              |
| P_230       | -10.72     | yes          | no              |
| C_231       | -5.69      | no           | no              |
| D_233       | -5.93      | no           | no              |
| R_234       | -6.32      | no           | no              |
| F_235       | -6.39      | no           | no              |
| K_236       | -6.24      | no           | no              |
| P_237       | -9.76      | yes          | no              |
| S_238       | -9.67      | yes          | no              |
| Y_239       | -13.25     | yes          | no              |
| V_240       | -6.68      | in_gray_zone | no              |
| H_241       | -13.29     | yes          | no              |
| F_243       | -13.51     | yes          | no              |
| G_244       | -10.36     | yes          | no              |
| L_245       | -9.19      | yes          | no              |
| T_246       | -13.46     | yes          | no              |
| P_247       | -6.86      | in_gray_zone | no              |
| N_248       | -12.55     | yes          | no              |
| Y_249       | -12.29     | yes          | no              |
| I_250       | -9.95      | yes          | no              |
| V_251       | -12.25     | yes          | no              |
| V_253       | -10.37     | yes          | no              |
| E_254       | -12.23     | yes          | no              |
| T_255       | -4.41      | no           | no              |
| P_256       | -12        | yes          | no              |
| V_257       | -6.71      | in_gray_zone | no              |
| K_258       | -10.44     | yes          | no              |
| I_259       | -9.1       | yes          | no              |
| N_260       | -9.97      | yes          | no              |
| L_261       | -9.44      | yes          | no              |
| K_263       | -10.04     | yes          | no              |
| F_264       | -7.16      | in_gray_zone | no              |
| L_265       | -8.58      | yes          | no              |
| S_266       | -5.59      | no           | no              |
| S_267       | -7.32      | in_gray_zone | no              |
| W_268       | -9.71      | yes          | no              |
| S_269       | -3.49      | no           | no              |
| L_270       | -8.56      | yes          | no              |
| W_271       | -5.63      | no           | no              |
| A_273       | -7.31      | in_gray_zone | no              |
| N_274       | -7.53      | in_gray_zone | no              |
| Y_275       | -9.39      | yes          | no              |
| M_276       | -7.54      | in_gray_zone | no              |
| D_277       | -10.76     | yes          | no              |
| C_278       | -11.01     | yes          | no              |
| F_279       | -10.56     | yes          | no              |
| E_280       | -8.79      | yes          | no              |
| S_281       | -4.43      | no           | no              |
| E_283       | -8.81      | yes          | no              |
| T_284       | -5.24      | no           | no              |
| M_285       | -3.9       | no           | no              |
| G_286       | -9.54      | yes          | no              |
| V_287       | -7.04      | in_gray_zone | no              |
| W_288       | -7.5       | in_gray_zone | no              |
| L_289       | -9.22      | yes          | no              |
| H_290       | -7.88      | in_gray_zone | no              |
| I_291       | -9.15      | yes          | no              |
| D_293       | -10.66     | yes          | no              |
| K_294       | -10.77     | yes          | no              |
| K_295       | -6.8       | in_gray_zone | no              |
| R_296       | -4.19      | no           | no              |
| K_297       | -6.56      | in_gray_zone | no              |
| K_298       | -4.82      | no           | no              |
| Y_299       | -2         | no           | no              |
| L_300       | -7.53      | in_gray_zone | no              |
| N_301       | -3.17      | no           | no              |
| N_302       | -1.1       | no           | no              |
| Y_304       | -9.36      | yes          | no              |
| R_305       | -3.93      | no           | no              |
| T_306       | -6.78      | in_gray_zone | no              |
| S_307       | -4.86      | no           | no              |
| P_308       | -8.58      | yes          | no              |
| F_309       | -10.22     | yes          | no              |
| N_310       | -8.13      | in_gray_zone | no              |
| L_311       | -7.24      | in_gray_zone | no              |
| F_312       | -11.6      | yes          | no              |
| H_314       | -9.89      | yes          | no              |
| I_315       | -13.02     | yes          | no              |
| N_316       | -12.72     | yes          | no              |
| T_317       | -10.5      | yes          | no              |
| Y_318       | -12.32     | yes          | no              |
| E_319       | -13.83     | yes          | no              |
| D_320       | -9.97      | yes          | no              |
| N_321       | -8.06      | in_gray_zone | no              |
| G_322       | -10.49     | yes          | no              |
| L_324       | -10.5      | yes          | no              |
| I_325       | -10.37     | yes          | no              |
| V_326       | -10.24     | yes          | no              |
| D_327       | -14.36     | yes          | no              |
| L_328       | -8.75      | yes          | no              |
| C_329       | -12.53     | yes          | no              |
| C_330       | -9.54      | yes          | no              |
| W_331       | -7.2       | in_gray_zone | no              |
| K_332       | -8.47      | in_gray_zone | no              |
| F_334       | -2.46      | no           | no              |
| E_335       | -7.71      | in_gray_zone | no              |
| F_336       | -6.28      | no           | no              |
| V_337       | -7.3       | in_gray_zone | no              |
| Y_338       | -6.49      | no           | no              |
| N_339       | -7.46      | in_gray_zone | no              |
| Y_340       | -7.15      | in_gray_zone | no              |
| L_341       | -9.18      | yes          | no              |
| Y_342       | -10.35     | yes          | no              |
| A_344       | -3.52      | no           | no              |
| N_345       | -10.13     | yes          | no              |
| L_346       | -11.58     | yes          | no              |
| R_347       | -7.88      | in_gray_zone | no              |
| E_348       | -6.92      | in_gray_zone | no              |
| N_349       | -8.36      | in_gray_zone | no              |
| W_350       | -6.88      | in_gray_zone | no              |
| E_351       | -7.9       | in_gray_zone | no              |
| E_352       | -6.64      | in_gray_zone | no              |
| K_354       | -6.7       | in_gray_zone | no              |
| K_355       | -6         | no           | no              |
| N_356       | -7.14      | in_gray_zone | no              |
| A_357       | -4.06      | no           | no              |
| R_358       | -4.82      | no           | no              |
| K_359       | -3.77      | no           | no              |
| A_360       | -3.39      | no           | no              |
| P_361       | -7.02      | in_gray_zone | no              |
| Q_362       | -4.21      | no           | no              |
| E_364       | -6.53      | in_gray_zone | no              |
| V_365       | -6.68      | in_gray_zone | no              |
| R_366       | -8.4       | in_gray_zone | no              |
| R_367       | -14.1      | yes          | no              |
| Y_368       | -11.56     | yes          | no              |
| V_369       | -12.26     | yes          | no              |
| L_370       | -11.85     | yes          | no              |
| P_371       | -12.48     | yes          | no              |
| L_372       | -11.28     | yes          | no              |
| I_374       | -8.6       | yes          | no              |
| D_375       | -6.42      | no           | no              |
| K_376       | -6.34      | no           | no              |
| A_377       | -5.3       | no           | no              |
| D_378       | -4.84      | no           | no              |
| T_379       | -5         | no           | no              |
| G_380       | -9.15      | yes          | no              |
| K_381       | -4.94      | no           | no              |
| N_382       | -11.74     | yes          | no              |
| V_384       | -10.38     | yes          | no              |
| T_385       | -6.03      | no           | no              |
| L_386       | -12.83     | yes          | no              |
| P_387       | -9.74      | yes          | no              |
| N_388       | -6.76      | in_gray_zone | no              |
| T_389       | -7.52      | in_gray_zone | no              |
| T_390       | -5.58      | no           | no              |
| A_391       | -9.83      | yes          | no              |
| T_392       | -7.83      | in_gray_zone | no              |
| I_394       | -6.9       | in_gray_zone | no              |
| L_395       | -3.18      | no           | no              |
| C_396       | -4.31      | no           | no              |
| S_397       | -4.79      | no           | no              |
| D_398       | -9.23      | yes          | no              |
| E_399       | -5.85      | no           | no              |
| T_400       | -6.66      | in_gray_zone | no              |
| I_401       | -8.89      | yes          | no              |
| W_402       | -9.11      | yes          | no              |
| L_403       | -8.48      | in_gray_zone | no              |
| P_405       | -9.06      | yes          | no              |
| E_406       | -11.04     | yes          | no              |
| V_407       | -5.63      | no           | no              |
| L_408       | -9.72      | yes          | no              |
| F_409       | -7.56      | in_gray_zone | no              |
| S_410       | -5.02      | no           | no              |
| G_411       | -8.57      | yes          | no              |
| P_412       | -8.01      | in_gray_zone | no              |
| R_413       | -5.12      | no           | no              |
| A_415       | -6.46      | no           | no              |
| F_416       | -8.85      | yes          | no              |
| E_417       | -12.72     | yes          | no              |
| F_418       | -9.3       | yes          | no              |
| P_419       | -12.07     | yes          | no              |
| Q_420       | -10.25     | yes          | no              |
| I_421       | -11.1      | yes          | no              |
| N_422       | -11.64     | yes          | no              |
| Y_423       | -12.36     | yes          | no              |
| K_425       | -7.94      | in_gray_zone | no              |
| Y_426       | -8.81      | yes          | no              |
| C_427       | -3.98      | no           | no              |
| G_428       | -9.26      | yes          | no              |
| K_429       | -8.72      | yes          | no              |
| P_430       | -6.97      | in_gray_zone | no              |
| Y_431       | -12.71     | yes          | no              |
| T_432       | -5.97      | no           | no              |
| Y_433       | -12.61     | yes          | no              |
| Y_435       | -15.2      | yes          | no              |
| G_436       | -10.54     | yes          | no              |
| L_437       | -7.93      | in_gray_zone | no              |
| G_438       | -9.6       | yes          | no              |
| L_439       | -7.72      | in_gray_zone | no              |
| N_440       | -7.66      | in_gray_zone | no              |
| H_441       | -9.72      | yes          | no              |
| F_442       | -7.97      | in_gray_zone | no              |
| V_443       | -8.78      | yes          | no              |
| D_445       | -9.94      | yes          | no              |
| R_446       | -5.45      | no           | no              |
| L_447       | -10.62     | yes          | no              |
| C_448       | -5.02      | no           | no              |
| K_449       | -10.88     | yes          | no              |
| L_450       | -6.58      | in_gray_zone | no              |
| N_451       | -9.06      | yes          | no              |
| V_452       | -8.37      | in_gray_zone | no              |
| K_453       | -6.31      | no           | no              |
| K_455       | -9.42      | yes          | no              |
| E_456       | -8.01      | in_gray_zone | no              |
| T_457       | -5.49      | no           | no              |
| W_458       | -1.18      | no           | no              |
| V_459       | -5.87      | no           | no              |
| W_460       | -11.75     | yes          | no              |
| Q_461       | -6.31      | no           | no              |
| E_462       | -8.53      | yes          | no              |
| P_463       | -6.45      | no           | no              |
| S_465       | -5.83      | no           | no              |
| Y_466       | -10.01     | yes          | no              |
| P_467       | -11.51     | yes          | no              |
| S_468       | -10.86     | yes          | no              |
| E_469       | -14.71     | yes          | no              |
| P_470       | -12.03     | yes          | no              |
| I_471       | -10.65     | yes          | no              |
| F_472       | -14.08     | yes          | no              |
| V_473       | -12.53     | yes          | no              |
| H_475       | -3.57      | no           | no              |
| P_476       | -10.32     | yes          | no              |
| D_477       | -7.31      | in_gray_zone | no              |
| A_478       | -8.92      | yes          | no              |
| L_479       | -4.19      | no           | no              |
| E_480       | -8.54      | yes          | no              |
| E_481       | -13.29     | yes          | no              |
| D_482       | -14.04     | yes          | no              |
| D_483       | -10.24     | yes          | no              |
| V_485       | -12.32     | yes          | no              |
| V_486       | -10.92     | yes          | no              |
| L_487       | -12.17     | yes          | no              |
| S_488       | -11.6      | yes          | no              |
| V_489       | -8.37      | in_gray_zone | no              |
| V_490       | -11.53     | yes          | no              |
| V_491       | -9.98      | yes          | no              |
| S_492       | -6.18      | no           | no              |
| P_493       | -8.48      | in_gray_zone | no              |
| A_495       | -3.98      | no           | no              |
| G_496       | -5.09      | no           | no              |
| Q_497       | -3.93      | no           | no              |
| K_498       | -6.99      | in_gray_zone | no              |
| P_499       | -7.5       | in_gray_zone | no              |
| A_500       | -7.71      | in_gray_zone | no              |
| Y_501       | -7.3       | in_gray_zone | no              |
| L_502       | -12.57     | yes          | no              |
| L_503       | -11.8      | yes          | no              |
| I_504       | -11.37     | yes          | no              |
| N_506       | -9.07      | yes          | no              |
| A_507       | -11.77     | yes          | no              |
| K_508       | -8.21      | in_gray_zone | no              |
| D_509       | -9.33      | yes          | no              |
| L_510       | -9.54      | yes          | no              |
| S_511       | -7.14      | in_gray_zone | no              |
| E_512       | -11.68     | yes          | no              |
| V_513       | -6.92      | in_gray_zone | no              |
| A_514       | -9.36      | yes          | no              |
| A_516       | -10.53     | yes          | no              |
| E_517       | -8.03      | in_gray_zone | no              |
| V_518       | -9.05      | yes          | no              |
| E_519       | -6.06      | no           | no              |
| I_520       | -4.65      | no           | no              |
| N_521       | -5.45      | no           | no              |
| I_522       | -9.93      | yes          | no              |
| P_523       | -10.69     | yes          | no              |
| V_524       | -7.82      | in_gray_zone | no              |
| F_526       | -10.62     | yes          | no              |
| H_527       | -13.88     | yes          | no              |
| G_528       | -13.02     | yes          | no              |
| L_529       | -9.68      | yes          | no              |
| F_530       | -12.78     | yes          | no              |
| K_531       | -7.74      | in_gray_zone | no              |
| K_532       | -7.55      | in_gray_zone | no              |
| S_533       | -4.92      | no           | no              |
