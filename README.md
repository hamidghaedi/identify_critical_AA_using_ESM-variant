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


The resulting data is available below and also for download as the file named "RPE65_critical_residues.csv" from this repo.

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
| K_13        | -6.18      | no           | no              |
| K_14        | -6.74      | in_gray_zone | no              |
| L_15        | -12.45     | yes          | no              |
| F_16        | -10.67     | yes          | no              |
| E_17        | -7.67      | in_gray_zone | no              |
| T_18        | -8.82      | yes          | no              |
| V_19        | -8.16      | in_gray_zone | no              |
| E_20        | -8.22      | in_gray_zone | no              |
| E_21        | -12.52     | yes          | no              |
| L_22        | -7.1       | in_gray_zone | no              |
| S_23        | -5.29      | no           | no              |
| S_24        | -5.16      | no           | no              |
| P_25        | -11.71     | yes          | no              |
| L_26        | -9.27      | yes          | no              |
| T_27        | -5.2       | no           | no              |
| A_28        | -11.78     | yes          | no              |
| H_29        | -1.83      | no           | no              |
| V_30        | -10.71     | yes          | no              |
| T_31        | -6.09      | no           | no              |
| G_32        | -13.1      | yes          | no              |
| R_33        | -5.56      | no           | no              |
| I_34        | -11.98     | yes          | no              |
| P_35        | -14.36     | yes          | no              |
| L_36        | -3.12      | no           | no              |
| W_37        | -13.88     | yes          | no              |
| L_38        | -13.61     | yes          | no              |
| T_39        | -5.54      | no           | no              |
| G_40        | -15.28     | yes          | no              |
| S_41        | -8.06      | in_gray_zone | no              |
| L_42        | -11.6      | yes          | no              |
| L_43        | -12.33     | yes          | no              |
| R_44        | -14.24     | yes          | no              |
| C_45        | -8.39      | in_gray_zone | no              |
| G_46        | -14.67     | yes          | no              |
| P_47        | -13.32     | yes          | no              |
| G_48        | -13.43     | yes          | no              |
| L_49        | -9.99      | yes          | no              |
| F_50        | -11.16     | yes          | no              |
| E_51        | -10.99     | yes          | no              |
| V_52        | -11.23     | yes          | no              |
| G_53        | -12.32     | yes          | no              |
| S_54        | -6.13      | no           | no              |
| E_55        | -10.79     | yes          | no              |
| P_56        | -7.74      | in_gray_zone | no              |
| F_57        | -11.36     | yes          | no              |
| Y_58        | -6.02      | no           | no              |
| H_59        | -13.2      | yes          | no              |
| L_60        | -10.57     | yes          | no              |
| F_61        | -14.2      | yes          | no              |
| D_62        | -14.95     | yes          | no              |
| G_63        | -13.15     | yes          | no              |
| Q_64        | -8.97      | yes          | no              |
| A_65        | -13.23     | yes          | no              |
| L_66        | -12.82     | yes          | no              |
| L_67        | -13.86     | yes          | no              |
| H_68        | -10.39     | yes          | no              |
| K_69        | -11.39     | yes          | no              |
| F_70        | -15.53     | yes          | no              |
| D_71        | -7.32      | in_gray_zone | no              |
| F_72        | -9.93      | yes          | no              |
| K_73        | -8.65      | yes          | no              |
| E_74        | -8.77      | yes          | no              |
| G_75        | -12.65     | yes          | no              |
| H_76        | -8.74      | yes          | no              |
| V_77        | -13.52     | yes          | no              |
| T_78        | -10.09     | yes          | no              |
| Y_79        | -14.17     | yes          | no              |
| H_80        | -5.9       | no           | no              |
| R_81        | -7.65      | in_gray_zone | no              |
| R_82        | -11.32     | yes          | no              |
| F_83        | -12.5      | yes          | no              |
| I_84        | -8.4       | in_gray_zone | no              |
| R_85        | -7.23      | in_gray_zone | no              |
| T_86        | -9.46      | yes          | no              |
| D_87        | -10.36     | yes          | no              |
| A_88        | -9.39      | yes          | no              |
| Y_89        | -14.08     | yes          | no              |
| V_90        | -6.17      | no           | no              |
| R_91        | -8.78      | yes          | no              |
| A_92        | -7.92      | in_gray_zone | no              |
| M_93        | -8.54      | yes          | no              |
| T_94        | -5.25      | no           | no              |
| E_95        | -7.68      | in_gray_zone | no              |
| K_96        | -8.42      | in_gray_zone | no              |
| R_97        | -11.73     | yes          | no              |
| I_98        | -11.39     | yes          | no              |
| V_99        | -13.53     | yes          | no              |
| I_100       | -9.06      | yes          | no              |
| T_101       | -8.89      | yes          | no              |
| E_102       | -12.08     | yes          | no              |
| F_103       | -11.64     | yes          | no              |
| G_104       | -13.39     | yes          | no              |
| T_105       | -12.93     | yes          | no              |
| C_106       | -3.26      | no           | no              |
| A_107       | -11.25     | yes          | yes             |
| F_108       | -8.96      | yes          | no              |
| P_109       | -12.77     | yes          | no              |
| D_110       | -14.65     | yes          | no              |
| P_111       | -14.27     | yes          | no              |
| C_112       | -15.2      | yes          | no              |
| K_113       | -10.76     | yes          | no              |
| N_114       | -12.51     | yes          | no              |
| I_115       | -13.89     | yes          | no              |
| F_116       | -15.52     | yes          | no              |
| S_117       | -9.78      | yes          | no              |
| R_118       | -11.13     | yes          | no              |
| F_119       | -14.51     | yes          | no              |
| F_120       | -10.78     | yes          | no              |
| S_121       | -11.93     | yes          | no              |
| Y_122       | -9.89      | yes          | no              |
| F_123       | -12.11     | yes          | no              |
| R_124       | -8.17      | in_gray_zone | no              |
| G_125       | -10.51     | yes          | yes             |
| V_126       | -8.64      | yes          | no              |
| E_127       | -9.61      | yes          | no              |
| V_128       | -6.13      | no           | no              |
| T_129       | -12.29     | yes          | no              |
| D_130       | -15.79     | yes          | no              |
| N_131       | -13.45     | yes          | no              |
| A_132       | -10.02     | yes          | no              |
| L_133       | -9.63      | yes          | no              |
| V_134       | -12.41     | yes          | no              |
| N_135       | -11.68     | yes          | no              |
| V_136       | -10.93     | yes          | no              |
| Y_137       | -10.51     | yes          | no              |
| P_138       | -8.09      | in_gray_zone | no              |
| V_139       | -10.02     | yes          | no              |
| G_140       | -8.81      | yes          | no              |
| E_141       | -8.21      | in_gray_zone | no              |
| D_142       | -10.56     | yes          | no              |
| Y_143       | -11.29     | yes          | no              |
| Y_144       | -14.56     | yes          | no              |
| A_145       | -11.05     | yes          | no              |
| C_146       | -7.9       | in_gray_zone | no              |
| T_147       | -11.25     | yes          | no              |
| E_148       | -15.08     | yes          | no              |
| T_149       | -11.7      | yes          | no              |
| N_150       | -10.41     | yes          | no              |
| F_151       | -9.46      | yes          | no              |
| I_152       | -13.5      | yes          | no              |
| T_153       | -9.27      | yes          | no              |
| K_154       | -11.24     | yes          | no              |
| I_155       | -12.96     | yes          | no              |
| N_156       | -10.53     | yes          | no              |
| P_157       | -9.48      | yes          | no              |
| E_158       | -7.8       | in_gray_zone | no              |
| T_159       | -10.51     | yes          | no              |
| L_160       | -13.41     | yes          | no              |
| E_161       | -10.75     | yes          | no              |
| T_162       | -11.13     | yes          | no              |
| I_163       | -7.27      | in_gray_zone | no              |
| K_164       | -10.41     | yes          | no              |
| Q_165       | -8.34      | in_gray_zone | no              |
| V_166       | -10.85     | yes          | no              |
| D_167       | -9.15      | yes          | no              |
| L_168       | -7.86      | in_gray_zone | no              |
| C_169       | -6.04      | no           | no              |
| N_170       | -7.87      | in_gray_zone | no              |
| Y_171       | -9.19      | yes          | no              |
| V_172       | -11.9      | yes          | no              |
| S_173       | -8.62      | yes          | no              |
| V_174       | -12.73     | yes          | no              |
| N_175       | -11.33     | yes          | no              |
| G_176       | -8.47      | in_gray_zone | no              |
| A_177       | -9.05      | yes          | no              |
| T_178       | -12.11     | yes          | no              |
| A_179       | -12.16     | yes          | no              |
| H_180       | -15.53     | yes          | yes             |
| P_181       | -11.97     | yes          | no              |
| H_182       | -12.99     | yes          | yes             |
| I_183       | -6.92      | in_gray_zone | no              |
| E_184       | -7.76      | in_gray_zone | no              |
| N_185       | -4.21      | no           | no              |
| D_186       | -13.19     | yes          | no              |
| G_187       | -12.02     | yes          | no              |
| T_188       | -9.94      | yes          | no              |
| V_189       | -8.83      | yes          | no              |
| Y_190       | -11.32     | yes          | no              |
| N_191       | -12.31     | yes          | no              |
| I_192       | -9.96      | yes          | no              |
| G_193       | -12.59     | yes          | no              |
| N_194       | -10.14     | yes          | no              |
| C_195       | -9.04      | yes          | no              |
| F_196       | -11.26     | yes          | no              |
| G_197       | -9.6       | yes          | no              |
| K_198       | -7.87      | in_gray_zone | no              |
| N_199       | -5.62      | no           | no              |
| F_200       | -5.13      | no           | no              |
| S_201       | -4.33      | no           | no              |
| I_202       | -5.6       | no           | no              |
| A_203       | -3.33      | no           | no              |
| Y_204       | -13.28     | yes          | no              |
| N_205       | -8.41      | in_gray_zone | no              |
| I_206       | -11.75     | yes          | no              |
| V_207       | -11.11     | yes          | no              |
| K_208       | -8.85      | yes          | no              |
| I_209       | -9.86      | yes          | no              |
| P_210       | -10.15     | yes          | no              |
| P_211       | -10.08     | yes          | no              |
| L_212       | -4.28      | no           | no              |
| Q_213       | -4.25      | no           | no              |
| A_214       | -4.73      | no           | no              |
| D_215       | -6.59      | in_gray_zone | no              |
| K_216       | -4.44      | no           | no              |
| E_217       | -5.89      | no           | no              |
| D_218       | -7.2       | in_gray_zone | no              |
| P_219       | -7.12      | in_gray_zone | no              |
| I_220       | -5.09      | no           | no              |
| S_221       | -5.01      | no           | no              |
| K_222       | -6.24      | no           | no              |
| S_223       | -5.12      | no           | no              |
| E_224       | -8.64      | yes          | no              |
| I_225       | -8.52      | yes          | no              |
| V_226       | -8.5       | yes          | no              |
| V_227       | -5.91      | no           | no              |
| Q_228       | -4.73      | no           | no              |
| F_229       | -7.64      | in_gray_zone | no              |
| P_230       | -11.28     | yes          | no              |
| C_231       | -5.99      | no           | no              |
| S_232       | -6.99      | in_gray_zone | no              |
| D_233       | -6.25      | no           | no              |
| R_234       | -6.65      | in_gray_zone | no              |
| F_235       | -6.72      | in_gray_zone | no              |
| K_236       | -6.56      | in_gray_zone | no              |
| P_237       | -10.27     | yes          | no              |
| S_238       | -10.18     | yes          | no              |
| Y_239       | -13.94     | yes          | no              |
| V_240       | -7.03      | in_gray_zone | no              |
| H_241       | -13.99     | yes          | no              |
| S_242       | -12.7      | yes          | no              |
| F_243       | -14.22     | yes          | no              |
| G_244       | -10.91     | yes          | no              |
| L_245       | -9.68      | yes          | no              |
| T_246       | -14.16     | yes          | no              |
| P_247       | -7.22      | in_gray_zone | no              |
| N_248       | -13.21     | yes          | no              |
| Y_249       | -12.94     | yes          | no              |
| I_250       | -10.47     | yes          | no              |
| V_251       | -12.9      | yes          | no              |
| F_252       | -12.56     | yes          | no              |
| V_253       | -10.91     | yes          | no              |
| E_254       | -12.88     | yes          | no              |
| T_255       | -4.64      | no           | no              |
| P_256       | -12.64     | yes          | no              |
| V_257       | -7.06      | in_gray_zone | no              |
| K_258       | -10.99     | yes          | no              |
| I_259       | -9.58      | yes          | no              |
| N_260       | -10.5      | yes          | no              |
| L_261       | -9.94      | yes          | no              |
| F_262       | -6.88      | in_gray_zone | no              |
| K_263       | -10.57     | yes          | no              |
| F_264       | -7.54      | in_gray_zone | no              |
| L_265       | -9.03      | yes          | no              |
| S_266       | -5.89      | no           | no              |
| S_267       | -7.7       | in_gray_zone | no              |
| W_268       | -10.22     | yes          | no              |
| S_269       | -3.68      | no           | no              |
| L_270       | -9.01      | yes          | no              |
| W_271       | -5.92      | no           | no              |
| G_272       | -9.73      | yes          | no              |
| A_273       | -7.69      | in_gray_zone | no              |
| N_274       | -7.93      | in_gray_zone | no              |
| Y_275       | -9.88      | yes          | no              |
| M_276       | -7.94      | in_gray_zone | no              |
| D_277       | -11.33     | yes          | no              |
| C_278       | -11.59     | yes          | no              |
| F_279       | -11.12     | yes          | no              |
| E_280       | -9.26      | yes          | no              |
| S_281       | -4.66      | no           | no              |
| N_282       | -4.84      | no           | no              |
| E_283       | -9.28      | yes          | no              |
| T_284       | -5.51      | no           | no              |
| M_285       | -4.11      | no           | no              |
| G_286       | -10.04     | yes          | no              |
| V_287       | -7.42      | in_gray_zone | no              |
| W_288       | -7.9       | in_gray_zone | no              |
| L_289       | -9.71      | yes          | no              |
| H_290       | -8.3       | in_gray_zone | no              |
| I_291       | -9.63      | yes          | no              |
| A_292       | -7.57      | in_gray_zone | no              |
| D_293       | -11.22     | yes          | no              |
| K_294       | -11.33     | yes          | no              |
| K_295       | -7.16      | in_gray_zone | no              |
| R_296       | -4.41      | no           | no              |
| K_297       | -6.9       | in_gray_zone | no              |
| K_298       | -5.08      | no           | no              |
| Y_299       | -2.11      | no           | no              |
| L_300       | -7.93      | in_gray_zone | no              |
| N_301       | -3.34      | no           | no              |
| N_302       | -1.16      | no           | no              |
| K_303       | -7.35      | in_gray_zone | no              |
| Y_304       | -9.86      | yes          | no              |
| R_305       | -4.14      | no           | no              |
| T_306       | -7.14      | in_gray_zone | no              |
| S_307       | -5.11      | no           | no              |
| P_308       | -9.03      | yes          | no              |
| F_309       | -10.76     | yes          | no              |
| N_310       | -8.56      | yes          | no              |
| L_311       | -7.63      | in_gray_zone | no              |
| F_312       | -12.21     | yes          | no              |
| H_313       | -13.86     | yes          | no              |
| H_314       | -10.41     | yes          | no              |
| I_315       | -13.7      | yes          | no              |
| N_316       | -13.39     | yes          | no              |
| T_317       | -11.06     | yes          | no              |
| Y_318       | -12.97     | yes          | no              |
| E_319       | -14.56     | yes          | no              |
| D_320       | -10.5      | yes          | no              |
| N_321       | -8.48      | in_gray_zone | no              |
| G_322       | -11.04     | yes          | no              |
| F_323       | -9.6       | yes          | no              |
| L_324       | -11.05     | yes          | no              |
| I_325       | -10.91     | yes          | no              |
| V_326       | -10.78     | yes          | no              |
| D_327       | -15.12     | yes          | no              |
| L_328       | -9.22      | yes          | no              |
| C_329       | -13.19     | yes          | no              |
| C_330       | -10.04     | yes          | no              |
| W_331       | -7.58      | in_gray_zone | no              |
| K_332       | -8.91      | yes          | no              |
| G_333       | -7.73      | in_gray_zone | no              |
| F_334       | -2.59      | no           | no              |
| E_335       | -8.11      | in_gray_zone | no              |
| F_336       | -6.61      | in_gray_zone | no              |
| V_337       | -7.68      | in_gray_zone | no              |
| Y_338       | -6.83      | in_gray_zone | no              |
| N_339       | -7.85      | in_gray_zone | no              |
| Y_340       | -7.52      | in_gray_zone | no              |
| L_341       | -9.66      | yes          | no              |
| Y_342       | -10.89     | yes          | no              |
| L_343       | -8.75      | yes          | no              |
| A_344       | -3.7       | no           | no              |
| N_345       | -10.67     | yes          | no              |
| L_346       | -12.19     | yes          | no              |
| R_347       | -8.3       | in_gray_zone | no              |
| E_348       | -7.28      | in_gray_zone | no              |
| N_349       | -8.8       | yes          | no              |
| W_350       | -7.24      | in_gray_zone | no              |
| E_351       | -8.31      | in_gray_zone | no              |
| E_352       | -6.99      | in_gray_zone | no              |
| V_353       | -7         | in_gray_zone | no              |
| K_354       | -7.05      | in_gray_zone | no              |
| K_355       | -6.31      | no           | no              |
| N_356       | -7.52      | in_gray_zone | no              |
| A_357       | -4.27      | no           | no              |
| R_358       | -5.08      | no           | no              |
| K_359       | -3.97      | no           | no              |
| A_360       | -3.57      | no           | no              |
| P_361       | -7.39      | in_gray_zone | no              |
| Q_362       | -4.43      | no           | no              |
| P_363       | -9.68      | yes          | no              |
| E_364       | -6.87      | in_gray_zone | no              |
| V_365       | -7.03      | in_gray_zone | no              |
| R_366       | -8.84      | yes          | no              |
| R_367       | -14.84     | yes          | no              |
| Y_368       | -12.17     | yes          | no              |
| V_369       | -12.9      | yes          | no              |
| L_370       | -12.48     | yes          | no              |
| P_371       | -13.14     | yes          | no              |
| L_372       | -11.87     | yes          | no              |
| N_373       | -6.77      | in_gray_zone | no              |
| I_374       | -9.05      | yes          | no              |
| D_375       | -6.75      | in_gray_zone | no              |
| K_376       | -6.67      | in_gray_zone | no              |
| A_377       | -5.58      | no           | no              |
| D_378       | -5.1       | no           | no              |
| T_379       | -5.26      | no           | no              |
| G_380       | -9.63      | yes          | no              |
| K_381       | -5.2       | no           | no              |
| N_382       | -12.36     | yes          | no              |
| L_383       | -12.06     | yes          | no              |
| V_384       | -10.93     | yes          | no              |
| T_385       | -6.35      | no           | no              |
| L_386       | -13.51     | yes          | no              |
| P_387       | -10.25     | yes          | no              |
| N_388       | -7.12      | in_gray_zone | no              |
| T_389       | -7.92      | in_gray_zone | no              |
| T_390       | -5.88      | no           | no              |
| A_391       | -10.35     | yes          | no              |
| T_392       | -8.24      | in_gray_zone | no              |
| A_393       | -10.61     | yes          | no              |
| I_394       | -7.27      | in_gray_zone | no              |
| L_395       | -3.35      | no           | no              |
| C_396       | -4.54      | no           | no              |
| S_397       | -5.04      | no           | no              |
| D_398       | -9.72      | yes          | no              |
| E_399       | -6.16      | no           | no              |
| T_400       | -7.01      | in_gray_zone | no              |
| I_401       | -9.36      | yes          | no              |
| W_402       | -9.59      | yes          | no              |
| L_403       | -8.93      | yes          | no              |
| E_404       | -7.57      | in_gray_zone | no              |
| P_405       | -9.54      | yes          | no              |
| E_406       | -11.62     | yes          | no              |
| V_407       | -5.92      | no           | no              |
| L_408       | -10.23     | yes          | no              |
| F_409       | -7.95      | in_gray_zone | no              |
| S_410       | -5.28      | no           | no              |
| G_411       | -9.02      | yes          | no              |
| P_412       | -8.43      | in_gray_zone | no              |
| R_413       | -5.39      | no           | no              |
| Q_414       | -7.73      | in_gray_zone | no              |
| A_415       | -6.8       | in_gray_zone | no              |
| F_416       | -9.31      | yes          | no              |
| E_417       | -13.39     | yes          | yes             |
| F_418       | -9.79      | yes          | no              |
| P_419       | -12.7      | yes          | no              |
| Q_420       | -10.79     | yes          | no              |
| I_421       | -11.68     | yes          | no              |
| N_422       | -12.25     | yes          | no              |
| Y_423       | -13.01     | yes          | no              |
| Q_424       | -4.7       | no           | no              |
| K_425       | -8.36      | in_gray_zone | no              |
| Y_426       | -9.27      | yes          | no              |
| C_427       | -4.19      | no           | no              |
| G_428       | -9.74      | yes          | no              |
| K_429       | -9.18      | yes          | no              |
| P_430       | -7.34      | in_gray_zone | no              |
| Y_431       | -13.38     | yes          | no              |
| T_432       | -6.29      | no           | no              |
| Y_433       | -13.28     | yes          | no              |
| A_434       | -7.3       | in_gray_zone | no              |
| Y_435       | -16        | yes          | no              |
| G_436       | -11.09     | yes          | no              |
| L_437       | -8.34      | in_gray_zone | no              |
| G_438       | -10.11     | yes          | no              |
| L_439       | -8.13      | in_gray_zone | no              |
| N_440       | -8.06      | in_gray_zone | no              |
| H_441       | -10.23     | yes          | no              |
| F_442       | -8.39      | in_gray_zone | no              |
| V_443       | -9.24      | yes          | no              |
| P_444       | -10.22     | yes          | no              |
| D_445       | -10.46     | yes          | no              |
| R_446       | -5.73      | no           | no              |
| L_447       | -11.18     | yes          | no              |
| C_448       | -5.28      | no           | no              |
| K_449       | -11.45     | yes          | no              |
| L_450       | -6.92      | in_gray_zone | no              |
| N_451       | -9.53      | yes          | no              |
| V_452       | -8.81      | yes          | no              |
| K_453       | -6.64      | in_gray_zone | no              |
| T_454       | -10.23     | yes          | no              |
| K_455       | -9.92      | yes          | no              |
| E_456       | -8.43      | in_gray_zone | no              |
| T_457       | -5.78      | no           | no              |
| W_458       | -1.24      | no           | no              |
| V_459       | -6.18      | no           | no              |
| W_460       | -12.37     | yes          | no              |
| Q_461       | -6.64      | in_gray_zone | no              |
| E_462       | -8.98      | yes          | no              |
| P_463       | -6.79      | in_gray_zone | no              |
| D_464       | -6.18      | no           | no              |
| S_465       | -6.14      | no           | no              |
| Y_466       | -10.54     | yes          | no              |
| P_467       | -12.11     | yes          | no              |
| S_468       | -11.43     | yes          | no              |
| E_469       | -15.49     | yes          | no              |
| P_470       | -12.66     | yes          | no              |
| I_471       | -11.21     | yes          | no              |
| F_472       | -14.83     | yes          | no              |
| V_473       | -13.19     | yes          | no              |
| S_474       | -6.87      | in_gray_zone | no              |
| H_475       | -3.76      | no           | no              |
| P_476       | -10.87     | yes          | no              |
| D_477       | -7.69      | in_gray_zone | no              |
| A_478       | -9.39      | yes          | no              |
| L_479       | -4.41      | no           | no              |
| E_480       | -8.98      | yes          | no              |
| E_481       | -13.99     | yes          | no              |
| D_482       | -14.77     | yes          | no              |
| D_483       | -10.77     | yes          | no              |
| G_484       | -14.45     | yes          | no              |
| V_485       | -12.97     | yes          | no              |
| V_486       | -11.5      | yes          | no              |
| L_487       | -12.81     | yes          | no              |
| S_488       | -12.21     | yes          | no              |
| V_489       | -8.81      | yes          | no              |
| V_490       | -12.14     | yes          | no              |
| V_491       | -10.51     | yes          | no              |
| S_492       | -6.51      | in_gray_zone | no              |
| P_493       | -8.93      | yes          | no              |
| G_494       | -7.02      | in_gray_zone | no              |
| A_495       | -4.19      | no           | no              |
| G_496       | -5.36      | no           | no              |
| Q_497       | -4.13      | no           | no              |
| K_498       | -7.36      | in_gray_zone | no              |
| P_499       | -7.89      | in_gray_zone | no              |
| A_500       | -8.11      | in_gray_zone | no              |
| Y_501       | -7.69      | in_gray_zone | no              |
| L_502       | -13.23     | yes          | no              |
| L_503       | -12.42     | yes          | no              |
| I_504       | -11.97     | yes          | no              |
| L_505       | -13.27     | yes          | no              |
| N_506       | -9.55      | yes          | no              |
| A_507       | -12.39     | yes          | no              |
| K_508       | -8.64      | yes          | no              |
| D_509       | -9.82      | yes          | no              |
| L_510       | -10.05     | yes          | no              |
| S_511       | -7.52      | in_gray_zone | no              |
| E_512       | -12.29     | yes          | no              |
| V_513       | -7.28      | in_gray_zone | no              |
| A_514       | -9.86      | yes          | no              |
| R_515       | -11.25     | yes          | no              |
| A_516       | -11.08     | yes          | no              |
| E_517       | -8.45      | in_gray_zone | no              |
| V_518       | -9.52      | yes          | no              |
| E_519       | -6.38      | no           | no              |
| I_520       | -4.89      | no           | no              |
| N_521       | -5.74      | no           | no              |
| I_522       | -10.45     | yes          | no              |
| P_523       | -11.25     | yes          | no              |
| V_524       | -8.23      | in_gray_zone | no              |
| T_525       | -9.85      | yes          | no              |
| F_526       | -11.18     | yes          | no              |
| H_527       | -14.61     | yes          | yes             |
| G_528       | -13.7      | yes          | no              |
| L_529       | -10.18     | yes          | no              |
| F_530       | -13.45     | yes          | no              |
| K_531       | -8.15      | in_gray_zone | no              |
| K_532       | -7.94      | in_gray_zone | no              |
| S_533       | -5.18      | no           | no              |



Steps in this guide

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
