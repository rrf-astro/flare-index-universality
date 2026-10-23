
## 1. Run header

Run started (UTC): 2026-10-23 21:37:26

Mode: full

Seed: 20261012

Package versions: python 3.13.16, numpy 2.5.3, pandas 3.0.5, scipy 1.18.1, astropy 8.0.2, astroquery 0.4.11, powerlaw 2.0.0, matplotlib 3.11.2

Data sources:

| catalog | rows | accessed_utc | source | sha256 |
|---|---|---|---|---|
| J/ApJS/241/29/table2 | 162262 | 2026-10-23 21:37:31 | local cache of an earlier VizieR download | 78a080286f8ed6b1ef3122342bd4a80ad599d93215ac64c4a5f181007d50e9b4 |
| J/ApJS/241/29/table1 | 3420 | 2026-10-23 21:37:32 | local cache of an earlier VizieR download | bf0dd52df3f5900f87be0dbabfdfe19fac01d7ac7ffd911ed5f93e26b323224d |
| J/ApJS/229/30/catalog | 197096 | 2026-10-23 21:37:32 | local cache of an earlier VizieR download | 607e868585f0e92e3ad56e34f98e6d836413ea5368b29d9151fe39326d2149e8 |


## 2. Data integrity

| check | status | detail |
|---|---|---|
| row count flares | ok | 162262 rows, expected 162262 |
| row count stars | ok | 3420 rows, expected 3420 |
| row count params | ok | 197096 rows, expected 197096 |
| logE finite | ok | 0 non finite values |
| logE range | ok | min 29.3066, max 38.5562, catalog range 29.3 to 38.6 |
| every flare KIC is in the star table | ok | 0 flare KIC values missing from table1 |
| Nf equals the number of flares per star | ok | 0 stars with inconsistent counts |
| duplicate stars in table1 | ok | 0 duplicated KIC values |
| duplicate stars in the Mathur table | ok | 0 duplicated KIC values |
| duplicate flares | ok | 8 exact duplicate flares (same KIC, begin, end and logE), reported and kept |
| flare to bolometric luminosity ratio reflects the 2020 erratum | ok | KIC 8935655 has Lum 1.90e-04. The corrected value is about 1.9e-4 and the value before the erratum is 5.2e-7. Lum is not used in the analysis |
| stars matched in the Mathur table | ok | 3302 of 3420 matched, 118 unmatched and excluded (unclassified targets are absent from the Mathur table) |
| Teff of table1 against Teff of Mathur | ok | median difference 0.0 K, standard deviation 13.5 K, 3 stars differ by more than 50 K |
| log g of table1 against log g of Mathur | ok | median difference 0.000, standard deviation 0.0036 |
| provenance of the Mathur parameters among flaring stars | ok | counts DSEP 3248, BTSL 54. BTSL means polynomial fits to BT-Settl models for cool dwarfs, and those stars have no distance or extinction |


## 3. Selection funnel (baseline configuration)

Configuration: log g at least 4.0, radius at most 1.5 times the median of its 200 K bin, class limits (3850.0, 5300.0, 5900.0, 7200.0) K.

| step | stars | flares |
|---|---|---|
| 1 all flares in table2 | 3420 | 162262 |
| 2 star has a counterpart in the Mathur table | 3302 | 161347 |
| 3 star passes the dwarf filter | 2894 | 138414 |
| 3b star not removed by the exclusion options | 2894 | 138414 |
| 4 Teff inside a class | 2872 | 138150 |
| 5 finite energy, analysed sample | 2872 | 138150 |

Classes of the analysed sample:

| class | stars | flares | logE min | logE median | logE max | Teff min (K) | Teff max (K) |
|---|---|---|---|---|---|---|---|
| M | 549 | 44151 | 29.31 | 32.54 | 35.18 | 2500 | 3849 |
| K | 1232 | 62553 | 31.38 | 33.62 | 38.18 | 3850 | 5297 |
| G | 616 | 22508 | 31.99 | 34.13 | 36.42 | 5300 | 5897 |
| F | 475 | 8938 | 33.05 | 34.34 | 36.81 | 5900 | 7199 |


## 4. Sample size under the filter and class variants

Stars and flares per class. These counts are upper bounds for the tail used in the fits.

| variant | M stars | M flares | K stars | K flares | G stars | G flares | F stars | F flares |
|---|---|---|---|---|---|---|---|---|
| baseline | 549 | 44151 | 1232 | 62553 | 616 | 22508 | 475 | 8938 |
| no filter | 551 | 44209 | 1398 | 74877 | 738 | 29878 | 572 | 11718 |
| log g 4.0 only | 551 | 44209 | 1236 | 62673 | 642 | 23248 | 491 | 9174 |
| log g 4.2 only | 551 | 44209 | 1232 | 62553 | 612 | 22308 | 375 | 6848 |
| radius factor 1.25 | 525 | 41677 | 1229 | 62539 | 585 | 21800 | 422 | 7717 |
| radius factor 2.0 | 551 | 44209 | 1233 | 62558 | 642 | 23248 | 491 | 9174 |
| strict 1 sigma interval | 545 | 43870 | 1224 | 62349 | 566 | 20929 | 376 | 6827 |
| Yang class limits | 669 | 50339 | 808 | 42489 | 1017 | 38877 | 382 | 6483 |
| limits minus 200 K | 360 | 29531 | 1191 | 66793 | 666 | 26720 | 647 | 14706 |
| limits plus 200 K | 685 | 51341 | 1322 | 64121 | 597 | 18242 | 271 | 4481 |
| Pecaut limits unrounded | 549 | 44151 | 1213 | 61713 | 654 | 23911 | 456 | 8375 |
| drop single flare stars | 548 | 44150 | 1225 | 62546 | 588 | 22480 | 447 | 8910 |
| drop neighbour stars | 532 | 43188 | 1196 | 61523 | 592 | 21745 | 437 | 8031 |
| drop stars without rotation | 507 | 42695 | 1190 | 62199 | 564 | 22271 | 384 | 8534 |
| drop non DSEP stars | 495 | 42022 | 1232 | 62553 | 616 | 22508 | 475 | 8938 |

Subdivision of the M class at 3200 K (Mhi between 3200 and 3850 K, Mfc below 3200 K):

| class | stars | flares | logE min | logE median | logE max | Teff min (K) | Teff max (K) |
|---|---|---|---|---|---|---|---|
| Mhi | 496 | 41926 | 30.34 | 32.58 | 35.18 | 3200 | 3849 |
| Mfc | 53 | 2225 | 29.31 | 31.37 | 33.93 | 2500 | 3193 |


## 5. Validation on synthetic data

Each estimator is tested on samples with known truth. Sample sizes depend on the run mode.

| test | status | detail |
|---|---|---|
| alpha agrees with powerlaw at the same xmin | ok | largest difference 2.52e-13 |
| grid scan agrees with the full scan | ok | largest difference 0.0080, smallest Hill sigma 0.0241 |
| recovery of alpha without bias | ok | largest mean bias 0.0017 over 200 samples per alpha |
| truncated fit recovers alpha and the pure fit is biased | ok | bias M1 0.0202, bias M2 -0.0024 |
| goodness of fit keeps its size under a true power law | ok | 100 samples, rejection rate 0.06 at 0.1 and 0.02 at 0.05, mean p 0.59 |
| goodness of fit power against the lognormal is above its size and grows with the sample size | ok | rejection rate 0.58 at 30000 points and 0.67 at 200000 points, size 0.06. Limited power is a documented limitation of the test at the KS xmin |
| power-law data are not wrongly assigned to the lognormal | ok | 0 of 100 samples with R negative and p below 0.1 (expected near 0.1 or lower) |
| likelihood ratio power against the lognormal grows with the sample size | ok | fraction 0.49 at 30000 points and 0.72 at 200000 points |
| cutoff data are separated from the pure power law | ok | 18 of 20 samples with the nested test p below 0.1 |

Recovery of alpha and xmin (V2):

| true alpha | mean bias | scatter | mean Hill error | coverage of 1 Hill error | mean log10 xmin minus 33 |
|---|---|---|---|---|---|
| 1.70 | -0.0004 | 0.0190 | 0.0144 | 0.60 | 0.133 |
| 2.00 | 0.0017 | 0.0269 | 0.0205 | 0.61 | 0.088 |
| 2.40 | 0.0002 | 0.0304 | 0.0284 | 0.68 | 0.058 |

Power of the goodness of fit test for lognormal data (V4):

| lognormal sample size | samples | mean points in the tail | rejection rate at 0.1 |
|---|---|---|---|
| 30000 | 12 | 1166 | 0.58 |
| 200000 | 12 | 2281 | 0.67 |

Likelihood ratio against the lognormal for lognormal data (V5):

| lognormal sample size | samples | lognormal favoured with p below 0.1 |
|---|---|---|
| 30000 | 100 | 0.49 |
| 200000 | 25 | 0.72 |

Truncated data (V3):

| model | mean alpha | scatter | bias |
|---|---|---|---|
| M1 pure | 1.8202 | 0.0081 | 0.0202 |
| M2 truncated | 1.7976 | 0.0110 | -0.0024 |


## 6. Baseline fits in energy as published

Maximum likelihood fits with xmin from the KS distance. M1 is the pure power law and M2 has an upper truncation at the largest energy of the class. The Hill error ignores the uncertainty on xmin and is only a reference. The goodness of fit p uses 100 synthetic samples, so its precision is limited. Positive R favours the power law over the alternative.

| class | model | points | tail | log10 xmin | alpha | Hill error | KS D | GoF p | R vs lognormal | p (Vuong) | R vs cutoff | p (nested) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| M | M1 | 44151 | 2490 | 33.58 | 2.352 | 0.027 | 0.0285 | 0.00 | -14.9 | 0.001 | -19.8 | 3.16e-10 |
| M | M2 | 44151 | 2490 | 33.58 | 2.297 | 0.026 | 0.0192 | 0.00 | 4.1 | 0.314 | -19.8 | 3.16e-10 |
| K | M1 | 62553 | 1873 | 34.82 | 2.731 | 0.040 | 0.0334 | 0.00 | -0.7 | 0.802 | 0.0 | 1 |
| K | M2 | 62553 | 1873 | 34.82 | 2.731 | 0.040 | 0.0334 | 0.00 | -0.7 | 0.803 | 0.0 | 0.998 |
| G | M1 | 22508 | 2734 | 34.70 | 2.516 | 0.029 | 0.0259 | 0.00 | -13.3 | 0.003 | -13.1 | 3.12e-07 |
| G | M2 | 22508 | 2734 | 34.70 | 2.491 | 0.029 | 0.0219 | 0.00 | -6.2 | 0.157 | -13.1 | 3.12e-07 |
| F | M1 | 8938 | 3506 | 34.46 | 2.072 | 0.018 | 0.0261 | 0.00 | -28.8 | 0.000 | -38.0 | 2.85e-18 |
| F | M2 | 8938 | 3364 | 34.48 | 2.060 | 0.018 | 0.0216 | 0.00 | -15.7 | 0.004 | -35.9 | 2.3e-17 |

Alpha of the pure power law (M1) for fixed values of xmin, with at least 300 flares above xmin. A plateau would indicate a power law over that range:

| class | log10 xmin | alpha | tail | Hill error | KS D |
|---|---|---|---|---|---|
| M | 33.0 | 2.097 | 10198 | 0.011 | 0.0341 |
| M | 33.5 | 2.301 | 3057 | 0.024 | 0.0310 |
| M | 34.0 | 2.551 | 708 | 0.058 | 0.0381 |
| K | 33.0 | 1.560 | 54928 | 0.002 | 0.1494 |
| K | 33.5 | 1.826 | 36628 | 0.004 | 0.0755 |
| K | 34.0 | 2.103 | 16147 | 0.009 | 0.0519 |
| K | 34.5 | 2.419 | 4852 | 0.020 | 0.0438 |
| K | 35.0 | 2.876 | 959 | 0.061 | 0.0496 |
| G | 33.0 | 1.375 | 22409 | 0.003 | 0.3002 |
| G | 33.5 | 1.604 | 20994 | 0.004 | 0.1799 |
| G | 34.0 | 1.985 | 13854 | 0.008 | 0.0841 |
| G | 34.5 | 2.411 | 4992 | 0.020 | 0.0286 |
| G | 35.0 | 2.697 | 1018 | 0.053 | 0.0394 |
| F | 33.0 | 1.323 | 8938 | 0.003 | 0.3459 |
| F | 33.5 | 1.479 | 8822 | 0.005 | 0.2420 |
| F | 34.0 | 1.798 | 7137 | 0.009 | 0.1015 |
| F | 34.5 | 2.081 | 3202 | 0.019 | 0.0278 |
| F | 35.0 | 2.267 | 1003 | 0.040 | 0.0748 |

Grid and full scan of xmin on random subsamples of 8000 flares, and alpha from the powerlaw package at the same xmin:

| class | points | alpha grid | alpha full scan | log10 xmin grid | log10 xmin full | alpha powerlaw | Hill error | full scan seconds |
|---|---|---|---|---|---|---|---|---|
| M | 8000 | 2.197 | 2.191 | 33.14 | 33.13 | 2.197 | 0.033 | 0.2 |
| K | 8000 | 2.345 | 2.345 | 34.32 | 34.32 | 2.345 | 0.042 | 0.3 |
| G | 8000 | 2.387 | 2.390 | 34.50 | 34.50 | 2.387 | 0.033 | 0.2 |
| F | 8000 | 2.059 | 2.058 | 34.44 | 34.44 | 2.059 | 0.019 | 0.2 |

Run time planning from the fit time on 60000 points (M1 0.06 s, M2 0.10 s):

| task | estimated minutes (serial) | estimated minutes (parallel) |
|---|---|---|
| bootstrap, one class, B = 1000 | 1.7 | 0.4 |
| all four classes, B = 1000 | 6.8 | 1.7 |
| equality test calibration, B_NULL = 500, four classes | 3.4 | 0.8 |


## 7. Luminosity-scaled energy

Peak of the differential flare energy distribution in energy E and in the luminosity-scaled energy E_L. The spread between classes is peak log E, 0.1 dex: 1.64 dex, peak log E, 0.3 dex: 1.94 dex, peak log E_L, 0.1 dex: 0.22 dex, peak log E_L, 0.3 dex: 0.32 dex.

| class | median log L | peak log E, 0.1 dex | peak log E, 0.3 dex | peak log E_L, 0.1 dex | peak log E_L, 0.3 dex |
|---|---|---|---|---|---|
| M | -1.76 | 32.66 | 32.46 | 34.11 | 34.11 |
| K | -0.78 | 33.53 | 33.63 | 34.27 | 34.17 |
| G | -0.24 | 34.14 | 33.94 | 34.32 | 34.42 |
| F | 0.17 | 34.30 | 34.40 | 34.13 | 34.13 |


## 8. Primary analysis in the luminosity-scaled energy

Lower limit of E_L: log10 34.5 (all classes), window of the truncated fit: log10 34.5 to 35.5. Bootstrap by star with B = 1000. Errors are bootstrap standard errors, and the intervals are percentile intervals.

Index of the power law with exponential cutoff (primary model):

| class | stars | flares | alpha | se | 16 to 84 low | 16 to 84 high | 99 low | 99 high |
|---|---|---|---|---|---|---|---|---|
| M | 549 | 14533 | 2.007 | 0.033 | 1.965 | 2.029 | 1.905 | 2.070 |
| K | 1232 | 25651 | 1.931 | 0.107 | 1.745 | 1.990 | 1.696 | 2.044 |
| G | 616 | 8267 | 1.947 | 0.057 | 1.885 | 2.000 | 1.793 | 2.087 |
| F | 475 | 1902 | 2.268 | 0.213 | 2.061 | 2.468 | 1.604 | 2.682 |

Index of the truncated power law in the window:

| class | stars | flares | alpha | se | 16 to 84 low | 16 to 84 high | 99 low | 99 high |
|---|---|---|---|---|---|---|---|---|
| M | 549 | 13590 | 1.951 | 0.022 | 1.930 | 1.973 | 1.900 | 2.008 |
| K | 1232 | 23593 | 1.836 | 0.020 | 1.817 | 1.857 | 1.789 | 1.888 |
| G | 616 | 7951 | 2.071 | 0.043 | 2.029 | 2.115 | 1.966 | 2.185 |
| F | 475 | 1859 | 2.383 | 0.141 | 2.267 | 2.553 | 2.033 | 2.751 |

Index of the pure power law above the lower limit:

| class | stars | flares | alpha | se | 16 to 84 low | 16 to 84 high | 99 low | 99 high |
|---|---|---|---|---|---|---|---|---|
| M | 549 | 14533 | 2.102 | 0.021 | 2.081 | 2.124 | 2.051 | 2.156 |
| K | 1232 | 25651 | 2.032 | 0.014 | 2.019 | 2.046 | 2.000 | 2.069 |
| G | 616 | 8267 | 2.250 | 0.035 | 2.216 | 2.286 | 2.168 | 2.349 |
| F | 475 | 1902 | 2.487 | 0.103 | 2.401 | 2.614 | 2.272 | 2.773 |

Tests of equality of the index between classes:

| estimator | classes | statistic | p bootstrap | p chi-squared | common index | largest difference | verdict |
|---|---|---|---|---|---|---|---|
| cutoff index (primary) | M, K, G | 1.11 | 0.6284 | 0.574 | 1.988 | 0.075 | H0 not rejected |
| cutoff index (primary) | M, K, G, F | 2.81 | 0.4715 | 0.421 | 1.993 | 0.337 | H0 not rejected |
| window truncated index | M, K, G | 31.18 | 0.0010 | 1.7e-07 | 1.908 | 0.235 | H0 rejected |
| window truncated index | M, K, G, F | 42.39 | 0.0010 | 3.32e-09 | 1.912 | 0.547 | H0 rejected |
| pure power law index | M, K, G | 35.85 | 0.0010 | 1.65e-08 | 2.072 | 0.218 | H0 rejected |
| pure power law index | M, K, G, F | 51.80 | 0.0010 | 3.3e-11 | 2.077 | 0.454 | H0 rejected |

Pairwise differences of the cutoff index, classes M, K and G:

| pair | difference | se | ci95 low | ci95 high | p | equivalent within margin | p Holm |
|---|---|---|---|---|---|---|---|
| M minus K | 0.075 | 0.112 | -0.051 | 0.309 | 0.5105 | False | 1.0000 |
| M minus G | 0.059 | 0.067 | -0.072 | 0.185 | 0.3746 | False | 1.0000 |
| K minus G | -0.016 | 0.120 | -0.281 | 0.134 | 0.9131 | False | 1.0000 |

Comparison with 3/2:

| class | estimator | alpha | alpha minus 3/2 | 99 low | 99 high | 3/2 rejected |
|---|---|---|---|---|---|---|
| M | cutoff | 2.007 | 0.507 | 1.905 | 2.070 | True |
| M | window | 1.951 | 0.451 | 1.900 | 2.008 | True |
| K | cutoff | 1.931 | 0.431 | 1.696 | 2.044 | True |
| K | window | 1.836 | 0.336 | 1.789 | 1.888 | True |
| G | cutoff | 1.947 | 0.447 | 1.793 | 2.087 | True |
| G | window | 2.071 | 0.571 | 1.966 | 2.185 | True |
| F | cutoff | 2.268 | 0.768 | 1.604 | 2.682 | True |
| F | window | 2.383 | 0.883 | 2.033 | 2.751 | True |

Nested likelihood ratio, separate against common index for the cutoff model, and the design effect of the clustering of flares in stars:

| classes | Lambda | p naive | alpha common |
|---|---|---|---|
| M, K, G | 19.59 | 5.57e-05 | 1.955 |
| M, K, G, F | 38.08 | 2.72e-08 | 1.958 |

| class | Hill error, independent flares | bootstrap error by star, pure power law | design effect, standard error ratio |
|---|---|---|---|
| M | 0.009 | 0.021 | 2.343 |
| K | 0.006 | 0.014 | 2.149 |
| G | 0.014 | 0.035 | 2.557 |
| F | 0.034 | 0.103 | 3.025 |

Goodness of fit of the cutoff model at the fixed lower limit, 200 synthetic samples (flares treated as independent):

| class | KS D | p | n |
|---|---|---|---|
| M | 0.0267 | 0.000 | 14533 |
| K | 0.0306 | 0.000 | 25651 |
| G | 0.0168 | 0.000 | 8267 |
| F | 0.0246 | 0.020 | 1902 |


## 9. Baseline in energy as published, with bootstrap by star

Truncated power law (M2) above the xmin of the KS distance, B = 1000 resamplings of stars with xmin estimated again in each. The goodness of fit comes from the fit of section 6.

| class | stars | tail | log10 xmin | se of log10 xmin | alpha | se | 16 to 84 low | 16 to 84 high | GoF p |
|---|---|---|---|---|---|---|---|---|---|
| M | 549 | 2490 | 33.58 | 0.19 | 2.297 | 0.087 | 2.200 | 2.365 | 0.00 |
| K | 1232 | 1873 | 34.82 | 0.24 | 2.731 | 0.222 | 2.303 | 2.816 | 0.00 |
| G | 616 | 2734 | 34.70 | 0.14 | 2.491 | 0.118 | 2.375 | 2.618 | 0.00 |
| F | 475 | 3364 | 34.48 | 0.10 | 2.060 | 0.165 | 1.916 | 2.238 | 0.00 |

| classes | statistic | p bootstrap | p chi-squared | common index | largest difference | verdict |
|---|---|---|---|---|---|---|
| M, K, G | 4.21 | 0.1808 | 0.122 | 2.398 | 0.433 | H0 not rejected |
| M, K, G, F | 7.81 | 0.0619 | 0.0501 | 2.351 | 0.671 | H0 not rejected |


## 10. Sensitivity analysis

Bootstrap by star with B = 200 resamplings in every variant. The verdict refers to the cutoff index of classes M, K and G.

| variant | alpha M | alpha K | alpha G | alpha F | p cutoff | p window | largest difference | min tail | verdict |
|---|---|---|---|---|---|---|---|---|---|
| baseline | 2.007 | 1.931 | 1.947 | 2.268 | 0.628 | 0.001 | 0.075 | 8267 | H0 not rejected |
| S5 window 34.6 to 35.4 | 2.064 | 2.002 | 2.032 | 2.414 | 0.796 | 0.005 | 0.062 | 6505 | H0 not rejected |
| S5 window 34.8 to 35.6 | 2.168 | 2.130 | 2.175 | 2.566 | 0.950 | 0.005 | 0.045 | 3787 | H0 not rejected |
| S5 window 35.0 to 36.0 | 2.271 | 2.267 | 2.221 | 2.851 | 0.876 | 0.005 | 0.050 | 2026 | H0 not rejected |
| S5 window 34.3 to 35.3 | 1.874 | 1.780 | 1.700 | 2.008 | 0.015 | 0.005 | 0.174 | 12174 | H0 rejected |
| S5 window 34.7 to 35.7 | 2.117 | 2.064 | 2.093 | 2.444 | 0.846 | 0.005 | 0.054 | 4996 | H0 not rejected |
| S2 Yang class limits | 1.991 | 1.751 | 2.152 | 2.109 | 0.005 | 0.005 | 0.401 | 15033 | H0 rejected |
| S2 limits minus 200 K | 2.025 | 1.785 | 2.141 | 2.123 | 0.005 | 0.020 | 0.355 | 9865 | H0 rejected |
| S2 limits plus 200 K | 1.994 | 1.948 | 2.044 | 1.912 | 0.791 | 0.005 | 0.097 | 5783 | H0 not rejected |
| S2 Pecaut limits unrounded | 2.007 | 1.933 | 1.931 | 2.245 | 0.512 | 0.005 | 0.076 | 8756 | H0 not rejected |
| S3 no dwarf filter | 2.007 | 1.907 | 1.918 | 2.206 | 0.328 | 0.005 | 0.101 | 10498 | H0 not rejected |
| S3 log g 4.0 only | 2.007 | 1.932 | 1.934 | 2.270 | 0.582 | 0.005 | 0.075 | 8533 | H0 not rejected |
| S3 log g 4.2 only | 2.007 | 1.931 | 1.946 | 2.188 | 0.617 | 0.005 | 0.076 | 8238 | H0 not rejected |
| S3 radius factor 1.25 | 2.015 | 1.931 | 1.937 | 2.257 | 0.527 | 0.005 | 0.084 | 8127 | H0 not rejected |
| S3 radius factor 2.0 | 2.007 | 1.931 | 1.934 | 2.270 | 0.552 | 0.005 | 0.076 | 8533 | H0 not rejected |
| S3 strict one sigma interval | 2.011 | 1.931 | 1.915 | 2.180 | 0.498 | 0.005 | 0.096 | 7888 | H0 not rejected |
| S4 energy of M times 1.25 | 1.943 | 1.931 | 1.947 | 2.268 | 0.990 | 0.005 | 0.016 | 8267 | H0 not rejected |
| S4 energy of M divided by 1.25 | 2.063 | 1.931 | 1.947 | 2.268 | 0.224 | 0.005 | 0.131 | 8267 | H0 not rejected |
| S4 energy of K times 1.25 | 2.007 | 1.854 | 1.947 | 2.268 | 0.338 | 0.005 | 0.153 | 8267 | H0 not rejected |
| S4 energy of K divided by 1.25 | 2.007 | 1.999 | 1.947 | 2.268 | 0.736 | 0.005 | 0.059 | 8267 | H0 not rejected |
| S4 energy of G times 1.25 | 2.007 | 1.931 | 1.836 | 2.268 | 0.040 | 0.005 | 0.170 | 10125 | H0 rejected |
| S4 energy of G divided by 1.25 | 2.007 | 1.931 | 2.032 | 2.268 | 0.701 | 0.005 | 0.100 | 6560 | H0 not rejected |
| S7 drop single flare stars | 2.007 | 1.931 | 1.946 | 2.264 | 0.547 | 0.005 | 0.076 | 8249 | H0 not rejected |
| S7 drop neighbour stars | 1.973 | 1.932 | 1.938 | 2.446 | 0.861 | 0.005 | 0.040 | 8095 | H0 not rejected |
| S7 drop stars without rotation | 2.003 | 1.930 | 1.944 | 2.162 | 0.642 | 0.005 | 0.073 | 8199 | H0 not rejected |
| S7 drop non DSEP stars | 2.007 | 1.931 | 1.947 | 2.268 | 0.677 | 0.005 | 0.076 | 8267 | H0 not rejected |
| S10 class F below 6500 K | 2.007 | 1.931 | 1.947 | 2.434 | 0.672 | 0.005 | 0.075 | 8267 | H0 not rejected |

Split of class M at 3200 K, upper part of M against K and G, with the index of the fully convective part:

| variant | alpha Mhi | alpha Mfc | alpha K | alpha G | alpha F | p cutoff | p window | largest difference | min tail | verdict |
|---|---|---|---|---|---|---|---|---|---|---|
| S2 M split at 3200 K, Mhi against K and G | 2.01 | 1.74 | 1.931 | 1.947 | 2.268 | 0.731 | 0.005 | 0.079 | 8267 | H0 not rejected |

Lower limits in energy as published (no luminosity scaling):

| variant | alpha M | alpha K | alpha G | alpha F | p cutoff | p window | largest difference | min tail | verdict | alpha pure power law M, K, G, F |
|---|---|---|---|---|---|---|---|---|---|---|
| S1 energy as published, lower limit log10 33.0 | 1.920 | 1.318 | 0.469 | 0.500 | 0.005 | 0.005 | 1.450 | 10198 | H0 rejected | 2.10, 1.56, 1.37, 1.31 |
| S1 energy as published, lower limit log10 33.5 | 2.045 | 1.726 | 0.946 | 0.872 | 0.005 | 0.005 | 1.098 | 3057 | H0 rejected | 2.30, 1.83, 1.60, 1.48 |
| S1 energy as published, lower limit log10 34.0 | 1.947 | 2.087 | 1.590 | 1.497 | 0.015 | 0.005 | 0.497 | 708 | H0 rejected | 2.55, 2.10, 1.98, 1.80 |

Noise in the bootstrap, centred on the bootstrap mean:

| variant | alpha M | mean bootstrap alpha M | alpha K | mean bootstrap alpha K | alpha G | mean bootstrap alpha G | alpha F | mean bootstrap alpha F | p cutoff | p window | largest difference | min tail | verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| S4 noise of 0.2 dex on the energy | 2.007 | 1.939 | 1.931 | 1.821 | 1.947 | 1.843 | 2.268 | 2.117 | 0.662 | 0.005 | 0.075 | 8267 | H0 not rejected |
| S4 noise on the luminosity from the Mathur intervals | 2.007 | 1.981 | 1.931 | 1.870 | 1.947 | 1.864 | 2.268 | 2.085 | 0.746 | 0.005 | 0.075 | 8267 | H0 not rejected |
| S4 both noises | 2.007 | 1.924 | 1.931 | 1.817 | 1.947 | 1.795 | 2.268 | 1.992 | 0.721 | 0.005 | 0.075 | 8267 | H0 not rejected |

Noise of 0.2 dex on the energy in the bootstrap of the baseline in energy as published:

| class | alpha | mean bootstrap alpha with noise | se with noise | se without noise |
|---|---|---|---|---|
| M | 2.297 | 2.293 | 0.094 | 0.087 |
| K | 2.731 | 2.632 | 0.166 | 0.222 |
| G | 2.491 | 2.522 | 0.126 | 0.118 |


## 11. Heterogeneity between stars and reproduction of the published index

Truncated index per star in the E_L window, stars with at least 50 flares in it.

| class | stars | mean index | Q | p Q | between star scatter | median star index |
|---|---|---|---|---|---|---|
| M | 112 | 1.916 | 264.6 | 1.46e-14 | 0.217 | 1.952 |
| K | 152 | 1.757 | 415.4 | 1.54e-26 | 0.251 | 1.806 |
| G | 44 | 1.871 | 273.8 | 2.34e-35 | 0.470 | 2.075 |
| F | 6 | 1.765 | 74.7 | 1.05e-14 | 0.675 | 2.139 |

Bisector of the least squares lines (Isobe et al. 1990) on the differential distribution in energy, with the class limits of Yang and Liu (4000, 5000, 6000 and 7500 K), no dwarf filter and the energy ranges of their Fig. 3. The index of the table is the negative slope. The last two variants use only bins with at least 20 flares and a line weighted by the Poisson counts:

Largest difference to the published index over all variants of each class: M 0.14, K 0.19, G 0.26, F 0.40. The last bin column uses all stars of table 1 of Yang and Liu, including those without Mathur parameters. The flare counts of the classes are compared with Fig. 3 of Yang and Liu in the columns that start with flares.

| class | published index | fit range | bin 0.2 dex | bin 0.25 dex | bin 0.3 dex | bin 0.25 dex, at least 20 flares per bin | bin 0.25 dex, Poisson weighted line | bin 0.25 dex, all stars of table 1 | flares, Mathur matched | flares, all of table 1 | flares, Fig. 3 of Yang | flares in range, 0.25 dex bins |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| M | 2.13 | 33.5 to 35.0 | 2.27 | 2.24 | 2.21 | 2.24 | 2.15 | 2.24 | 50397 | 50439 | 50439 | 4612 |
| K | 1.78 | 33.5 to 35.0 | 1.62 | 1.61 | 1.64 | 1.61 | 1.59 | 1.61 | 47130 | 47088 | 47112 | 24804 |
| G | 1.96 | 34.5 to 36.0 | 1.74 | 1.73 | 1.73 | 1.73 | 1.70 | 1.72 | 54360 | 55275 | 55259 | 18391 |
| F | 2.11 | 34.5 to 36.0 | 1.89 | 1.83 | 1.82 | 1.83 | 1.71 | 1.83 | 8877 | 8877 | 8869 | 4318 |


## 12. Verdict

Primary analysis: H0 not rejected. Bootstrap p value 0.628, common cutoff index 1.988, largest pairwise difference 0.075.

Variants with at least 1000 flares above the lower limit in each class: 26. Variants with another verdict: 4. The verdict depends on the choices of the analysis:

| variant | p cutoff | verdict |
|---|---|---|
| S5 window 34.3 to 35.3 | 0.015 | H0 rejected |
| S2 Yang class limits | 0.005 | H0 rejected |
| S2 limits minus 200 K | 0.005 | H0 rejected |
| S4 energy of G times 1.25 | 0.040 | H0 rejected |

The value 3/2 is outside the 99 per cent interval of the cutoff index in 3 of 3 classes.


## 14. Index in narrow bands of effective temperature

Twelve bands with the same number of stars. Lower limit of E_L log10 34.5, window 34.5 to 35.5, bootstrap by star with B = 1000. A band with fewer than 1000 flares in the tail is flagged and stays out of the tests.

| band | Teff from | Teff to | Teff median | stars | tail | valid | alpha cutoff | se cutoff | alpha window | se window | alpha pure | se pure |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 2500 | 3464 | 3317 | 238 | 6448 | yes | 2.000 | 0.036 | 1.911 | 0.034 | 2.051 | 0.038 |
| 1 | 3464 | 3770 | 3646 | 239 | 6172 | yes | 1.949 | 0.043 | 2.018 | 0.035 | 2.177 | 0.026 |
| 2 | 3770 | 4137 | 3900 | 240 | 5512 | yes | 1.861 | 0.048 | 1.924 | 0.042 | 2.066 | 0.028 |
| 3 | 4137 | 4441 | 4308 | 237 | 4472 | yes | 1.838 | 0.055 | 1.889 | 0.044 | 2.061 | 0.030 |
| 4 | 4441 | 4650 | 4538 | 242 | 5662 | yes | 1.647 | 0.048 | 1.753 | 0.041 | 1.985 | 0.025 |
| 5 | 4650 | 4959 | 4820 | 240 | 4708 | yes | 1.729 | 0.056 | 1.789 | 0.041 | 2.013 | 0.032 |
| 6 | 4959 | 5231 | 5105 | 239 | 4608 | yes | 2.043 | 0.147 | 1.893 | 0.047 | 2.045 | 0.037 |
| 7 | 5231 | 5422 | 5313 | 239 | 4812 | yes | 1.709 | 0.077 | 1.863 | 0.053 | 2.104 | 0.035 |
| 8 | 5422 | 5631 | 5522 | 240 | 3190 | yes | 1.938 | 0.099 | 2.070 | 0.069 | 2.265 | 0.053 |
| 9 | 5631 | 5892 | 5767 | 239 | 2821 | yes | 1.963 | 0.110 | 2.138 | 0.087 | 2.299 | 0.081 |
| 10 | 5892 | 6130 | 6021 | 239 | 1091 | yes | 2.538 | 0.136 | 2.581 | 0.096 | 2.611 | 0.087 |
| 11 | 6130 | 7200 | 6341 | 240 | 857 | no | 1.851 | 0.453 | 2.155 | 0.249 | 2.355 | 0.170 |

Equality of the index in the valid bands, and the slope with Teff:

| estimator | valid bands | equality statistic | p bootstrap | p chi-squared | slope per 1000 K | se | 95 low | 95 high | p slope |
|---|---|---|---|---|---|---|---|---|---|
| cutoff | 11 | 75.66 | 0.0010 | 3.54e-12 | -0.061 | 0.024 | -0.112 | -0.019 | 0.0120 |
| window | 11 | 95.20 | 0.0010 | 4.95e-16 | 0.021 | 0.018 | -0.012 | 0.060 | 0.2408 |
| pure | 11 | 92.81 | 0.0010 | 1.48e-15 | 0.025 | 0.015 | -0.003 | 0.057 | 0.1139 |

Shape of the index as a function of the band, weighted chi-squared with a normal reference p:

| estimator | model | chi-squared | dof | reference p |
|---|---|---|---|---|
| cutoff | constant | 75.7 | 10 | 3.54e-12 |
| cutoff | straight line | 69.0 | 9 | 2.39e-11 |
| cutoff | one value per class | 21.5 | 7 | 0.0031 |
| window | constant | 95.2 | 10 | 4.95e-16 |
| window | straight line | 93.8 | 9 | 2.76e-16 |
| window | one value per class | 27.2 | 7 | 0.000309 |
| pure | constant | 92.8 | 10 | 1.48e-15 |
| pure | straight line | 90.1 | 9 | 1.54e-15 |
| pure | one value per class | 23.2 | 7 | 0.00156 |

Completeness diagnostic: peak of the E_L distribution in each band. Slope of the peak -0.002 (error 0.017) dex per 1000 K, range 34.08 to 34.41. The median rotation period is given as context only.

| band | Teff median | peak log E_L | se of the peak | median Prot (d) |
|---|---|---|---|---|
| 0 | 3317 | 34.16 | 0.00 | 8.7 |
| 1 | 3646 | 34.23 | 0.03 | 10.9 |
| 2 | 3900 | 34.22 | 0.05 | 9.8 |
| 3 | 4308 | 34.24 | 0.07 | 9.6 |
| 4 | 4538 | 34.37 | 0.03 | 7.5 |
| 5 | 4820 | 34.25 | 0.04 | 5.8 |
| 6 | 5105 | 34.25 | 0.06 | 5.2 |
| 7 | 5313 | 34.35 | 0.06 | 4.3 |
| 8 | 5522 | 34.41 | 0.10 | 4.9 |
| 9 | 5767 | 34.17 | 0.04 | 3.7 |
| 10 | 6021 | 34.24 | 0.09 | 3.3 |
| 11 | 6341 | 34.08 | 0.06 | 2.7 |

Variants of the bands, B = 200 resamplings each:

| variant | valid bands | bands | p equality cutoff | slope cutoff | p slope cutoff | p equality window | slope window | p slope window | verdict |
|---|---|---|---|---|---|---|---|---|---|
| baseline, 12 bands | 11 | 12 | 0.001 | -0.061 | 0.012 | 0.001 | 0.021 | 0.241 | H0T rejected |
| 8 quantile bands | 8 | 8 | 0.005 | -0.081 | 0.005 | 0.005 | 0.001 | 0.965 | H0T rejected |
| 16 quantile bands | 14 | 16 | 0.005 | -0.074 | 0.010 | 0.005 | -0.012 | 0.502 | H0T rejected |
| bands of 250 K, merged | 13 | 13 | 0.005 | -0.079 | 0.005 | 0.005 | -0.002 | 0.930 | H0T rejected |
| lower limit 34.6 | 10 | 12 | 0.005 | -0.092 | 0.010 | 0.005 | -0.011 | 0.637 | H0T rejected |
| lower limit 34.8 | 10 | 12 | 0.005 | -0.006 | 0.811 | 0.005 | -0.015 | 0.587 | H0T rejected |
| lower limit at the peak of each band plus 0.3 | 11 | 12 | 0.010 | -0.077 | 0.015 | 0.005 | 0.014 | 0.463 | H0T rejected |
| lower limit at the peak of each band plus 0.5 | 10 | 12 | 0.005 | 0.007 | 0.831 | 0.005 | 0.017 | 0.483 | H0T rejected |
| no dwarf filter | 11 | 12 | 0.005 | -0.069 | 0.005 | 0.005 | 0.005 | 0.756 | H0T rejected |
| without neighbour stars | 11 | 12 | 0.005 | -0.101 | 0.005 | 0.005 | -0.008 | 0.711 | H0T rejected |
| without stars lacking rotation | 11 | 12 | 0.005 | -0.084 | 0.005 | 0.005 | -0.002 | 0.891 | H0T rejected |
| without non DSEP stars | 11 | 12 | 0.005 | -0.055 | 0.040 | 0.005 | 0.013 | 0.458 | H0T rejected |

Measurement noise in Teff, 200 replicates, centred on the bootstrap mean:

| estimator | mean change of the index in a band | p equality | slope per 1000 K | se | p slope |
|---|---|---|---|---|---|
| cutoff | -0.047 | 0.005 | -0.076 | 0.032 | 0.020 |
| window | -0.019 | 0.005 | -0.003 | 0.020 | 0.886 |

Verdict on the equality of the cutoff index across the bands: H0T rejected (p 0.001). Slope -0.061 per 1000 K with error 0.024 (p 0.012). 0 of 11 variants with at least 5 valid bands change the verdict, and the sign of the slope changes in 1. The slope of the cutoff index is significant in 9 of 11 variants, so the linear trend is not robust.


## 15. Additional checks

Common minimum duration. The shortest flare of the catalog has three points above the threshold. Share of the flares by number of points in the analysed sample:

| points | share of flares |
|---|---|
| 3 | 0.410 |
| 4 | 0.242 |
| 5 | 0.129 |
| 6 | 0.075 |
| 7 | 0.046 |
| 8 | 0.028 |

Share of the flares above the lower limit with at least 4 points, by class:

| class | share with at least 4 points |
|---|---|
| M | 0.864 |
| K | 0.851 |
| G | 0.845 |
| F | 0.876 |

Classes with a common minimum duration and with the scale that uses Teff alone, B = 200 resamplings:

| variant | alpha M | alpha K | alpha G | alpha F | p cutoff | p window | largest difference | min tail | verdict |
|---|---|---|---|---|---|---|---|---|---|
| baseline | 2.007 | 1.931 | 1.947 | 2.268 | 0.642 | 0.005 | 0.075 | 8267 | H0 not rejected |
| duration at least 4 points | 1.921 | 1.836 | 1.797 | 2.146 | 0.358 | 0.005 | 0.123 | 6989 | H0 not rejected |
| duration at least 5 points | 1.777 | 1.698 | 1.568 | 1.942 | 0.030 | 0.005 | 0.209 | 5095 | H0 rejected |
| scale by Teff only | 1.937 | 1.967 | 2.000 | 2.447 | 0.622 | 0.005 | 0.063 | 7986 | H0 not rejected |
| scale by Teff only, duration at least 4 points | 1.840 | 1.881 | 1.853 | 2.311 | 0.925 | 0.005 | 0.041 | 6762 | H0 not rejected |

Twelve bands with the same variants, B = 200 resamplings:

| variant | valid bands | bands | p equality cutoff | slope cutoff | p slope cutoff | p equality window | slope window | p slope window | verdict | alpha min | Teff at min | alpha max | Teff at max | statistic | dof |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| baseline | 11 | 12 | 0.005 | -0.065 | 0.010 | 0.005 | 0.015 | 0.428 | H0T rejected | 1.647 | 4538 | 2.538 | 6021 | 77.4 | 10 |
| duration at least 4 points | 11 | 12 | 0.005 | -0.127 | 0.005 | 0.005 | -0.007 | 0.741 | H0T rejected | 1.550 | 4530 | 2.216 | 6014 | 53.8 | 10 |
| duration at least 5 points | 10 | 12 | 0.005 | -0.173 | 0.005 | 0.005 | -0.038 | 0.050 | H0T rejected | 1.346 | 4753 | 1.880 | 5274 | 49.1 | 9 |
| scale by Teff only | 11 | 12 | 0.005 | 0.089 | 0.005 | 0.005 | 0.086 | 0.005 | H0T rejected | 1.650 | 4538 | 2.514 | 6021 | 104.8 | 10 |
| scale by Teff only, duration at least 4 points | 11 | 12 | 0.005 | -0.037 | 0.199 | 0.005 | 0.056 | 0.005 | H0T rejected | 1.548 | 4530 | 2.145 | 6014 | 35.1 | 10 |

Equality of the index across the bands after dropping the hottest valid band:

| bands | estimator | dropped band median Teff | dropped band tail | valid bands | statistic | dof | p bootstrap | p chi-squared |
|---|---|---|---|---|---|---|---|---|
| 12 bands | cutoff | 6021 | 1091 | 10 | 51.7 | 9 | 0.0010 | 5.07e-08 |
| 12 bands | window | 6021 | 1091 | 10 | 46.3 | 9 | 0.0010 | 5.27e-07 |
| 12 bands | pure | 6021 | 1091 | 10 | 55.2 | 9 | 0.0010 | 1.12e-08 |
| 8 bands | cutoff | 6214 | 1209 | 7 | 49.8 | 6 | 0.0050 | 5.11e-09 |
| 8 bands | window | 6214 | 1209 | 7 | 50.2 | 6 | 0.0050 | 4.28e-09 |
| 8 bands | pure | 6214 | 1209 | 7 | 46.7 | 6 | 0.0050 | 2.2e-08 |
| 16 bands | cutoff | 5926 | 1550 | 13 | 41.1 | 12 | 0.0050 | 4.8e-05 |
| 16 bands | window | 5926 | 1550 | 13 | 57.6 | 12 | 0.0050 | 6.18e-08 |
| 16 bands | pure | 5926 | 1550 | 13 | 78.9 | 12 | 0.0050 | 6.83e-12 |

Permutation test of the whole distribution of E_L above the lower limit, 2000 permutations of the labels among stars. It tests the equality of the entire distribution and not only of the index. Pairwise p values are corrected with Holm.

| test | lower limit | stat | p | P | stars | p Holm |
|---|---|---|---|---|---|---|
| classes M, K, G | 34.5 | 24.5098 | 0.0005 | 2000 | 2027 | nan |
| pair M and K | 34.5 | 6.5164 | 0.0010 | 2000 | 1518 | 0.0020 |
| pair M and G | 34.5 | 6.2183 | 0.0010 | 2000 | 938 | 0.0020 |
| pair K and G | 34.5 | 23.1030 | 0.0005 | 2000 | 1598 | 0.0015 |
| valid bands of Teff | 34.5 | 50.2443 | 0.0005 | 2000 | 2179 | nan |
| valid bands without the hottest | 34.5 | 32.4864 | 0.0005 | 2000 | 2024 | nan |
| classes M, K, G | 34.8 | 14.1508 | 0.0005 | 2000 | 1624 | nan |
| pair M and K | 34.8 | 2.4473 | 0.0100 | 2000 | 1263 | 0.0100 |
| pair M and G | 34.8 | 5.2922 | 0.0010 | 2000 | 734 | 0.0020 |
| pair K and G | 34.8 | 13.7425 | 0.0005 | 2000 | 1251 | 0.0015 |
| valid bands of Teff | 34.8 | 20.9483 | 0.0005 | 2000 | 1701 | nan |
| valid bands without the hottest | 34.8 | 17.8781 | 0.0005 | 2000 | 1622 | nan |

Holm correction over the main tests (equality by classes, equality by bands, and the comparison with 3/2 in each class):

| test | p | p Holm |
|---|---|---|
| equality of the cutoff index in M, K and G | 0.6284 | 0.6284 |
| equality of the cutoff index across the valid bands | 0.0010 | 0.0050 |
| index of class M equal to 3/2 | 0.0010 | 0.0050 |
| index of class K equal to 3/2 | 0.0010 | 0.0050 |
| index of class G equal to 3/2 | 0.0010 | 0.0050 |

Range of the cutoff index over the five windows and over all variants of the sensitivity analysis:

| class | baseline | window min | window max | all variants min | all variants max |
|---|---|---|---|---|---|
| M | 2.01 | 1.87 | 2.27 | 1.87 | 2.27 |
| K | 1.93 | 1.78 | 2.27 | 1.75 | 2.27 |
| G | 1.95 | 1.70 | 2.22 | 1.70 | 2.22 |

Duration subsets and the comparison with 3/2. Cutoff index of each class for all flares and for the flares with at least 4 or 5 points, at three lower limits, with the bootstrap by star (B = 1000) and the 99 per cent interval:

| subset | lower limit | class | tail | alpha | se | 99 low | 99 high | 3/2 rejected | p 3/2 |
|---|---|---|---|---|---|---|---|---|---|
| all flares | 34.5 | M | 14533 | 2.007 | 0.033 | 1.896 | 2.075 | True | 0.0010 |
| all flares | 34.5 | K | 25651 | 1.931 | 0.107 | 1.690 | 2.042 | True | 0.0010 |
| all flares | 34.5 | G | 8267 | 1.947 | 0.060 | 1.774 | 2.100 | True | 0.0010 |
| all flares | 34.8 | M | 7249 | 2.168 | 0.035 | 2.050 | 2.230 | True | 0.0010 |
| all flares | 34.8 | K | 13685 | 2.130 | 0.145 | 1.773 | 2.222 | True | 0.0010 |
| all flares | 34.8 | G | 3787 | 2.175 | 0.092 | 1.931 | 2.385 | True | 0.0010 |
| all flares | 35.0 | M | 4299 | 2.271 | 0.044 | 2.149 | 2.361 | True | 0.0010 |
| all flares | 35.0 | K | 8449 | 2.267 | 0.196 | 1.756 | 2.333 | True | 0.0010 |
| all flares | 35.0 | G | 2026 | 2.221 | 0.127 | 1.893 | 2.553 | True | 0.0010 |
| at least 4 points | 34.5 | M | 12553 | 1.921 | 0.042 | 1.801 | 2.009 | True | 0.0010 |
| at least 4 points | 34.5 | K | 21830 | 1.836 | 0.129 | 1.559 | 1.964 | True | 0.0010 |
| at least 4 points | 34.5 | G | 6989 | 1.797 | 0.065 | 1.613 | 1.941 | True | 0.0010 |
| at least 4 points | 34.8 | M | 6610 | 2.126 | 0.041 | 1.999 | 2.215 | True | 0.0010 |
| at least 4 points | 34.8 | K | 12428 | 2.077 | 0.163 | 1.684 | 2.172 | True | 0.0010 |
| at least 4 points | 34.8 | G | 3450 | 2.083 | 0.096 | 1.808 | 2.311 | True | 0.0010 |
| at least 4 points | 35.0 | M | 4009 | 2.251 | 0.047 | 2.115 | 2.352 | True | 0.0010 |
| at least 4 points | 35.0 | K | 7903 | 2.239 | 0.200 | 1.711 | 2.312 | True | 0.0010 |
| at least 4 points | 35.0 | G | 1906 | 2.132 | 0.133 | 1.802 | 2.461 | True | 0.0010 |
| at least 5 points | 34.5 | M | 9728 | 1.777 | 0.049 | 1.638 | 1.876 | True | 0.0010 |
| at least 5 points | 34.5 | K | 16443 | 1.698 | 0.157 | 1.362 | 1.863 | False | 0.3756 |
| at least 5 points | 34.5 | G | 5095 | 1.568 | 0.074 | 1.372 | 1.730 | False | 0.3526 |
| at least 5 points | 34.8 | M | 5590 | 2.058 | 0.049 | 1.900 | 2.146 | True | 0.0010 |
| at least 5 points | 34.8 | K | 10229 | 1.983 | 0.190 | 1.528 | 2.097 | True | 0.0020 |
| at least 5 points | 34.8 | G | 2833 | 1.893 | 0.103 | 1.614 | 2.119 | True | 0.0010 |
| at least 5 points | 35.0 | M | 3521 | 2.226 | 0.050 | 2.079 | 2.344 | True | 0.0010 |
| at least 5 points | 35.0 | K | 6840 | 2.181 | 0.225 | 1.603 | 2.266 | True | 0.0010 |
| at least 5 points | 35.0 | G | 1671 | 1.971 | 0.145 | 1.595 | 2.310 | True | 0.0040 |

Bands against 3/2 for the same subsets, with the 99 per cent interval of each valid band:

| subset | valid bands | CI above 3/2 | CI contains 3/2 | CI below 3/2 | lowest index |
|---|---|---|---|---|---|
| all flares | 11 | 11 | 0 | 0 | 1.647 |
| at least 4 points | 11 | 6 | 5 | 0 | 1.550 |
| at least 5 points | 10 | 2 | 7 | 1 | 1.346 |

| subset | band | Teff median | tail | alpha | se | 99 low | 99 high |
|---|---|---|---|---|---|---|---|
| all flares | 0 | 3317 | 6448 | 2.000 | 0.036 | 1.902 | 2.072 |
| all flares | 1 | 3646 | 6172 | 1.949 | 0.043 | 1.839 | 2.071 |
| all flares | 2 | 3900 | 5512 | 1.861 | 0.048 | 1.731 | 1.990 |
| all flares | 3 | 4308 | 4472 | 1.838 | 0.055 | 1.678 | 1.953 |
| all flares | 4 | 4538 | 5662 | 1.647 | 0.048 | 1.506 | 1.767 |
| all flares | 5 | 4820 | 4708 | 1.729 | 0.056 | 1.552 | 1.849 |
| all flares | 6 | 5105 | 4608 | 2.043 | 0.147 | 1.605 | 2.150 |
| all flares | 7 | 5313 | 4812 | 1.709 | 0.077 | 1.515 | 1.923 |
| all flares | 8 | 5522 | 3190 | 1.938 | 0.099 | 1.654 | 2.163 |
| all flares | 9 | 5767 | 2821 | 1.963 | 0.110 | 1.619 | 2.207 |
| all flares | 10 | 6021 | 1091 | 2.538 | 0.136 | 2.077 | 2.819 |
| at least 4 points | 0 | 3315 | 5497 | 1.928 | 0.040 | 1.813 | 2.011 |
| at least 4 points | 1 | 3635 | 5211 | 1.835 | 0.050 | 1.696 | 1.946 |
| at least 4 points | 2 | 3892 | 4559 | 1.765 | 0.048 | 1.634 | 1.878 |
| at least 4 points | 3 | 4298 | 4050 | 1.709 | 0.059 | 1.558 | 1.852 |
| at least 4 points | 4 | 4530 | 4836 | 1.550 | 0.058 | 1.406 | 1.679 |
| at least 4 points | 5 | 4804 | 3882 | 1.567 | 0.066 | 1.343 | 1.699 |
| at least 4 points | 6 | 5094 | 4125 | 1.980 | 0.178 | 1.478 | 2.068 |
| at least 4 points | 7 | 5308 | 4104 | 1.577 | 0.083 | 1.357 | 1.776 |
| at least 4 points | 8 | 5517 | 2645 | 1.795 | 0.108 | 1.503 | 2.026 |
| at least 4 points | 9 | 5752 | 2273 | 1.629 | 0.123 | 1.303 | 1.932 |
| at least 4 points | 10 | 6014 | 1061 | 2.216 | 0.192 | 1.722 | 2.686 |
| at least 5 points | 0 | 3310 | 4080 | 1.812 | 0.055 | 1.625 | 1.920 |
| at least 5 points | 1 | 3611 | 4020 | 1.667 | 0.053 | 1.517 | 1.788 |
| at least 5 points | 2 | 3880 | 3227 | 1.627 | 0.070 | 1.450 | 1.794 |
| at least 5 points | 3 | 4273 | 3371 | 1.555 | 0.056 | 1.415 | 1.701 |
| at least 5 points | 4 | 4519 | 3518 | 1.356 | 0.077 | 1.147 | 1.550 |
| at least 5 points | 5 | 4753 | 2754 | 1.346 | 0.068 | 1.148 | 1.478 |
| at least 5 points | 6 | 5044 | 3162 | 1.444 | 0.075 | 1.244 | 1.635 |
| at least 5 points | 7 | 5274 | 3027 | 1.880 | 0.327 | 1.061 | 1.951 |
| at least 5 points | 8 | 5480 | 2018 | 1.785 | 0.133 | 1.374 | 2.054 |
| at least 5 points | 9 | 5700 | 1737 | 1.459 | 0.154 | 1.011 | 1.793 |

Permutation test of the whole distribution with and without the scale of each group removed, 2000 permutations. The scale is the median of log E_L of the group, recomputed in every permutation. Pairwise p values are corrected with Holm.

| test | lower limit | scale | stat | p | P | stars | p Holm |
|---|---|---|---|---|---|---|---|
| classes M, K, G | 34.5 | not scaled | 24.5069 | 0.0005 | 2000 | 2027 | nan |
| pair M and K | 34.5 | not scaled | 6.5370 | 0.0005 | 2000 | 1518 | 0.0015 |
| pair M and G | 34.5 | not scaled | 6.1860 | 0.0010 | 2000 | 938 | 0.0015 |
| pair K and G | 34.5 | not scaled | 23.1045 | 0.0005 | 2000 | 1598 | 0.0015 |
| valid bands of Teff | 34.5 | not scaled | 50.2382 | 0.0005 | 2000 | 2179 | nan |
| valid bands without the hottest | 34.5 | not scaled | 32.5057 | 0.0005 | 2000 | 2024 | nan |
| classes M, K, G | 34.5 | scaled | 13.9931 | 0.0005 | 2000 | 2027 | nan |
| pair M and K | 34.5 | scaled | 3.9711 | 0.0020 | 2000 | 1518 | 0.0020 |
| pair M and G | 34.5 | scaled | 3.9049 | 0.0010 | 2000 | 938 | 0.0020 |
| pair K and G | 34.5 | scaled | 12.6592 | 0.0005 | 2000 | 1598 | 0.0015 |
| valid bands of Teff | 34.5 | scaled | 26.8775 | 0.0005 | 2000 | 2179 | nan |
| valid bands without the hottest | 34.5 | scaled | 17.9055 | 0.0005 | 2000 | 2024 | nan |
| classes M, K, G | 34.8 | not scaled | 14.1190 | 0.0005 | 2000 | 1624 | nan |
| pair M and K | 34.8 | not scaled | 2.4525 | 0.0130 | 2000 | 1263 | 0.0130 |
| pair M and G | 34.8 | not scaled | 5.2536 | 0.0010 | 2000 | 734 | 0.0020 |
| pair K and G | 34.8 | not scaled | 13.7199 | 0.0005 | 2000 | 1251 | 0.0015 |
| valid bands of Teff | 34.8 | not scaled | 20.8774 | 0.0005 | 2000 | 1701 | nan |
| valid bands without the hottest | 34.8 | not scaled | 17.8502 | 0.0005 | 2000 | 1622 | nan |
| classes M, K, G | 34.8 | scaled | 10.5071 | 0.0005 | 2000 | 1624 | nan |
| pair M and K | 34.8 | scaled | 2.4686 | 0.0075 | 2000 | 1263 | 0.0075 |
| pair M and G | 34.8 | scaled | 3.6765 | 0.0005 | 2000 | 734 | 0.0015 |
| pair K and G | 34.8 | scaled | 9.5865 | 0.0005 | 2000 | 1251 | 0.0015 |
| valid bands of Teff | 34.8 | scaled | 15.8933 | 0.0005 | 2000 | 1701 | nan |
| valid bands without the hottest | 34.8 | scaled | 13.5588 | 0.0005 | 2000 | 1622 | nan |

Simulation of the level and the power of the equality test, 1000 replicates for each scatter of the index between stars. Window estimator, three classes of 300 stars with a mean of 60 flares, bootstrap with B = 300. Power refers to an index 0.2 higher in one class.

| scatter between stars | size, bootstrap by star | size, Wald with independent flares | power, bootstrap by star | power, Wald with independent flares |
|---|---|---|---|---|
| 0.150 | 0.064 | 0.241 | 1.000 | 1.000 |
| 0.250 | 0.070 | 0.515 | 0.999 | 1.000 |
| 0.450 | 0.047 | 0.750 | 0.923 | 1.000 |

Variants with the additional checks included: 30 variants of the classes, 5 with another verdict (S5 window 34.3 to 35.3, S2 Yang class limits, S2 limits minus 200 K, S4 energy of G times 1.25, duration at least 5 points). Planned variants alone: 26, 4 with another verdict. Bands: 16 variants, 0 with another verdict, slope significant in 13, sign of the slope changes in 2. Planned variants of the bands alone: 11.


## 17. Figures

Files in the folder figures, in EPS and PDF. Figure 1 shows the complementary distribution of E_L above the common lower limit with the cutoff fits. Figure 2 shows the indices of each class for the primary analysis, the window estimate, the baseline in energy, the values of Yang and Liu and the range over the sensitivity variants. Figure 3 shows the bootstrap distributions of the differences between classes and of the test statistic. Figure 4 shows the index in the narrow bands of effective temperature for the cutoff and window estimators, with the class limits, and the median rotation period of each band. Figure 5 shows the cutoff index in the bands for three lower limits. The check of the fonts found only Courier Prime in the files.

| file | fonts | width mm |
|---|---|---|
| figure1.eps | CourierPrime-Italic-0, CourierPrime-Regular-0 | 80.00 |
| figure2.eps | CourierPrime-Italic-0, CourierPrime-Regular-0 | 80.00 |
| figure3.eps | CourierPrime-Italic-0, CourierPrime-Regular-0 | 165.00 |
| figure4.eps | CourierPrime-Italic-0, CourierPrime-Regular-0 | 80.00 |
| figure5.eps | CourierPrime-Italic-0, CourierPrime-Regular-0 | 80.00 |


## 18. Limitations found by the run

* The cutoff model is rejected by the goodness of fit test for class M (p below 0.005, KS D 0.027). The test treats flares as independent and rejects small departures when the sample is large.

* The cutoff model is rejected by the goodness of fit test for class K (p below 0.005, KS D 0.031). The test treats flares as independent and rejects small departures when the sample is large.

* The cutoff model is rejected by the goodness of fit test for class G (p below 0.005, KS D 0.017). The test treats flares as independent and rejects small departures when the sample is large.

* The cutoff model is rejected by the goodness of fit test for class F (p 0.020, KS D 0.025). The test treats flares as independent and rejects small departures when the sample is large.

* The index varies between stars of class M (between star scatter 0.22).

* The index varies between stars of class K (between star scatter 0.25).

* The index varies between stars of class G (between star scatter 0.47).

* The verdict changes in 5 variants of the sensitivity analysis.

* Bands [11] of the twelve have fewer than 1000 flares in the tail and stay out of the band tests.

* The linear trend of the cutoff index with Teff is not robust: the slope is significant in 13 of 16 variants and depends on the lower limit of E_L. The window estimate has a slope compatible with zero.

* The cutoff model trades the index against the cutoff scale, which raises the error of the index in some bands. The window estimate is more precise but is not the primary model.

* The dependence of the index on rotation is not analysed, and the rotation period distribution changes with Teff.

* The luminosity used in E_L is bolometric, from the radius and temperature of the Mathur catalog, and not the Kepler band luminosity. The typical radius uncertainty is 27 per cent.

* The catalog has no injection and recovery test, so completeness is handled through the lower limit of E_L.


## 19. Validation summary

| test | status | detail |
|---|---|---|
| alpha agrees with powerlaw at the same xmin | ok | largest difference 2.52e-13 |
| grid scan agrees with the full scan | ok | largest difference 0.0080, smallest Hill sigma 0.0241 |
| recovery of alpha without bias | ok | largest mean bias 0.0017 over 200 samples per alpha |
| truncated fit recovers alpha and the pure fit is biased | ok | bias M1 0.0202, bias M2 -0.0024 |
| goodness of fit keeps its size under a true power law | ok | 100 samples, rejection rate 0.06 at 0.1 and 0.02 at 0.05, mean p 0.59 |
| goodness of fit power against the lognormal is above its size and grows with the sample size | ok | rejection rate 0.58 at 30000 points and 0.67 at 200000 points, size 0.06. Limited power is a documented limitation of the test at the KS xmin |
| power-law data are not wrongly assigned to the lognormal | ok | 0 of 100 samples with R negative and p below 0.1 (expected near 0.1 or lower) |
| likelihood ratio power against the lognormal grows with the sample size | ok | fraction 0.49 at 30000 points and 0.72 at 200000 points |
| cutoff data are separated from the pure power law | ok | 18 of 20 samples with the nested test p below 0.1 |
| grid and full scan agree on real data | ok | largest difference in alpha 0.0068 |
| cutoff fit from per star sums recovers the index | ok | mean alpha 1.897 against 1.9, scatter 0.014, 40 samples of 20000 flares |
| equality test keeps its size with correlated flares, scatter 0.15 | ok | rejection rate 0.064 at 0.05 with the bootstrap by star against 0.241 for the independent flares Wald test, 1000 replicates |
| equality test detects a difference of 0.2, scatter 0.15 | ok | rejection rate 1.000 against 1.000 for the independent flares Wald test |
| equality test keeps its size with correlated flares, scatter 0.25 | ok | rejection rate 0.070 at 0.05 with the bootstrap by star against 0.515 for the independent flares Wald test, 1000 replicates |
| equality test detects a difference of 0.2, scatter 0.25 | ok | rejection rate 0.999 against 1.000 for the independent flares Wald test |
| equality test keeps its size with correlated flares, scatter 0.45 | ok | rejection rate 0.047 at 0.05 with the bootstrap by star against 0.750 for the independent flares Wald test, 1000 replicates |
| equality test detects a difference of 0.2, scatter 0.45 | ok | rejection rate 0.923 against 1.000 for the independent flares Wald test |
| band equality test keeps its size without a trend | ok | rejection rate 0.04 at 0.05, 100 replicates of 12 bands |
| band trend test keeps its size without a trend | ok | rejection rate 0.07 at 0.05 |
| band equality test detects a trend of 0.1 per 1000 K | ok | rejection rate 1.00, and 1.00 for the trend test |
| slope of the index in Teff is recovered | ok | mean bias 0.001 without trend and -0.001 with a trend of 0.1 per 1000 K |
