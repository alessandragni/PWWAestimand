# PWWAestimand: Patient-Weighted While-Alive Estimand


This repository contains R demo scripts for the estimation 
of the Patient-Weighted While-Alive Estimand, 
as proposed by Schmidli et al. (2023) [1], within a randomized treatment setting. 

The implementation builds on theoretical developments from

> **Nonparametric estimation of the Patient Weighted While-Alive Estimand**  
> [arXiv preprint](https://arxiv.org/abs/2412.03246)
> by A. Ragni, T. Martinussen, T. Scheike  

and leverages the implementation included within the R package `mets` 
(specifically `WA_recurrent()`, by Thomas Scheike).

The development version of the package can be installed through
```R
remotes::install_github("kkholst/mets",ref="develop")
```

## data folder
This folder includes the colorectal cancer data [2], used in the paper
in the Web Appendix D.

The HF-ACTION study data [3] analyzed in Section 6 of the paper, 
cannot be shared due to privacy and confidentiality considerations.
For our paper, they were kindly made available 
by the BioLINCC of the National Heart, Lung, and Blood Institute.


## R folder
This folder provides R scripts to reproduce the analyses from Sections 5 and 6 of the paper:

* The script `DEMO_CaseStudy.R` contains the code to estimate the PWWA Estimand
in the colorectal cancer case study setting, whose results are presented in the Web Appendix D of the paper.
Also includes code to compare the PWWA estimand with the EWWA estimand [4], and to generate relevant plots.

* The script `DEMO_SimulationStudy.R` performs simulations under various settings to assess PWWA estimator performance,
and compare it to the EWWA estimator.
Computes key performance metrics such as mean, bias, standard error, coverage, and power.
Although based here on Copenhagen hospital data, 
similar simulations produced Tables 1 and 2 and Web Table 4 in the paper using HF-ACTION data.

* The script `utils.R` provides functions useful for the DEMO scripts.

* The scripts `ProduceTab1.R`, `ProduceTab2.R`, `ProduceWebTab4.R` respectively produce Tables 1, 2 and Web Table 4 in the paper,
from the outputs of the simulations using HF-ACTION data, that are available at the following 
[link](https://polimi365-my.sharepoint.com/:f:/g/personal/10500163_polimi_it/Et-1kZPX-chDviMUAh2qxXoBXfh1xRw-OrEyZKcQAGhC8A?e=qRFJ1D)


## References
[1] Schmidli, Heinz, James H. Roger, and Mouna Akacha. "Estimands for recurrent event endpoints in the presence of a terminal event." Statistics in Biopharmaceutical Research 15.2 (2023): 238-248.

[2] Ducreux, Michel, et al. "Sequential versus combination chemotherapy for the treatment of advanced colorectal cancer (FFCD 2000–05): an open-label, randomised, phase 3 trial." The lancet oncology 12.11 (2011): 1032-1044

[3] O’Connor, Christopher M., et al. "Efficacy and safety of exercise training in patients 
with chronic heart failure: HF-ACTION randomized controlled trial." Jama 301.14 (2009): 1439-1450.

[4] Mao, Lu. "Nonparametric inference of general while‐alive estimands for recurrent events." Biometrics 79.3 (2023): 1749-1760.

