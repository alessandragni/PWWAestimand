# PWWAestimand: Patient-Weighted While-Alive Estimand

Github repository that provides demo scripts implemented in R
allowing for the estimation of the Patient-Weighted While-Alive Estimand, 
as suggested in Schmidli et al. (2023) [1],
within a randomized treatment setting, 
relying on the theory developed in https://doi.org/10.48550/arXiv.2412.03246
and the implementation included within the R package `mets`.

## data folder
In the folder data, we provide the colorectal cancer data [2], employed in the paper
in the Web Appendix D.

The HF-ACTION study data [3] analyzed in Section 6 of the paper,
that were kindly made available 
by the Biologic Specimen and Data Repository Information Coordinating Center (BioLINCC) 
of the National Heart, Lung, and Blood Institute
cannot be shared due to privacy and confidentiality considerations.


## R folder
In the folder R, we provide the R scripts that 
would allow to reproduce the results of the paper (Sections 5 and 6) if the data from HF-ACTION were available.

More specifically, the script `CaseStudy_DEMO.R` contains the code to estimate the Patient-Weighted While-Alive Estimand
in the colorectal cancer case study setting, whose results are presented in the Web Appendix D of the paper.
The script allows to produce plots to compare the PWWA estimand with the EWWA estimand [4].

The script `SimulationStudy_DEMO.R` contains the code to estimate the Patient-Weighted While-Alive Estimand
within different simulated setting, allowing to compare different estimators 
and computing the main indicators, such as mean, bias, standard error, coverage and power.
The simulation is based on patents data from a Copenhagen hospital; 
in the paper, Tables 1 and 2 were obtained through similar code applied to the HF-ACTION data instead.


## References
[1] Schmidli, Heinz, James H. Roger, and Mouna Akacha. "Estimands for recurrent event endpoints in the presence of a terminal event." Statistics in Biopharmaceutical Research 15.2 (2023): 238-248.

[2] Ducreux, Michel, et al. "Sequential versus combination chemotherapy for the treatment of advanced colorectal cancer (FFCD 2000–05): an open-label, randomised, phase 3 trial." The lancet oncology 12.11 (2011): 1032-1044

[3] O’Connor, Christopher M., et al. "Efficacy and safety of exercise training in patients 
with chronic heart failure: HF-ACTION randomized controlled trial." Jama 301.14 (2009): 1439-1450.

[4] Mao, Lu. "Nonparametric inference of general while‐alive estimands for recurrent events." Biometrics 79.3 (2023): 1749-1760.

