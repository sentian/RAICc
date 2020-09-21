### Selection of Regression Models under Linear Restrictions for Fixed and Random Designs

This project studies the use of information criteria, either squared-error based (e.g Cp, RCp, FPE and Sp) or Kullback-Leibler (KL) based (e.g. AICc), on variable selection and general restriction problems, under either fixed-X or random-X design. A novel KL-based criterion for random-X, RAICc is proposed and justified. KL-based criteria (AICc and RAICc) provide better predictive performance compared to the squared-error based criteria and cross-validation. For more details, see our paper [Tian, S., Hurvich, C. and Simonoff, J. (2020): "Selection of Regression Models under Linear Restrictions for Fixed and Random Designs"](https://github.com/sentian/RAICc/blob/master/paper/ms.pdf).

This repository provides the code to reproduce the simulation results and the figures in the paper. The complete set of simulation results is presented at the [supplemental material](https://github.com/sentian/RAICc/blob/master/paper/simuresults-compressed.pdf)).

The structure of the *'code'* directory is shown below. *'plots.R'* creates the figures in the paper and the supplemental material, which are saved in the folder *paper/figures*. it relies on the simulation results in the folder *'code/run_model/results'*. *'.sh'* are bash files that submit the corresponding R code to a Linux server for running in parallels.

```
├── plots_tables
│   └── plots.R                         ## code that generates the figures in the paper
├── run_model
│   ├── para_forhpc                     ## parameters for the simulation configurations
│   ├── results
│   │   ├── fixedx                      ## simulation results for fixed-X
│   │   └── randomx                     ## simulation results for random-X
│   ├── run.R                           ## code that fits and evaluates the methods
│   ├── run_n1000_fixedx.sh             ## bash file for HPC submission
│   ├── run_n1000_randomx.sh            ## bash file for HPC submission
│   ├── run_n10n40n200_fixedx.sh        ## bash file for HPC submission
│   └── run_n10n40n200_randomx.sh       ## bash file for HPC submission
└── utils.R                             ## code that contains all the functions required 
```
