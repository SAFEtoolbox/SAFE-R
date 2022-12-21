# SAFE-R
R version of the Sensitivity Analysis For Everybody (SAFE) toolbox

### BEFORE STARTING

An introduction to the SAFE Toolbox is provided in the paper:

Pianosi, F., Sarrazin, F., Wagener, T. (2015), A Matlab toolbox for Global Sensitivity Analysis, Environmental Modelling & Software, 70, 80-85. The paper is freely available at: www.sciencedirect.com/science/article/pii/S1364815215001188

We recommend reading this (short) paper before getting started. Other reading materials, including general introductions to Sensitivity Analysis and case study applications, can be found at: www.safetoolbox.info

### INSTALLING SAFE

Before installing SAFE, you must install these other packages (unless you have installed them already):
- caTools
- ggplot2
- cowplot
- gridExtra
- FRACTION
- matrixStats
- devtools
- calibrater 

The first six are available on CRAN, whereas calibrater is available at: https://people.maths.bris.ac.uk/~mazjcr/calibrater_0.51.tar.gz

To install a package from CRAN, open the R command line and type:

    install.packages("caTools")

To install a package from a web page, use:

    install.packages("http://www.maths.bris.ac.uk/~mazjcr/calibrater_0.51.tar.gz", repos = NULL, type = "source")

Once you have installed the above packages, install SAFE with the commands:

    library(devtools)
    install_github("FrancescaPianosi/SAFE-R")
 
### GETTING STARTED

To get started using SAFE, we suggest opening one of the workflow scripts in the 'demo' folder and running the code step by step. The header of each workflow script gives a short description of the method and case study model, and of the main steps and purposes of that workflow, as well as references for further reading. The name of each workflow is composed as: workflow_method_model

Implemented models are:
- the hydrological Hymod model 
- the hydrological HBV model 
- the Ishigami and Homma test function 
- the Sobol' g-function 

Implemented methods are:
- eet (elementary effects test, or method of Morris)
- fast (Fourier amplitude sensitivity test)
- pawn
- rsa (regional sensitivity analysis)
- vbsa (variance-based sensitivity analysis, or method of Sobol')

If the user still does not have a clear idea of what method(s) to start with, we suggest one of the three most widely used methods: eet (e.g. workflow_eet_hymod), rsa (workflow_rsa_hymod) or vbsa (workflow_vbsa_hymod).

