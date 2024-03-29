% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SAFER-package.R
\docType{package}
\name{SAFER-package}
\alias{SAFER-package}
\alias{SAFER}
\title{Sensitivity analysis for everybody with R}
\description{
Range of tools for Global Sensitivity Analysis (GSA).
}
\details{
It implements several established GSA methods, including method of Morris, regional sensitivity analysis, variance-based sensitivity analysis (Sobol') and FAST. It also includes new approaches and visualization tools to complement these established methods.

\code{SAFER} is also available in Matlab/Octave version through a Toolbox called \code{SAFE}.
\code{SAFE} and \code{SAFER} are open source and freely available from the following website: \url{http://bristol.ac.uk/cabot/resources/safe-toolbox/}.

BEFORE STARTING

An introduction to the SAFE Toolbox is provided in the paper: Pianosi, F., Sarrazin, F., Wagener, T. (2015), A Matlab toolbox for Global Sensitivity Analysis, Environmental Modelling & Software, 70, 80-85. The paper is freely available at: \url{http://www.sciencedirect.com/science/article/pii/S1364815215001188}.
We recommend reading this (short) paper before getting started.

TO GET STARTED

To get started using SAFER, we suggest opening one of the workflow scripts and running the code step by step. The header of each workflow script gives a short description of the method and case study model, and of the main steps and purposes of that workflow, as well as references for further reading. The name of each workflow is composed as: workflow_<method>_<model>

Implemented models are:
\itemize{
\item the hydrological Hymod model (see documentation and help of functions \code{help("hymod")} and \code{system.file("docs", "Hymod_structure.pdf", package = "SAFER")})
\item the hydrological HBV model (see documentation and help of functions in \code{help("hbv")} and \code{system.file("docs", "HBV_structure.pdf", package = "SAFER")})
\item the Ishigami and Homma test function (see help of functions in \code{help("ishigami_homma")}
\item the Sobol' g-function (see help of functions in \code{help("sobol_g_function")}
}

Implemented methods are:
\itemize{
\item eet (elementary effects test, or method of Morris)
\item fast (Fourier amplitude sensitivity test)
\item rsa (regional sensitivity analysis)
\item vbsa (variance-based sensitivity analysis, or method of Sobol')
\item pawn (distribution-based sensitivity analysis)
}



There are a 9 workflows in the demo folder of the package:
\itemize{
  \item \code{workflow_eet_hbv} This script provides an application example of the  Elementary Effects Test to the HBV rainfall-runoff model.
  \item \code{workflow_eet_hymod} This script provides a basic application example of the Elementary Effects Test. The application example is the rainfall-runoff Hymod model.
  \item \code{workflow_fast_gsobol} This script provides an application example of the Fourier Amplitude Sensitivity Test (FAST). The application example is the Sobol g-function.
  \item \code{workflow_fast_hymod} This script provides an application example of the Fourier Amplitude Sensitivity Test (FAST). The application example is the rainfall-runoff Hymod model.
  \item \code{workflow_rsa_hymod} This script provides an application example of regional sensitivity analysis. The application example is the rainfall-runoff Hymod model.
  \item \code{workflow_vbsa_hymod} This script provides an application example of Variance Based Sensitivity Analysis (VBSA) The application example is the rainfall-runoff Hymod model
  \item \code{workflow_vbsa_ishigami_homma}  This script applies Variance-Based Sensitivity Analysis to the Ishigami-Homma function.
  \item \code{workflow_PAWN_ishigami_homma}	 This script applies the PAWN distribution-based Sensitivity Analysis method to the Ishigami-Homma function.
  \item \code{workflow_PAWN_hymod}	 This script provides an application example of the PAWN distribution-based Sensitivity Analysis method. The application example is the rainfall-runoff Hymod model.
}

The command to find the list of demo in the package is
\code{demo(package = "SAFER")}

The command to start a demo is
\code{demo("workflow_visual_ishigami_homma")}

The command to find the location of the demo folder is
\code{system.file("demo", package = "SAFER")}

\code{SAFER} depends on the package \code{calibrater} that is available on Jonty Rougier's webpage
\url{http://www.maths.bris.ac.uk/~mazjcr/#software}.
}
