# PSPM-Disease

Program and data repository for:

**Pathogens stabilize or destabilize depending on host stage structure**

Jessica Hite & Andre M. de Roos

Mathematical Biosciences and Engineering 20(11): XXXX-XXXX

## Repository contents

- **EBT/PSPM-Disease.c** \
Model implementation for numerical simulation of the size-structured consumer-resource model with SI dynamics using the *Escalator Boxcar Train* (EBT) approach. See the [EBTtool webpage](https://staff.fnwi.uva.nl/a.m.deroos/EBT/index.html) for a description and the EBT software package.

- **EBT/PSPM-Disease.h** \
Header file defining the dimensions of the model, including the number of environmental variables, the number of populations and the number of state variables characterising individuals

- **EBT/Makefile** \
  Simple Makefile to produce the executables `PSPM-Disease` (for computing model dynamics) and `PSPM-Disease_bif` (for performing numerical simulations for bifurcation analysis)

- **EBT/Default.cvf** \
  File with default values of the model parameters and the numerical settings needed for simulations.

- **EBT/Default.isf** \
  File with the default initial state of the model for simulations.

- **Equi/PSPM-Disease.h** \
  Model implementation for computing ecological steady states of the model using the 
  R package [FindCurve](https://doi.org/10.5281/zenodo.5642759).

- **Figures/Figure2.R** \
  R script to produce Figure 2 in the main text of the manuscript

- **Figures/Figure3.R** \
  R script to produce Figure 3 in the main text of the manuscript

- **Figures/Figure4.R** \
  R script to produce Figure 4 in the main text of the manuscript

- **Figures/Figure5.R** \
  R script to produce Figure 5 in the main text of the manuscript

- **Figures/Figure6.R** \
  R script to produce Figure 4 in the main text of the manuscript

- **Figures/Figure7.R** \
  R script to produce Figure 5 in the main text of the manuscript

- **Figures/Figure1.pdf** \
  PDF file of Figure 1 in the main text of the manuscript

- **Figures/Figure2.pdf** \
  PDF file of Figure 2 in the main text of the manuscript

- **Figures/Figure3.pdf** \
  PDF file of Figure 3 in the main text of the manuscript

- **Figures/Figure4.pdf** \
  PDF file of Figure 4 in the main text of the manuscript

- **Figures/Figure5.pdf** \
  PDF file of Figure 5 in the main text of the manuscript

- **Figures/Figure6.pdf** \
  PDF file of Figure 6 in the main text of the manuscript

- **Figures/Figure7.pdf** \
  PDF file of Figure 7 in the main text of the manuscript

- **Figures/EBToutput/** \
  This directory contains the files with the values of the model parameters and numerical settings (files with `.cvf` extension) and the initial state of the model (files with `.isf` extension) that are used by the EBT programs `PSPM-Disease`, for computing results in Figure 7 in the main text and by `PSPM-Disease_bif`, for computing the bifurcation results (minimum and maximum values during cylces) in Figure 4, 5, and 6. Because these model simulations take quite a long time also the all the output files are included (files with `.out` extension) that are generated during these model simulations and that are used for producing the figures included in the manuscript.

