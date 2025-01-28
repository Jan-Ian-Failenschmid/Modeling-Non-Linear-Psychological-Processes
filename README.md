<!---
# Filename: /home/jan_failenschmid/research/project_1/paper_1/README.txt
# Created Date: Thursday, October 10th 2024, 12:57:24 pm
# Author: Jan Ian Failenschmid
# 
# Copyright (c) 2024 by Jan Ian Failenschmid
# E-mail: J.I.Failenschmid@tilburguniveristy.edu
# 
# License: GNU General Public License v3.0 or later
# License URL: https://www.gnu.org/licenses/gpl-3.0-standalone.html
# License File: gpl-3.0.txt
#
-->

# Modeling Non-Linear Psychological Processes: Reviewing and Evaluating Non-parametric Approaches and Their Applicability to Intensive Longitudinal Data

This directory contains the suplimentary material to the article titled: 
"Modeling Non-Linear Psychological Processes: Reviewing and Evaluating 
Non-parametric Approaches and Their Applicability to Intensive Longitudinal Data".

[Preprint](https://osf.io/preprints/psyarxiv/26mde)

[OSF Directoy](https://osf.io/3ytdx/wiki/home/)


## Abstract 

Psychological concepts are increasingly understood as complex dynamic systems that change
over time. To study these complex systems, researchers are increasingly gathering intensive
longitudinal data (ILD), revealing non-linear phenomena such as asymptotic growth, mean-level
switching, and regulatory oscillations. However, psychological researchers currently lack
advanced statistical methods that are flexible enough to capture these non-linear processes
accurately, which hinders theory development. While methods such as local polynomial
regression, Gaussian processes, and generalized additive models (GAMs) exist outside of
psychology, they are rarely applied within the field because they have not yet been reviewed
accessibly and evaluated within the context of ILD. To address this important gap, this article
introduces these three methods for an applied psychological audience. We further conducted a
simulation study, which demonstrates that all three methods infer non-linear processes that have
been found in ILD more accurately than polynomial regression. Particularly, GAMs closely
captured the underlying processes, performing almost as well as the data-generating parametric
models. Finally, we illustrate how GAMs can be applied to explore idiographic processes and
identify potential phenomena in ILD. This comprehensive analysis empowers psychological
researchers to model non-linear processes accurately and select a method that aligns with their
data and research goals.

## Running Web App

[Shiny APP](https://jan-ian-failenschmid.shinyapps.io/modeling_non_linearity_app/)

## File Overview

Here is an overview of the files in this directory:

- main.txt is the main latex file that was used to build the article.

- apa7.cls contains the latex class necessary for building apa7 journal 
articles. 

- bibliography.bib and R_bib.bib are bibliographies that contain all used 
journal articles and R packages. 

- the R subdirectory includes all the R-code files that were used to simulate
the data, perform the analyses, generate the tables and figures that are
described in the article. This directory contains the main-files to replicate 
the analyses presented in the article. Each of the code files is writen 
to be sourceable (e.g., they can be executed top to bottom to replicate 
the analyses and results presented in the article).
    - simulation.R runs the simulation and saves the simulated data and results 
    into the data subdirectory. 
    - analysis.R runs the analysis on the simulation results, generates 
    the results summary figures and saves them in the figures subdirectory, and 
    it runs the associated ANOVAs and saves the effect sizes in the tables 
    subdirectory. 
    - demonstration.R runs the analysis for the empirical example. This is the 
    only code file that cannot be sourced, without adding the code files from 
    the Leuven Clinical Study to the data subdirectory
    - exemplar_plot.R generates the figures for the exampler processes and 
    saves them in the figures subdirectory. 
    - method_illustration.R generates the figures that illustrate how each 
    method works and saves them in the figures subdirectory. 
    - sessionInfo.txt contains replicability information of the exact R 
    computing environment that was used for the simulation and analysis. 
    - helper subdirectory contains R files with class, method, and function 
    definitions that are sourced by the main R files. 
    - stan_files subdirectory contains stan files and executables to run 
    the GP regression. (Executables may need to be recompiled on different 
    machines)
    - data subdirectory contains all simulated data and results

- figures directory contains all figures that are imported into the 
article.

- tables directory contains tex files for all tables that have been 
copied into main.tex as per journal requirements. 

## Any comments or feedback is very welcome!
