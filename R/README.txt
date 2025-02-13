This directory contains the suplimentary material to the article:

Modeling Non-Linear Psychological Processes: Reviewing and Evaluating
  Non-parametric Approaches and Their Applicability to Intensive
  Longitudinal Data

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