# information-theoretic-sdfs
R Code to extract "information theoretic" SDFs from Asset Pricing Theory

This repository is designed to enable students (and practitioners) to extract "information theoretic" stochastic discount factors (SDFs) from asset pricing theory.
There are two key functions - the exponential_tilting_fn.R and the exponential_tilting_grad.R which are used in the main script file.  
The code leverages the amazing lgfgs package in R which is able to solve optimizations with an l1-penalty quickly and efficiently.
The repository is designed so that the user can simply execute the extract_sdf.R script file and the SDF will be extracted from there.

The current iteration of the code extracts the SDF using the exponential tilting method of Kitamura and Stutzer (1997).  This method is operationlised empirically by 
Julliard, Ghosh and Taylor (2019) and is augmented with an l1-penalty by Qui and Otsu (2020).  The code allows for the l1-penalty.

If anyone has any suggestions/comments, please feel free to contact me.
