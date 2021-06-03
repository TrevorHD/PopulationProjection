# Overview

This code projects population growth given a transition matrix **A** and initial abundances of each age/stage/class. Graphical outputs include graphs of abundance and proportional abundance. Numeric outputs include abundances and proportional abundances at the final time step, the stable age distribution (right eigenvector of **A**), and the stable growth rate (dominant eigenvalue of **A**). Raw outputs of abundance and proportional abundance at each time step can be called as well.

<br/>

# Files

## Scripts

**PopProjection** *(.R)* - R code used to specify initial conditions, run the matrix models, print output, and generate simple visualisations of total and proportional abundance.

## Figures

**ExAbundance** *(.jpeg)* - Example output showing abundance for each age/stage/class as a function of time.

**ExProportions** *(.jpeg)* - Example output showing proportional abundance for each age/stage/class as a function of time.
