## Code and data in support of: 'Mechanisms underpinning community stability along a latitudinal gradient: insights from a niche-based approach'

Luke Christopher Evans, Yolanda Melero, Reto Schmucki, Philipp H. Boersch-Supan, Lluís Brotons, Colin Fontaine, Frédéric Jiguet, Mikko Kuussaari, Dario Massimino, Robert A Robinson, David B. Roy, Oliver Schweiger, Josef Settele, Constanti Stefanescu, Chris A.M. van Turnhout, Tom Henry Oliver

Code by Luke Christopher Evans
[![DOI](https://zenodo.org/badge/495495841.svg)](https://zenodo.org/badge/latestdoi/495495841)

Note: the manuscript and the accompanying code is currently under peer review.

### Layout
The project has separate data and code folders. The master data sheet 'CommunityButterflyData.csv' is stored in MechanismsUnderpinningCommunityStability/Data, along with a list that is used for sending data to stan 'CommunityList.rds'.

The code is organised into the different stages of modelling
1) Structural Equation models
2) D separation tests
3) Correlation estimates

Each stage has one folder for the R code and one for the stan code e.g. 'Rcode Structural Equation Model' and 'Stan code Structural Equation model'

Files within the folders finish with numbers e.g. ...1.R to suggest the running order.

The code relies on the use of R stan and the rethinking package. Installation of these packages is more involved than typical R packages so follow the below links to install:
https://mc-stan.org/users/interfaces/rstan
https://github.com/rmcelreath/rethinking

note: the code was written in stan 2.21.0 and more recent versions of stan may show warnings due
to new variable declarations 
