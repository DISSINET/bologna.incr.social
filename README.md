# bologna.incr.social
R syntax of the Gender, kinship, and other social predictors of incrimination in the inquisition register of Bologna (1291–1310): Results from an exponential random graph model manuscript

## this repo
contains the R-syntax for analyzing the social patterns of the inquisition register of Bologna (1291–1310).

## related paper
is available at [this hyperlink](https://www.dissinet.cz/). 
The two tsv datasets used in the analysis are in the `data` subdirectory.

## citation
If you use this code in your research, please cite this repo and/or the paper as:
...

## dependencies
The code is written in R 4.3.1.
Following packages are required:
* parallel 4.3.1 (R Core Team (2023), <https://stat.ethz.ch/R-manual/R-devel/library/parallel/html/00Index.html>.)
* Matrix 1.6.1 (Bates D, Maechler M, Jagan M (2023). _Matrix: Sparse and Dense Matrix Classes and Methods_. R package version 1.6-1,
  <https://CRAN.R-project.org/package=Matrix>.)
* igraph 1.5.1 (Csardi G, Nepusz T (2006). The igraph software package for complex network research. _InterJournal_, *Complex Systems*, 1695. <https://igraph.org>.)
* netUtils 0.8.2 (David Schoch (2023). netUtils: A Collection of Tools for Network Analysis <https://CRAN.R-project.org/package=netUtils>) 
* netseq 1.0.2 (Bojanowski M (2021). _Measures of Network Segregation and Homophily_. R package version 1.0-1, <https://mbojan.github.io/netseg/>)
* poweRlaw 0.70.6 (Colin S. Gillespie (2015). Fitting Heavy Tailed Distributions: The poweRlaw Package. Journal of Statistical Software, 64(2), 1-16. URL
  <http://www.jstatsoft.org/v64/i02/>)
* network 1.18.1 (Butts C (2015). _network: Classes for Relational Data_. The Statnet Project (<http://www.statnet.org>). R package version 1.13.0.1)
* ergm 4.5.0 (Handcock MS, Hunter DR, Butts CT, Goodreau SM, Krivitsky PN, Morris M (2023). _ergm: Fit, Simulate and Diagnose Exponential-Family Models for Networks_.
  The Statnet Project (<https://statnet.org>). R package version 4.5.0, <https://CRAN.R-project.org/package=ergm>)
* statnet.common 4.9.0 (Krivitsky PN (2023). _statnet.common: Common R Scripts and Utilities Used by the Statnet Project Software_. The Statnet Project (<https://statnet.org>).
  R package version 4.9.0, <https://CRAN.R-project.org/package=statnet.common>)  
* sna 2.7-1 (  Butts CT (2023). _sna: Tools for Social Network Analysis_. R package version 2.7-1, <https://CRAN.R-project.org/package=sna> )
* tergm 4.2.0 (Krivitsky PN, Handcock MS (2023). _tergm: Fit, Simulate and Diagnose Models for Network Evolution Based on Exponential-Family Random Graph Models_. The
  Statnet Project (<https://statnet.org>). R package version 4.2.0, <https://CRAN.R-project.org/package=tergm>)
* networkDynamic 0.11.3 (Butts C, Leslie-Cook A, Krivitsky P, Bender-deMoll S (2023). _networkDynamic: Dynamic Extensions for Network Objects_. R package version 0.11.3,
  <https://CRAN.R-project.org/package=networkDynamic>)
* intergraph 2.0-3 (Bojanowski M (2023). _intergraph: Coercion Routines for Network Data Objects_. R package version 2.0-3, <https://mbojan.github.io/intergraph/>)  
* ergMargins 0.1.3.1 ( Duxbury S (2023). _ergMargins: Process Analysis for Exponential Random Graph Models_. R package version 0.1.3.1,
  <https://CRAN.R-project.org/package=ergMargins>)
