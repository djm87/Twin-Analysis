# Twin-Analysis
A twin analysis code built on MTEX and Matlab's graph toolbox

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://github.com/djm87/Twin-Analysis/blob/master/LICENSE)


## Summary
To perform quantitative twin analysis each twinned grain fragment reconstructed from EBSD scans needs to have properties of twin type, generation, and the parent grain (pretwinned). For small volume fractions of twins, this is a trivial task; however, for highly twinned or deformed microstructures the task becomes non-trival.

The objective of this project is to recreate and improve on the descriptions of twin analysis presented in [graph theory](https://link.springer.com/article/10.1007/s40192-018-0106-y) and [hierarchy](https://onlinelibrary.wiley.com/doi/full/10.1111/j.1365-2818.2009.03343.x) based twin analysis methodologies using an open source toolset. MTEX was chosen for its extensive high level toolsets in orientation analysis, EBSD, and texture. Unfortunately this also means you need Matlab to use the source code provided here.

Please cite our recent paper on hiearchial twinning in Ti if you use this code.
https://doi.org/10.1016/j.matchar.2020.110808 

## Getting started
Look at the [wiki](https://github.com/djm87/Twin-Analysis/wiki) for examples and technical documentation.

## Compatibility 
The source code is known to work with Matlab 2020b and MTEX version 5.8.
