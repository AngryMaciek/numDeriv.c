# numDeriv.c: Reimplementation of the numDeriv R library into C language
*Maciej Bak*  
*Swiss Institute of Bioinformatics*

This repository contains a reimplementation of the well-known numDeriv package from R into C.
Created as a single-file library to be included in other C/C++ programs; dedicated to scientific computing.
Although the implementation is based on numDeriv version: 2016.8-1.1 the its functionality is limitted:  
* Precise calculations are done using Richardson's extrapolation only. The parameters for Richardson's method are explicit arguments for the library functions but the accuracy of the approximations depends on the objective function and may vary for a given set of parameters values. Therefore it is strongly advised to always simulate a dataset and adjust the parameters to minimize the error for every model being built.  
* Side derrivatives are not implemented yet.  
* Objective functions can have a vector argument but scalar value is a requirement: f:R^N -> R. Therefore jacobian function is not implemented.  


Original numDeriv R package:  
[https://cran.r-project.org/web/packages/numDeriv/index.html][link1]  

## Compilation
The library has been tested on macOS 10.13.6 with gcc version 4.2.1 and Ubuntu 14.04.4 with gcc 4.8.4
```
$ gcc tests.c -std=c99 -o exe -lm
$ ./exe
```
The results have been compared to exact solutions as well as approximations calculated in R.

## Repository
This repository contains four files:

| File | Description |
| ------ | ------ |
| README.md | (this file) |
| numDeriv.c | Main files: gradient, hessian, genD |
| tests.c | Source file with all the tests |
| LICENSE | The GNU General Public License v3.0 |

## License
GNU General Public License 


[link1]: <https://cran.r-project.org/web/packages/numDeriv/index.html>