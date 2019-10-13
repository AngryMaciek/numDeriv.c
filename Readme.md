# numDeriv.c: Reimplementation of the numDeriv R library into C language
*Maciej Bak*  
*Swiss Institute of Bioinformatics*

This repository contains a reimplementation of the well-known numDeriv package (dedicated to scientific computing) from R into C static library that can be included in C or C++ code. Although the implementation is based on numDeriv version: 2016.8-1.1 the functionality of this library is limited:  
* Precise calculations are done using Richardson's extrapolation only. The parameters for Richardson's method are explicit arguments for the library functions but the accuracy of the approximations depends on the objective function and may vary for a given set of parameters values. Therefore it is strongly advised to always simulate a dataset, adjust the parameters and rebuild the library in order to minimize the error for every model being built.  
* Side derrivatives are not implemented yet.  
* Objective functions can have a vector argument but scalar value is a requirement: f:R^N -> R. Therefore jacobian function is not implemented.  

Original numDeriv R package:  
[https://cran.r-project.org/web/packages/numDeriv/index.html][link1]  

## Installation

#### Quick and dirty:
Just copy the library header (numDeriv.h) and implementation (numDeriv.c) files into the directory with your source code and include the header in your main compilation unit as: `#include "numDeriv.h"`. As the library requires math.h it is necessary to provide `-lm` to the linker.  Also, the library is written in C99 standard, therefore compile as:
```bash
$ {COMPILER} {SOURCE CODE} numDeriv.c -o {OUTFILE} -std=c99 -lm
```

#### Fancy and elegant:
To build and install a static library on macOS and Linux systems with the prepared Makefile one would require `make` and `gcc` installed. 

 1. Run `./configure` to generate a full Makefile. By default the library will be installed under `~/lib` (compiled file under `~/lib/lib`, header file under `~/lib/include`). The prefix path might be adjusted in the configure file if needed.
 2. Run `make` to build a static library.
 3. Run `make install` to copy the files into `~/lib`.
 4. Run `make clean` to remove the build files from the working directory.  
 
In order to remove the library run: `make uninstall`.

Once installed the library can be tested by :
```bash
$ cd test
$ gcc tests.c -I ~/lib/include/ ~/lib/lib/numDeriv.a -o exe -std=c99 -lm
$ ./exe
```
In general, while compiling another source code that utilizes numDeriv functions one has to include the header file as `<#include numDeriv.h>`  (and provide information to the linker about its path) and explicitly state the full path to the static library.

The library has been tested on macOS 10.14.6 with clang version 4.0.1 and Ubuntu 14.04.4 with gcc 4.8.4. The results have been compared to exact solutions as well as approximations calculated in R.

## License
GNU General Public License v3.0

[link1]: <https://cran.r-project.org/web/packages/numDeriv/index.html>
