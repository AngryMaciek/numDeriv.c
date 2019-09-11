//
// numDeriv R library re-written into C
//
// Maciej Bak
// Swiss Institute of Bioinformatics
// 11.09.2019
//

// include guard:
// directives to mark that we included this header
#ifndef NUMDERIV_H
#define NUMDERIV_H

// use C's linkage and naming for these functions
// regardless of the compiler
#ifdef __cplusplus
extern "C" {
#endif

void gradient(double(*f)(unsigned short n_args, double* args, void* box), 
              unsigned short n_args, double* args, double* grad, void* box,
              double eps, double d, double zero_tol, unsigned short r, unsigned short v);

void genD(double(*f)(unsigned short n_args, double* args, void* box), 
              unsigned short n_args, double* args, double* D, void* box,
              double eps, double d, double zero_tol, unsigned short r, unsigned short v);

void hessian(double(*f)(unsigned short n_args, double* args, void* box), 
          unsigned short n_args, double* args, double* H, void* box,
          double eps, double d, double zero_tol, unsigned short r, unsigned short v);

#ifdef __cplusplus
}
#endif

#endif