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

/*
Calculates gradient of a given objective function by Richardson approximation
Arguments:
 - objective function (N*double -> double)
 - number of arguments of the objective function
 - array of the arguments at which gradient should be evaluated
 - array for the gradient to be updated
 - pointer to a structure with all input/parameters for the objective function
 - Richardsons: eps value for the initial approximation of arguments less then zero_tol, def=1e-4
 - Richardsons: d gives the fraction of argument to use for the initial approximation, def=1e-4
 - Richardsons: zero.tol is the tolerance for deciding which elements of arguments are zero, def=2.220446e-16/7e-7
 - Richardsons: r gives the number of Richardson improvement iterations, def=4
 - Richardsons: v is the reduction factor, def=2 
*/
void gradient(double(*f)(unsigned short n_args, double* args, void* box), 
              unsigned short n_args, double* args, double* grad, void* box,
              double eps, double d, double zero_tol, unsigned short r, unsigned short v);

/*
 Generate an array of function derivative information.
 Arguments:
 - objective function (N*double -> double)
 - number of arguments of the objective function (n_args)
 - array of the arguments at which gradient should be evaluated
 - array for the first der. + lower triangle of Hessian: size = n_args*(n_args+3)/2
 - pointer to a structure with all input/parameters for the objective function
 - Richardsons: eps value for the initial approximation of arguments less then zero_tol, def=1e-4
 - Richardsons: d gives the fraction of argument to use for the initial approximation, def=1e-4
 - Richardsons: zero.tol is the tolerance for deciding which elements of arguments are zero, def=2.220446e-16/7e-7
 - Richardsons: r gives the number of Richardson improvement iterations, def=4
 - Richardsons: v is the reduction factor, def=2 
 */
void genD(double(*f)(unsigned short n_args, double* args, void* box), 
              unsigned short n_args, double* args, double* D, void* box,
              double eps, double d, double zero_tol, unsigned short r, unsigned short v);

/*
 Generate a Hessian matrix.
 Arguments:
 - objective function (N*double -> double)
 - number of arguments of the objective function (n_args)
 - array of the arguments at which gradient should be evaluated
 - pointer to the 2D array for the Hessian matrix: size = n_args*n_args
 - pointer to a structure with all input/parameters for the objective function
 - Richardsons: eps value for the initial approximation of arguments less then zero_tol, def=1e-4
 - Richardsons: d gives the fraction of argument to use for the initial approximation, def=1e-4
 - Richardsons: zero.tol is the tolerance for deciding which elements of arguments are zero, def=2.220446e-16/7e-7
 - Richardsons: r gives the number of Richardson improvement iterations, def=4
 - Richardsons: v is the reduction factor, def=2 
 */
void hessian(double(*f)(unsigned short n_args, double* args, void* box), 
          unsigned short n_args, double* args, double* H, void* box,
          double eps, double d, double zero_tol, unsigned short r, unsigned short v);

#ifdef __cplusplus
}
#endif

#endif