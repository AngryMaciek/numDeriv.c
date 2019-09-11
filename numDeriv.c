//
// numDeriv R library re-written into C
//
// Maciej Bak
// Swiss Institute of Bioinformatics
// 11.09.2019
//

// These are the definitions of library functions.
// Test the code with compilation into an object file:
// $ gcc numDeriv.c -c -std=c99


// This file defines the interface declared in the following header.
// It is good to put that information here to avoid potential errors
// (eg. inconsistency between .h and .c)
#include "numDeriv.h"


// Headers from C standard library that I need
// for the implementation of numDeriv functions:
// (they are protected by include guards)
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


void gradient(double(*f)(unsigned short n_args, double* args, void* box), 
              unsigned short n_args, double* args, double* grad, void* box,
              double eps, double d, double zero_tol, unsigned short r, unsigned short v){

  // create the 2D array of doubles
  unsigned short len = sizeof(double *) * r + sizeof(double) * n_args * r; 
  double **a = (double **)malloc(len); 
  // ptr is now pointing to the first element in of 2D array 
  double *ptr = (double *)(a + r); //int?
  // for loop to point rows pointer to appropriate location in 2D array 
  for(unsigned short i = 0; i < r; i+=1) a[i] = (ptr + n_args * i); 
  
  // first order derivatives are stored in the matrix a[k,i], 
  // where the indexing variables k for rows, i for columns,
  // r is the number of iterations, and n_args is the number of variables.

  double* h = (double *) malloc(sizeof(double) * n_args);
  for (unsigned short i=0; i<n_args; i+=1) {
    h[i] = fabs(d*args[i]) + eps * ((fabs(args[i])<zero_tol) ? 1 : 0); // ? has lower order than * (!)
  }

  double* move_plus = (double *) malloc(sizeof(double) * n_args);
  double* move_minus = (double *) malloc(sizeof(double) * n_args);
  for (unsigned short k=0; k<r; k+=1) { // successively reduce h	
    for (unsigned short i=0; i<n_args; i+=1) {
      if ((k!=0) && (fabs(a[k-1][i])<1e-20)){
        a[k][i] = 0; // some func are unstable near zero
      }
      else{
        for (unsigned short xx=0; xx<n_args; xx+=1) {
          move_plus[xx] = args[xx];
          move_minus[xx] = args[xx];
        }
        move_plus[i] += h[i];
        move_minus[i] -= h[i];
        a[k][i] = ((*f)(n_args, move_plus, box) - (*f)(n_args, move_minus, box)) / (2*h[i]);
      }
    }
    
    for (unsigned short i=0; i<n_args; i+=1) {
      h[i] = h[i]/v; // Reduced h by 1/v.
    }
  }
  
  free(move_plus);
  free(move_minus);
  free(h);
  
  /*
  #------------------------------------------------------------------------
  # 1 Applying Richardson Extrapolation to improve the accuracy of 
  #   the first and second order derivatives. The algorithm as follows:
  #
  #   --  For each column of the derivative matrix a,
  #	  say, A1, A2, ..., Ar, by Richardson Extrapolation, to calculate a
  #	  new sequence of approximations B1, B2, ..., Br used the formula
  #
  #	     B(i) =( A(i+1)*4^m - A(i) ) / (4^m - 1) ,  i=1,2,...,r-m
  #
  #		N.B. This formula assumes v=2.
  #
  #   -- Initially m is taken as 1  and then the process is repeated 
  #	 restarting with the latest improved values and increasing the 
  #	 value of m by one each until m equals r-1
  #
  #  Return the final improved  derivative vector.
  #-------------------------------------------------------------------------
  */
  
  for (unsigned short m=0; m<r-1; m+=1){
    for (unsigned short i = 0; i < r-m-1; i+=1){
      for (unsigned short j = 0; j < n_args; j+=1){
        a[i][j] = (a[i+1][j] * pow(4.0,m+1) - a[i][j]) / (pow(4.0,m+1)-1);
      }
    }
  }
  
  // Update the array with the gradient
  for (unsigned short i=0; i<n_args; i+=1) grad[i] = a[0][i];
  free(a);
}


void genD(double(*f)(unsigned short n_args, double* args, void* box), 
              unsigned short n_args, double* args, double* D, void* box,
              double eps, double d, double zero_tol, unsigned short r, unsigned short v){
  
  if (v!=2){
    printf("Base R code assumed v=2");
    return;
  }
  
  double f0, f1, f2;
  // f0 is the value of the objective function at the given argument
  f0 = f(n_args,args,box);
  
  double* h0 = (double *) malloc(sizeof(double) * n_args);
  for (unsigned short i=0; i<n_args; i+=1) {
    h0[i] = fabs(d*args[i]) + eps * ((fabs(args[i])<zero_tol) ? 1 : 0);
  }
  double* h = (double *) malloc(sizeof(double) * n_args);
  
  double* Daprox = (double *) malloc(sizeof(double) * r);
  for (unsigned short i=0; i<r; i+=1) Daprox[i] = 0;
  
  double* Hdiag = (double *) malloc(sizeof(double) * n_args);
  for (unsigned short i=0; i<n_args; i+=1) Hdiag[i] = 0;
  
  double* Haprox = (double *) malloc(sizeof(double) * r);
  for (unsigned short i=0; i<r; i+=1) Haprox[i] = 0;  

  double* move_plus = (double *) malloc(sizeof(double) * n_args);
  double* move_minus = (double *) malloc(sizeof(double) * n_args);
  
  for (unsigned short i=0; i<n_args; i+=1) { //each parameter  - first deriv. & hessian diagonal
    
    for (unsigned short xx=0; xx<n_args; xx+=1) h[xx] = h0[xx];
    
    for (unsigned short k=0; k<r; k+=1){ // successively reduce h 
      
      for (unsigned short xx=0; xx<n_args; xx+=1){
        move_plus[xx] = args[xx];
        move_minus[xx] = args[xx];
      }
      move_plus[i] += h[i];
      move_minus[i] -= h[i];
      
      f1 = (*f)(n_args, move_plus, box);
      f2 = (*f)(n_args, move_minus, box);
      
      Daprox[k] = (f1 - f2)  / (2*h[i]); // F'(i)
      Haprox[k] = (f1-2*f0+f2) / pow(h[i],2); // F''(i,i) hessian diagonal

      for (unsigned short i=0; i<n_args; i+=1) {
        h[i] = h[i]/v; // Reduced h by 1/v.
      }
    }

    for (unsigned short m=0; m<r-1; m+=1) {
      for (unsigned short k=0; k<r-m-1; k+=1) {
        Daprox[k] = (Daprox[k+1] * pow(4,m+1.0) - Daprox[k]) / (pow(4,m+1.0)-1);
        Haprox[k] = (Haprox[k+1] * pow(4,m+1.0) - Haprox[k]) / (pow(4,m+1.0)-1);
      }
    }
    D[i] = Daprox[0];
    Hdiag[i] = Haprox[0];
  }

  unsigned short u = n_args-1;
  
  for (unsigned short i=0; i<n_args; i+=1) { // 2nd derivative  - do lower half of hessian only
    for (unsigned short j=0; j<=i; j+=1) {
      
      u = u + 1;
      if (i==j){
        D[u] = Hdiag[i];
      }
      else {
        
        for (unsigned short xx=0; xx<n_args; xx+=1) h[xx] = h0[xx];
        
        for (unsigned short k=0; k<r; k+=1){ // successively reduce h 
          
          for (unsigned short xx=0; xx<n_args; xx+=1){
            move_plus[xx] = args[xx];
            move_minus[xx] = args[xx];
          }
          move_plus[i] += h[i];
          move_minus[i] -= h[i];
          move_plus[j] += h[j];
          move_minus[j] -= h[j];
          
          f1 = (*f)(n_args, move_plus, box);
          f2 = (*f)(n_args, move_minus, box);
          
          // F''(i,j)          
          Daprox[k] = (f1 - 2*f0 + f2 -
                      Hdiag[i]*pow(h[i],2) -
                      Hdiag[j]*pow(h[j],2)) / (2*h[i]*h[j]);
          
          for (unsigned short xx=0; xx<n_args; xx+=1) {
            h[xx] = h[xx]/v; // Reduced h by 1/v.
          }
        }

        for (unsigned short m=0; m<r-1; m+=1) {
          for (unsigned short k=0; k<r-m-1; k+=1) {
            Daprox[k] = (Daprox[k+1] * pow(4,m+1) - Daprox[k]) / (pow(4,m+1)-1);
          }
        }
        D[u] = Daprox[0];
      }
    }
  }

  free(h0);
  free(h);
  free(move_plus);
  free(move_minus);
  free(Daprox);
  free(Hdiag);
  free(Haprox);
}


void hessian(double(*f)(unsigned short n_args, double* args, void* box), 
          unsigned short n_args, double* args, double* H, void* box,
          double eps, double d, double zero_tol, unsigned short r, unsigned short v){
  
  // call the function to generate matrix of function derivatives information
  int N = n_args*(n_args+3)/2;
  double* D = (double *) malloc(sizeof(double) * N );
  for (unsigned short i=0; i<N; i+=1) D[i] = 0;
  genD(f, n_args, args, D, box, eps, d, zero_tol, r, v); 
  
  // reshape the 1D array into a Hessian matrix
  int u = n_args-1;
  for (unsigned short i=0; i<n_args; i+=1){
    for (unsigned short j=0; j<=i; j+=1){
      u = u+1;
      H[i*n_args+j] = D[u];
    }
  }
  for (unsigned short i=0; i<n_args; i+=1){
    for (unsigned short j=0; j<i; j+=1){
      H[j*n_args+i] = H[i*n_args+j];
    }
  }
  
  free(D);
}
