//
// numDeriv R library re-written into C
//
// Maciej Bak
// Swiss Institute of Bioinformatics
// 11.09.2019
//

// include the library header
#include "numDeriv.h"

// includes from standard C lubrary required for this file
#include <math.h>
#include <stdio.h>

// define a constant 
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

//=======================================================================================

//=======================================================================================

// Gradient objective function:
// F(a,b,c) = a^1 + b^2 + c^3
double x1x2x3(unsigned short n_args, double* args, void* box){
  return pow(*args,1) + pow(*(args+1),2) + pow(*(args+2),3);
}

// Gradient objective function:
// F(x) = arcsin(x)
double arcsin(unsigned short n_args, double* args, void* box){
  return asin(*args);
}

// Gradient objective function:
// F(x) = sin(1/x)
double sin1overX(unsigned short n_args, double* args, void* box){
  return sin(1 / *args);
}

// Gradient objective function:
// F(x) = (x-100)^2 + 1.e-06 * (x-300)^3
double polynomial(unsigned short n_args, double* args, void* box){
  return pow(*args-100,2) + 1.e-06 * pow(*args-300,3);
}

// Gradient objective function:
// F(x) = SUM_i[ sin(x_i) ]
double sum_sin(unsigned short n_args, double* args, void* box){
  double result = 0.0;
  for (unsigned short i=0; i<n_args; i+=1) result += sin( *(args+i) );
  return result;
}

// A very simple data strucutre to hold the parameters required
// to calculate values of the objective function
struct BOX_gradient{
  double* A;
};

// Gradient objective function:
// F(x) = SUM_i[A_i*x_i]
double implementation_test(unsigned short n_args, double* args, void* box){
  double result = 0.0;
  // de-reference the void pointer to a pointer to BOX_gradient instance
  struct BOX_gradient *box_with_parameters = (struct BOX_gradient*)box;
  double* parameters = box_with_parameters->A;
  // calculate the value of the function:
  for (unsigned short i=0; i<n_args; i+=1) result += parameters[i] * args[i];
  return result;
}

// A very simple data strucutre to hold the parameters required
// to calculate values of the objective function
struct BOX_genD{
  double X;
  double Y;
};

// genD objective function:
// F(th1,th2) = (th1*X / (th2+X)) - Y
double puromycin_test(unsigned short n_args, double* args, void* box){
  // de-reference the void pointer to a pointer to BOX_genD instance
  struct BOX_genD *box_with_parameters = (struct BOX_genD*)box;
  double X = box_with_parameters->X;
  double Y = box_with_parameters->Y;
  // calculate the value of the function:
  return (args[0]*X / (args[1]+X)) - Y;
}

// Hessian objective function:
// F(x) = sin(x)
double sin_test(unsigned short n_args, double* args, void* box){
  return sin(*args);
}

// Hessian objective function:
// F(x) = SUM_i[exp(2*x_i)]
double exp_test(unsigned short n_args, double* args, void* box){
  double result = 0.0;
  for (unsigned short i=0; i<n_args; i+=1) result += exp(2*args[i]);
  return result;
}

// A very simple data strucutre to hold the parameters required
// to calculate values of the objective function
struct BOX_hessian{
  double A;
};

// hessian objective function:
// F(x,y) = A*x*y;
double hessian_test(unsigned short n_args, double* args, void* box){
  // de-reference the void pointer to a pointer to BOX_hessian instance
  struct BOX_hessian *box_with_parameters = (struct BOX_hessian*)box;
  double A = box_with_parameters->A;
  // calculate the value of the function:
  return A*args[0]*args[1];
}

//=======================================================================================

int main(int argc, const char **argv){
  
    
  // ===== Simple gradient test =====
  // F(a,b,c) = a^1 + b^2 + c^3
  double args_test1[] = {10.0,20.0,30.0};
  double grad_test1[] = {0.0,0.0,0.0};
  gradient(x1x2x3, 3, args_test1, grad_test1, 0, 1e-4, 1e-4, 2.220446e-16/7e-7, 4, 2); 
  printf("%f %f %f\n",*grad_test1,*(grad_test1+1),*(grad_test1+2));
  //1 40 2700
  
  
  // ===== Arcsin gradient test =====
  // F(x) = arcsin(x)
  double args_test2a[] = {0.9};
  double grad_test2a[] = {0.0}; 
  gradient(arcsin, 1, args_test2a, grad_test2a, 0, 1e-4, 1e-4, 2.220446e-16/7e-7, 4, 2); 
  printf("%f\n",*grad_test2a);
  // 2.294157
  double args_test2b[] = {0.99};
  double grad_test2b[] = {0.0}; 
  gradient(arcsin, 1, args_test2b, grad_test2b, 0, 1e-4, 1e-4, 2.220446e-16/7e-7, 4, 2); 
  printf("%f\n",*grad_test2b);
  // 7.088812
  double args_test2c[] = {0.999};
  double grad_test2c[] = {0.0}; 
  gradient(arcsin, 1, args_test2c, grad_test2c, 0, 1e-4, 1e-4, 2.220446e-16/7e-7, 4, 2); 
  printf("%f\n",*grad_test2c);
  // 22.366272
  
  
  // ===== Sin(1/x) gradient test =====
  // F(x) = sin(1/x)
  double args_test3a[] = {0.1};
  double grad_test3a[] = {0.0}; 
  gradient(sin1overX, 1, args_test3a, grad_test3a, 0, 1e-4, 1e-4, 2.220446e-16/7e-7, 4, 2); 
  printf("%f\n",*grad_test3a);
  // 8.390715e+01
  double args_test3b[] = {0.01};
  double grad_test3b[] = {0.0}; 
  gradient(sin1overX, 1, args_test3b, grad_test3b, 0, 1e-4, 1e-4, 2.220446e-16/7e-7, 4, 2); 
  printf("%f\n",*grad_test3b);
  // -8.623189e+03
  double args_test3c[] = {0.001};
  double grad_test3c[] = {0.0}; 
  gradient(sin1overX, 1, args_test3c, grad_test3c, 0, 1e-4, 1e-4, 2.220446e-16/7e-7, 4, 2); 
  printf("%f\n",*grad_test3c);
  // -5.623791e+05
  double args_test3d[] = {0.0001};
  double grad_test3d[] = {0.0}; 
  gradient(sin1overX, 1, args_test3d, grad_test3d, 0, 1e-4, 1e-4, 2.220446e-16/7e-7, 4, 2); 
  printf("%f\n",*grad_test3d);
  // 9.521554e+07
  
  
  // ===== Polynomial gradient test =====
  // F(x) = (x-100)^2 + 1.e-06 * (x-300)^3
  double args_test4a[] = {100.001};
  double grad_test4a[] = {0.0}; 
  gradient(polynomial, 1, args_test4a, grad_test4a, 0, 1e-4, 1e-4, 2.220446e-16/7e-7, 4, 2); 
  printf("%f\n",*grad_test4a);
  // 0.1219988
  double args_test4b[] = {300.001};
  double grad_test4b[] = {0.0}; 
  gradient(polynomial, 1, args_test4b, grad_test4b, 0, 1e-4, 1e-4, 2.220446e-16/7e-7, 4, 2); 
  printf("%f\n",*grad_test4b);
  // 400.0020000
  
  
  // ===== SUM[sin(x)] gradient test =====
  // F(x) = SUM_i[ sin(x_i) ]
  double args_test5[11];
  for (unsigned short i=0; i<11; i+=1) args_test5[i] = i*2.0*M_PI/10.0;
  double grad_test5[] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}; 
  gradient(sum_sin, 11, args_test5, grad_test5, 0, 1e-4, 1e-4, 2.220446e-16/7e-7, 4, 2); 
  printf("%f %f %f %f\n",*grad_test5,*(grad_test5+1),*(grad_test5+2),*(grad_test5+3));
  printf("%f %f %f %f\n",*(grad_test5+4),*(grad_test5+5),*(grad_test5+6),*(grad_test5+7));
  printf("%f %f %f\n",*(grad_test5+8),*(grad_test5+9),*(grad_test5+10));
  // 1.000000 0.809017 0.309017 -0.309017
  // -0.809017 -1.000000 -0.809017 -0.309017
  // 0.309017 0.809017 1.000000

  
  // ===== Simple gradient implementation test =====
  // F(x) = SUM_i[A_i*x_i]
  double args_test6[] = {1.0,1.0,1.0,1.0};
  double params_test6[] = {10.0,-5.0,0.99,0.5};
  struct BOX_gradient gradient_box_for_parameters;
  gradient_box_for_parameters.A = params_test6;
  void* gradient_void_pointer_box = &gradient_box_for_parameters;
  double grad_test6[] = {0.0,0.0,0.0,0.0};
  gradient(implementation_test, 4, args_test6, grad_test6, gradient_void_pointer_box, 1e-4, 1e-4, 2.220446e-16/7e-7, 4, 2); 
  printf("%f %f %f %f\n",*grad_test6,*(grad_test6+1),*(grad_test6+2),*(grad_test6+3));
  // gradient is given by params:
  // 10.000000 -5.000000 0.990000 0.500000


  // ===== Simple genD test =====
  // F(a,b,c) = a^1 + b^2 + c^3
  double args_test7[] = {10.0,20.0,30.0};
  double D_test7[] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  genD(x1x2x3, 3, args_test7, D_test7, 0, 1e-4, 1e-4, 2.220446e-16/7e-7, 4, 2); 
  printf("%.10f %.10f %.10f\n",*D_test7,*(D_test7+1),*(D_test7+2));
  printf("%.10f %.10f %.10f\n",*(D_test7+3),*(D_test7+4),*(D_test7+5));
  printf("%.10f %.10f %.10f\n",*(D_test7+6),*(D_test7+7),*(D_test7+8));
  // 0.9999999826 40.0000000026 2699.9999999738
  // 0.0000000000 0.0000840983 1.9999376037
  // 0.0000000000 0.0000000000 179.9999994949
  

  // ===== R:numDeriv puromycin genD test =====
  // F(th1,th2) = (th1*X / (th2+X)) - Y
  double args_test8[] = {212.7000,0.0641};
  double D_test8[] = {0.0,0.0,0.0,0.0,0.0};
  double X_test8[] = {0.02,0.02,0.06,0.06,0.11,0.11,0.22,0.22,0.56,0.56,1.10,1.10};
  double Y_test8[] = {76,47,97,107,123,139,159,152,191,201,207,200};
  double numerical_D_test8[5][12];
  struct BOX_genD genD_box_for_parameters;
  void* genD_void_pointer_box = &genD_box_for_parameters;
  
  // Analytical solution
  double analytical_D_test8[5][12];
  for (unsigned short i=0; i<12; i+=1){
    analytical_D_test8[0][i] = X_test8[i] / (args_test8[1]+X_test8[i]);
    analytical_D_test8[1][i] = -args_test8[0]*X_test8[i] / pow((args_test8[1]+X_test8[i]),2);
    analytical_D_test8[2][i] = 0.0;
    analytical_D_test8[3][i] = -X_test8[i] / pow((args_test8[1]+X_test8[i]),2);
    analytical_D_test8[4][i] = 2*args_test8[0]*X_test8[i] / pow(args_test8[1]+X_test8[i],3);
    printf("%.15f %.15f %.15f %.15f %.15f\n",
      analytical_D_test8[0][i],analytical_D_test8[1][i],analytical_D_test8[2][i],analytical_D_test8[3][i],analytical_D_test8[4][i]);
  }
  printf("\n");

  // genD standard approximation
  for (unsigned short i=0; i<12; i+=1){
    for (unsigned short xx=0; xx<5; xx+=1) D_test8[xx] = 0.0; // reset D array
    // update the parameters
    genD_box_for_parameters.X = X_test8[i];
    genD_box_for_parameters.Y = Y_test8[i];
    // calculate the der. information array
    genD(puromycin_test, 2, args_test8, D_test8, genD_void_pointer_box, 1e-4, 1e-4, 2.220446e-16/7e-7, 4, 2);
    printf("%.15f %.15f %.15f %.15f %.15f\n",*D_test8,*(D_test8+1),*(D_test8+2),*(D_test8+3),*(D_test8+4));
    numerical_D_test8[0][i] = D_test8[0];
    numerical_D_test8[1][i] = D_test8[1];
    numerical_D_test8[2][i] = D_test8[2];
    numerical_D_test8[3][i] = D_test8[3];
    numerical_D_test8[4][i] = D_test8[4];
  }
  printf("\n");
  // print the difference between the approximation
  // and the analytical solution
  for (unsigned short i=0; i<12; i+=1){
    printf("%.15f %.15f %.15f %.15f %.15f\n",
        analytical_D_test8[0][i] - numerical_D_test8[0][i],
        analytical_D_test8[1][i] - numerical_D_test8[1][i],
        analytical_D_test8[2][i] - numerical_D_test8[2][i],
        analytical_D_test8[3][i] - numerical_D_test8[3][i],
        analytical_D_test8[4][i] - numerical_D_test8[4][i]);
  }
  printf("\n");
  
  // genD approximation with d=0.01 - closer approximation!
  for (unsigned short i=0; i<12; i+=1){
    for (unsigned short xx=0; xx<5; xx+=1) D_test8[xx] = 0.0; // reset D array
    // update the parameters
    genD_box_for_parameters.X = X_test8[i];
    genD_box_for_parameters.Y = Y_test8[i];
    // calculate the der. information array
    genD(puromycin_test, 2, args_test8, D_test8, genD_void_pointer_box, 1e-4, 1e-2, 2.220446e-16/7e-7, 4, 2);
    printf("%.15f %.15f %.15f %.15f %.15f\n",*D_test8,*(D_test8+1),*(D_test8+2),*(D_test8+3),*(D_test8+4));
    numerical_D_test8[0][i] = D_test8[0];
    numerical_D_test8[1][i] = D_test8[1];
    numerical_D_test8[2][i] = D_test8[2];
    numerical_D_test8[3][i] = D_test8[3];
    numerical_D_test8[4][i] = D_test8[4];
  }
  printf("\n");
  // print the difference between the approximation
  // and the analytical solution
  for (unsigned short i=0; i<12; i+=1){
    printf("%.15f %.15f %.15f %.15f %.15f\n",
        analytical_D_test8[0][i] - numerical_D_test8[0][i],
        analytical_D_test8[1][i] - numerical_D_test8[1][i],
        analytical_D_test8[2][i] - numerical_D_test8[2][i],
        analytical_D_test8[3][i] - numerical_D_test8[3][i],
        analytical_D_test8[4][i] - numerical_D_test8[4][i]);
  }
  printf("\n");
  
  // genD approximation with d=0.01 and r=6 - further approximation!
  for (unsigned short i=0; i<12; i+=1){
    for (unsigned short xx=0; xx<5; xx+=1) D_test8[xx] = 0.0; // reset D array
    // update the parameters
    genD_box_for_parameters.X = X_test8[i];
    genD_box_for_parameters.Y = Y_test8[i];
    // calculate the der. information array
    genD(puromycin_test, 2, args_test8, D_test8, genD_void_pointer_box, 1e-4, 1e-2, 2.220446e-16/7e-7, 6, 2);
    printf("%.15f %.15f %.15f %.15f %.15f\n",*D_test8,*(D_test8+1),*(D_test8+2),*(D_test8+3),*(D_test8+4));
    numerical_D_test8[0][i] = D_test8[0];
    numerical_D_test8[1][i] = D_test8[1];
    numerical_D_test8[2][i] = D_test8[2];
    numerical_D_test8[3][i] = D_test8[3];
    numerical_D_test8[4][i] = D_test8[4];
  }
  printf("\n");
  // print the difference between the approximation
  // and the analytical solution
  for (unsigned short i=0; i<12; i+=1){
    printf("%.15f %.15f %.15f %.15f %.15f\n",
        analytical_D_test8[0][i] - numerical_D_test8[0][i],
        analytical_D_test8[1][i] - numerical_D_test8[1][i],
        analytical_D_test8[2][i] - numerical_D_test8[2][i],
        analytical_D_test8[3][i] - numerical_D_test8[3][i],
        analytical_D_test8[4][i] - numerical_D_test8[4][i]);
  }
  printf("\n");

  
  // ===== Simple hessian test =====
  // F(a,b,c) = a^1 + b^2 + c^3
  double args_test9[] = {10.0,20.0,30.0};
  double H_test9[3][3];
  hessian(x1x2x3, 3, args_test9, (double *)H_test9, 0, 1e-4, 1e-4, 2.220446e-16/7e-7, 4, 2); 
  printf("%.2f %.2f %.2f\n",H_test9[0][0],H_test9[0][1],H_test9[0][2]);
  printf("%.2f %.2f %.2f\n",H_test9[1][0],H_test9[1][1],H_test9[1][2]);
  printf("%.2f %.2f %.2f\n",H_test9[2][0],H_test9[2][1],H_test9[2][2]);
  // 0.00 0.00 0.00
  // 0.00 2.00 0.00
  // 0.00 0.00 180.00
  
  
  // ===== Sin hessian test =====
  // F(x) = sin(x)
  double args_test10[] = {0.25*M_PI};
  double H_test10[1][1];
  // Analytical solution
  printf("%.10f\n", sin(*args_test10+M_PI));
  // -0.7071067812
  // Numerical approximation: eps=1e-8 and d=1e-2 results in a close approximation
  hessian(sin_test, 1, args_test10, (double *)H_test10, 0, 1e-8, 1e-2, 2.220446e-16/7e-7, 4, 2); 
  printf("%.10f\n",H_test10[0][0]);
  // -0.7071067811

  
  // ===== Exp hessian test =====
  // F(x) = SUM_i[exp(2*x_i)]
  double args_test11[] = {1,3,5};
  double analytical_H[3][3] = {
    {4*exp(2*args_test11[0]),0,0}, {0,4*exp(2*args_test11[1]),0}, {0,0,4*exp(2*args_test11[2])}
  };

  // Analytical solution
  printf("%.7f %.7f %.7f\n",analytical_H[0][0],analytical_H[0][1],analytical_H[0][2]);
  printf("%.7f %.7f %.7f\n",analytical_H[1][0],analytical_H[1][1],analytical_H[1][2]);
  printf("%.7f %.7f %.7f\n",analytical_H[2][0],analytical_H[2][1],analytical_H[2][2]);
  printf("\n");
  // 29.5562244 0.0000000 0.0000000
  // 0.0000000 1613.7151740 0.0000000
  // 0.0000000 0.0000000 88105.8631792
  
  // Numerical approximation: eps=1e-4 and d=1e-2 results in a close approximation
  double H_test11[3][3];
  hessian(exp_test, 3, args_test11, (double *)H_test11, 0, 1e-4, 1e-2, 2.220446e-16/7e-7, 4, 2); 
  printf("%.7f %.7f %.7f\n",H_test11[0][0],H_test11[0][1],H_test11[0][2]);
  printf("%.7f %.7f %.7f\n",H_test11[1][0],H_test11[1][1],H_test11[1][2]);
  printf("%.7f %.7f %.7f\n",H_test11[2][0],H_test11[2][1],H_test11[2][2]);
  printf("\n");
  // 29.5562244 0.0000000 0.0000000
  // 0.0000000 1613.7151740 0.0000000
  // 0.0000000 0.0000000 88105.8631792  
  
  // print the difference between the approximation
  // and the analytical solution
  for (unsigned short i=0; i<3; i+=1){
    printf("%.15f %.15f %.15f\n",
        analytical_H[0][i] - H_test11[0][i],
        analytical_H[1][i] - H_test11[1][i],
        analytical_H[2][i] - H_test11[2][i]);
  }
  printf("\n");

  
  // ===== Simple hessian implementation test =====
  // F(x,y) = A*x*y;
  double args_test12[] = {1.0,1.0};
  struct BOX_hessian hessian_box_for_parameters;
  hessian_box_for_parameters.A = -100.0;
  void* hessian_void_pointer_box = &hessian_box_for_parameters;
  double H_test12[2][2];
  hessian(hessian_test, 2, args_test12, (double *)H_test12, hessian_void_pointer_box, 1e-4, 1e-4, 2.220446e-16/7e-7, 4, 2); 
  printf("%.15f %.15f\n",H_test12[0][0],H_test12[0][1]);
  printf("%.15f %.15f\n",H_test12[1][0],H_test12[1][1]);
  // -0.000010947621410 -100.000049023600340
  // -100.000049023600340 -0.000010947621410 
  

  return 0;
}

