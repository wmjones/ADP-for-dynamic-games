#ifndef fitting_hpp
#define fitting_hpp

#include <stdio.h>
#include <math.h>
#include "parameters.hpp"

double predict(double *state);
void fitting(double **xy_data, double *v_hat, double *dz_data0, double *dz_data1, double *b, double **z_knots);

inline double cheb_T(int n, double x){
    if(n==-1){
        return 0.0;
    }
    else{
        return cos(n*acos(x));
    }
}

inline double cheb_U(int n, double x){
    if(n==-1){
        return 0.0;
    }
    else{
        double a = acos(x);
        return sin((n+1)*a)/sin(a);
    }
}

inline double T_alpha(double *z, int *alpha){
    double product = 1;
    for(size_t i=0; i<d; i++){
        product *= cheb_T(alpha[i], z[i]);
    }
    return product;
}

inline double dT_alpha(double *z, int del, int *alpha){
    if(del==0){
        return alpha[0]*cheb_U(alpha[0]-1, z[0])*cheb_T(alpha[1], z[1])*2/(xmax[0]-xmin[0]);
    }
    else{
        return cheb_T(alpha[0], z[0])*alpha[1]*cheb_U(alpha[1]-1, z[1])*2/(xmax[1]-xmin[1]);
    }
}

inline double V_hat(double *x, double *coef, int **alpha, size_t num_of_alpha){
    double sum = 0;
    double z[d];
    for(size_t j=0; j<d; j++){
        z[j] = (2*x[j] - xmin[j] - xmax[j])/(xmax[j] - xmin[j]);
    }
    for(size_t i=0; i<num_of_alpha; i++){
        sum += coef[i]*T_alpha(z, alpha[i]);
    }
    return sum;
}

inline double dV_hat(double *x, double *b, int del, int **alpha, size_t num_of_alpha){
    double sum = 0;
    double z[d];
    for(size_t j=0; j<d; j++){
        z[j] = (2*x[j] - xmin[j] - xmax[j])/(xmax[j] - xmin[j]);
    }
    for(size_t i=0; i<num_of_alpha; i++){
        sum += b[i]*dT_alpha(z, del, alpha[i]);
    }
    return sum;
}

#endif
