#include "parameters.hpp"
#include <math.h>
#include <stdio.h>

size_t d;
size_t num_of_actions;
size_t Nmax;
size_t num_of_knots;
size_t S;
size_t num_of_test;
size_t S_test;
size_t num_of_coef;
size_t coef_degree;
double M;
double delta;
double alpha;
double v;
double beta;
double c;
double eta;
double theta;
double *xmin;
double *xmax;
double *value_coef;
double *policy_coef;
double *price_coef;
int **alpha_coef;
std::string approx_type = "cai";
std::string cai("cai");
std::string cheb("cheb");
std::string ann("ann");


const size_t choose(size_t n, size_t k){
    if(k==0){
        return 1;
    }
    else{
        return n*choose(n-1, k-1)/k;
    }
}

__attribute__((constructor)) void initialize(){
    d = 2;
    Nmax = 100;
    num_of_knots = 10;
    S = num_of_knots*num_of_knots;
    num_of_test = 150;
    S_test = num_of_test*num_of_test;
    coef_degree=0;
    M = 5.0;
    beta = 0.925;
    delta = 0.7;
    alpha = 3.0;
    c = 5.0;
    eta = 1;
    v = 1;
    theta = 1.5;
    xmin = new double[d];
    xmin[0] = 1.0; xmin[1] = 1.0;
    xmax = new double[d];
    xmax[0] = 18.0; xmax[1] = 18.0;

    if(approx_type.compare(cai)==0 || approx_type.compare(cheb)==0){
	while(choose(coef_degree+1 + d, d)<pow(num_of_knots, d)){
	    coef_degree += 1;
	}
	coef_degree = coef_degree/5*2;
	coef_degree = 12;
	printf("coef_degree = %zd\n", coef_degree);
	num_of_coef = choose(coef_degree + d, d);
	printf("num_of_coef = %zd\n", num_of_coef);
	value_coef = new double[num_of_coef];
	policy_coef = new double[num_of_coef];
	price_coef = new double[num_of_coef];
	alpha_coef = new int*[num_of_coef];
	for(size_t i=0; i<num_of_coef; i++){
	    alpha_coef[i] = new int[d];
	    policy_coef[i] = 0.0;
	    price_coef[i] = 0.0;
	    value_coef[i] = 0.0;
	}
	size_t k=0;
	for(int i=0; i<=(int)coef_degree; i++){
	    for(int j=0; j<=(int)coef_degree; j++){
		if(i+j<=(int)coef_degree){
		    alpha_coef[k][0] = i;
		    alpha_coef[k][1] = j;
		    k+=1;
		}
	    }
	}
    }
    else{
	printf("Value Function Approximation with ANN\n");
    }

    // value_coef_degree = 2*num_of_knots - 1;
    // value_coef_degree = policy_coef_degree;
    // printf("value_coef_degree = %zd\n", value_coef_degree);
    // num_of_value_coef = choose(value_coef_degree+d, d);
    // value_coef = new double[num_of_value_coef];
    // alpha_value_coef = new int*[num_of_value_coef];
    // for(size_t i=0; i<num_of_value_coef; i++){
    //     alpha_value_coef[i] = new int[d];
    //     value_coef[i] = 0.0;
    // }
    // k=0;
    // for(int i=0; i<=(int)value_coef_degree; i++){
    //     for(int j=0; j<=(int)value_coef_degree; j++){
    //         if(i+j<=(int)value_coef_degree){
    //             alpha_value_coef[k][0] = i;
    //             alpha_value_coef[k][1] = j;
    //             k+=1;
    //         }
    //     }
    // }
}
