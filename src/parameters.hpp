#ifndef parameters_hpp
#define parameters_hpp

#include <string>
#include <stdio.h>
#include <dlib/dnn.h>
#include <dlib/data_io.h>
using namespace std;
using namespace dlib;

extern size_t d;
extern size_t Nmax;
extern size_t num_of_knots;
extern size_t S;
extern size_t num_of_test;
extern size_t S_test;
extern size_t num_of_coef;
extern size_t num_of_alpha;
extern size_t coef_degree;
extern double M;
extern double delta;
extern double alpha;
extern double v;
extern double beta;
extern double c;
extern double eta;
extern double theta;
extern double *xmin;
extern double *xmax;
extern double *value_coef;
extern double *policy_coef;
extern double *price_coef;
extern int **alpha_coef;
extern int approx_type;
extern int func_count;

typedef struct {
    double *state;
    double p_other;
    double policy_other;
} objective_data;

void initalize();
#endif /* parameters_hpp */
