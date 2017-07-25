#ifndef model_hpp
#define model_hpp

#include <stdio.h>
#include <vector>

//double g(double x);
//double demand(int j, double *state, double *actions);
//double objective(double *state, double *actions);
//double D_objective_d_p0(double *state, double *actions);
//double objective_nlopt(const std::vector<double> &action, std::vector<double> &grad, void *data);
//double g_qual(double x);
double objective(const std::vector<double> &actions, std::vector<double> &grad, void *data);
// double D_objective_d_p0(double *state, const std::vector<double> &actions);
// double D_objective_d_x0(double *state, const std::vector<double> &actions);


#endif
