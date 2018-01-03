#ifndef model_hpp
#define model_hpp

#include <stdio.h>
#include <vector>

double objective(const std::vector<double> &actions, std::vector<double> &grad, void *data);

#endif
