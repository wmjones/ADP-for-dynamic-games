#include "model.hpp"
#include "fitting.hpp"
#include "parameters.hpp"
#include <math.h>
#include <vector>


// inline double g(double x){
//     return 12.6/(1 + exp(4-x));
// }

// inline double D_g(double x){
//     double a1 = exp(4-x);
//     return 12.6*a1/((1+a1)*(1+a1));
// }

inline double g(double x){
    return 14.7/(1 + exp(4-1.2*x))-2;
}

// inline double D_g(double x){
//     double a1 = exp(4-x);
//     return 14.7*a1/((1+a1)*(1+a1));
// }


inline double profit(double *state, double p_other, const std::vector<double> &actions){
    return M*exp(g(state[0]) - actions[0])/(1 + exp(g(state[0])-actions[0]) + exp(g(state[1]) - p_other))*(actions[0] - c);
}

inline double expectation(double *state, double policy_other, const std::vector<double> &actions){
    double p[9];
    double x0 = actions[1];
    double x1 = policy_other;
    double a1 = delta/(1 + alpha*x0);
    double b1 = delta/(1 + alpha*x1);
    double a2 = (1 - delta + delta*alpha*x0)/(1 + alpha*x0);
    double b2 = (1 - delta + delta*alpha*x1)/(1 + alpha*x1);
    double a3 = (alpha*x0 - delta*alpha*x0)/(1 + alpha*x0);
    double b3 = (alpha*x1 - delta*alpha*x1)/(1 + alpha*x1);
    double a4 = (1 + delta*alpha*x0)/(1 + alpha*x0);
    double b4 = (1 + delta*alpha*x1)/(1 + alpha*x1);
    double a5 = (1 - delta + alpha*x0)/(1 + alpha*x0);
    double b5 = (1 - delta + alpha*x1)/(1 + alpha*x1);
    double sprime0[2] = {state[0]-1, state[1]-1};
    double sprime1[2] = {state[0], state[1]-1};
    double sprime2[2] = {state[0]+1, state[1]-1};
    double sprime3[2] = {state[0]-1, state[1]};
    double sprime4[2] = {state[0], state[1]};
    double sprime5[2] = {state[0]+1, state[1]};
    double sprime6[2] = {state[0]-1, state[1]+1};
    double sprime7[2] = {state[0], state[1]+1};
    double sprime8[2] = {state[0]+1, state[1]+1};
    if(state[0]-1<xmin[0] && state[1]-1<xmin[1]){
        double v4 = predict(sprime4);
        double v5 = predict(sprime5);
        double v7 = predict(sprime7);
        double v8 = predict(sprime8);
        p[0] = 0;
        p[1] = 0;
        p[2] = 0;
        p[3] = 0;
        p[4] = a4*b4*v4;
        p[5] = a3*b4*v5;
        p[6] = 0;
        p[7] = a4*b3*v7;
        p[8] = a3*b3*v8;
    }
    else if(state[0]-1<xmin[0] && state[1]-1>=xmin[1] && state[1]+1<=xmax[1]){
        double v1 = predict(sprime1);
        double v2 = predict(sprime2);
        double v4 = predict(sprime4);
        double v5 = predict(sprime5);
        double v7 = predict(sprime7);
        double v8 = predict(sprime8);
        p[0] = 0;
        p[1] = a4*b1*v1;
        p[2] = a3*b1*v2;
        p[3] = 0;
        p[4] = a4*b2*v4;
        p[5] = a3*b2*v5;
        p[6] = 0;
        p[7] = a4*b3*v7;
        p[8] = a3*b3*v8;
    }
    else if(state[0]-1<xmin[0] && state[1]+1>xmax[1]){
        double v1 = predict(sprime1);
        double v2 = predict(sprime2);
        double v4 = predict(sprime4);
        double v5 = predict(sprime5);
        p[0] = 0;
        p[1] = a4*b1*v1;
        p[2] = a3*b1*v2;
        p[3] = 0;
        p[4] = a4*b5*v4;
        p[5] = a3*b5*v5;
        p[6] = 0;
        p[7] = 0;
        p[8] = 0;
    }
    else if(state[0]-1>=xmin[0] && state[0]+1<=xmax[0] && state[1]-1<xmin[1]){
        double v3 = predict(sprime3);
        double v4 = predict(sprime4);
        double v5 = predict(sprime5);
        double v6 = predict(sprime6);
        double v7 = predict(sprime7);
        double v8 = predict(sprime8);
        p[0] = 0;
        p[1] = 0;
        p[2] = 0;
        p[3] = a1*b4*v3;
        p[4] = a2*b4*v4;
        p[5] = a3*b4*v5;
        p[6] = a1*b3*v6;
        p[7] = a2*b3*v7;
        p[8] = a3*b3*v8;
    }
    else if(state[0]-1>=xmin[0] && state[0]+1<=xmax[0] && state[1]-1>=xmin[1] && state[1]+1<=xmax[1]){
        double v0 = predict(sprime0);
        double v1 = predict(sprime1);
        double v2 = predict(sprime2);
        double v3 = predict(sprime3);
        double v4 = predict(sprime4);
        double v5 = predict(sprime5);
        double v6 = predict(sprime6);
        double v7 = predict(sprime7);
        double v8 = predict(sprime8);
        p[0] = a1*b1*v0;
        p[1] = a2*b1*v1;
        p[2] = a3*b1*v2;
        p[3] = a1*b2*v3;
        p[4] = a2*b2*v4;
        p[5] = a3*b2*v5;
        p[6] = a1*b3*v6;
        p[7] = a2*b3*v7;
        p[8] = a3*b3*v8;
    }
    else if(state[0]-1>=xmin[0] && state[0]+1<=xmax[0] && state[1]+1>xmax[1]){
        double v0 = predict(sprime0);
        double v1 = predict(sprime1);
        double v2 = predict(sprime2);
        double v3 = predict(sprime3);
        double v4 = predict(sprime4);
        double v5 = predict(sprime5);
        p[0] = a1*b1*v0;
        p[1] = a2*b1*v1;
        p[2] = a3*b1*v2;
        p[3] = a1*b5*v3;
        p[4] = a2*b5*v4;
        p[5] = a3*b5*v5;
        p[6] = 0;
        p[7] = 0;
        p[8] = 0;
    }
    else if(state[0]+1>xmax[0] && state[1]-1<xmin[1]){
        double v3 = predict(sprime3);
        double v4 = predict(sprime4);
        double v6 = predict(sprime6);
        double v7 = predict(sprime7);
        p[0] = 0;
        p[1] = 0;
        p[2] = 0;
        p[3] = a1*b4*v3;
        p[4] = a5*b4*v4;
        p[5] = 0;
        p[6] = a1*b3*v6;
        p[7] = a5*b3*v7;
        p[8] = 0;
    }
    else if(state[0]+1>xmax[0] && state[1]-1>=xmin[1] && state[1]+1<=xmax[1]){
        double v0 = predict(sprime0);
        double v1 = predict(sprime1);
        double v3 = predict(sprime3);
        double v4 = predict(sprime4);
        double v6 = predict(sprime6);
        double v7 = predict(sprime7);
        p[0] = a1*b1*v0;
        p[1] = a5*b1*v1;
        p[2] = 0;
        p[3] = a1*b2*v3;
        p[4] = a5*b2*v4;
        p[5] = 0;
        p[6] = a1*b3*v6;
        p[7] = a5*b3*v7;
        p[8] = 0;
    }
    else if(state[0]+1>xmax[0] && state[1]+1>xmax[1]){
        double v0 = predict(sprime0);
        double v1 = predict(sprime1);
        double v3 = predict(sprime3);
        double v4 = predict(sprime4);
        p[0] = a1*b1*v0;
        p[1] = a5*b1*v1;
        p[2] = 0;
        p[3] = a1*b5*v3;
        p[4] = a5*b5*v4;
        p[5] = 0;
        p[6] = 0;
        p[7] = 0;
        p[8] = 0;
    }
    else{
        printf("ERROR: expectation reaching region it shouldn't\n");
	printf("state = (%f, %f\n)", state[0], state[1]);
    }
    double out = 0;
    for(size_t i=0; i<9; i++){
        out += p[i];
    }
    return out;
}

inline double D_expectation_d_x0(double *state, double policy_other, const std::vector<double> &actions){
    double p[9];
    double x0 = actions[1];
    double x1 = policy_other;
    double c1 = ((1+alpha*x0)*(1+alpha*x0));
    double a1 = -alpha*delta/c1;
    double b1 = delta/(1 + alpha*x1);
    double a2 = (alpha*(-1+2*delta))/c1;
    double b2 = (1 - delta + delta*alpha*x1)/(1 + alpha*x1);
    double a3 = (alpha - alpha*delta)/c1;
    double b3 = (alpha*x1 - delta*alpha*x1)/(1 + alpha*x1);
    double a4 = (alpha*delta - alpha)/c1;
    double b4 = (1 + delta*alpha*x1)/(1 + alpha*x1);
    double a5 = alpha*delta/c1;
    double b5 = (1 - delta + alpha*x1)/(1 + alpha*x1);
    double sprime0[2] = {state[0]-1, state[1]-1};
    double sprime1[2] = {state[0], state[1]-1};
    double sprime2[2] = {state[0]+1, state[1]-1};
    double sprime3[2] = {state[0]-1, state[1]};
    double sprime4[2] = {state[0], state[1]};
    double sprime5[2] = {state[0]+1, state[1]};
    double sprime6[2] = {state[0]-1, state[1]+1};
    double sprime7[2] = {state[0], state[1]+1};
    double sprime8[2] = {state[0]+1, state[1]+1};
    if(state[0]-1<xmin[0] && state[1]-1<xmin[1]){
        double v4 = predict(sprime4);
        double v5 = predict(sprime5);
        double v7 = predict(sprime7);
        double v8 = predict(sprime8);
        p[0] = 0;
        p[1] = 0;
        p[2] = 0;
        p[3] = 0;
        p[4] = a4*b4*v4;
        p[5] = a3*b4*v5;
        p[6] = 0;
        p[7] = a4*b3*v7;
        p[8] = a3*b3*v8;
    }
    else if(state[0]-1<xmin[0] && state[1]-1>=xmin[1] && state[1]+1<=xmax[1]){
        double v1 = predict(sprime1);
        double v2 = predict(sprime2);
        double v4 = predict(sprime4);
        double v5 = predict(sprime5);
        double v7 = predict(sprime7);
        double v8 = predict(sprime8);
        p[0] = 0;
        p[1] = a4*b1*v1;
        p[2] = a3*b1*v2;
        p[3] = 0;
        p[4] = a4*b2*v4;
        p[5] = a3*b2*v5;
        p[6] = 0;
        p[7] = a4*b3*v7;
        p[8] = a3*b3*v8;
    }
    else if(state[0]-1<xmin[0] && state[1]+1>xmax[1]){
        double v1 = predict(sprime1);
        double v2 = predict(sprime2);
        double v4 = predict(sprime4);
        double v5 = predict(sprime5);
        p[0] = 0;
        p[1] = a4*b1*v1;
        p[2] = a3*b1*v2;
        p[3] = 0;
        p[4] = a4*b5*v4;
        p[5] = a3*b5*v5;
        p[6] = 0;
        p[7] = 0;
        p[8] = 0;
    }
    else if(state[0]-1>=xmin[0] && state[0]+1<=xmax[0] && state[1]-1<xmin[1]){
        double v3 = predict(sprime3);
        double v4 = predict(sprime4);
        double v5 = predict(sprime5);
        double v6 = predict(sprime6);
        double v7 = predict(sprime7);
        double v8 = predict(sprime8);
        p[0] = 0;
        p[1] = 0;
        p[2] = 0;
        p[3] = a1*b4*v3;
        p[4] = a2*b4*v4;
        p[5] = a3*b4*v5;
        p[6] = a1*b3*v6;
        p[7] = a2*b3*v7;
        p[8] = a3*b3*v8;
    }
    else if(state[0]-1>=xmin[0] && state[0]+1<=xmax[0] && state[1]-1>=xmin[1] && state[1]+1<=xmax[1]){
        double v0 = predict(sprime0);
        double v1 = predict(sprime1);
        double v2 = predict(sprime2);
        double v3 = predict(sprime3);
        double v4 = predict(sprime4);
        double v5 = predict(sprime5);
        double v6 = predict(sprime6);
        double v7 = predict(sprime7);
        double v8 = predict(sprime8);
        p[0] = a1*b1*v0;
        p[1] = a2*b1*v1;
        p[2] = a3*b1*v2;
        p[3] = a1*b2*v3;
        p[4] = a2*b2*v4;
        p[5] = a3*b2*v5;
        p[6] = a1*b3*v6;
        p[7] = a2*b3*v7;
        p[8] = a3*b3*v8;
    }
    else if(state[0]-1>=xmin[0] && state[0]+1<=xmax[0] && state[1]+1>xmax[1]){
        double v0 = predict(sprime0);
        double v1 = predict(sprime1);
        double v2 = predict(sprime2);
        double v3 = predict(sprime3);
        double v4 = predict(sprime4);
        double v5 = predict(sprime5);
        p[0] = a1*b1*v0;
        p[1] = a2*b1*v1;
        p[2] = a3*b1*v2;
        p[3] = a1*b5*v3;
        p[4] = a2*b5*v4;
        p[5] = a3*b5*v5;
        p[6] = 0;
        p[7] = 0;
        p[8] = 0;
    }
    else if(state[0]+1>xmax[0] && state[1]-1<xmin[1]){
        double v3 = predict(sprime3);
        double v4 = predict(sprime4);
        double v6 = predict(sprime6);
        double v7 = predict(sprime7);
        p[0] = 0;
        p[1] = 0;
        p[2] = 0;
        p[3] = a1*b4*v3;
        p[4] = a5*b4*v4;
        p[5] = 0;
        p[6] = a1*b3*v6;
        p[7] = a5*b3*v7;
        p[8] = 0;
    }
    else if(state[0]+1>xmax[0] && state[1]-1>=xmin[1] && state[1]+1<=xmax[1]){
        double v0 = predict(sprime0);
        double v1 = predict(sprime1);
        double v3 = predict(sprime3);
        double v4 = predict(sprime4);
        double v6 = predict(sprime6);
        double v7 = predict(sprime7);
        p[0] = a1*b1*v0;
        p[1] = a5*b1*v1;
        p[2] = 0;
        p[3] = a1*b2*v3;
        p[4] = a5*b2*v4;
        p[5] = 0;
        p[6] = a1*b3*v6;
        p[7] = a5*b3*v7;
        p[8] = 0;
    }
    else if(state[0]+1>xmax[0] && state[1]+1>xmax[1]){
        double v0 = predict(sprime0);
        double v1 = predict(sprime1);
        double v3 = predict(sprime3);
        double v4 = predict(sprime4);
        p[0] = a1*b1*v0;
        p[1] = a5*b1*v1;
        p[2] = 0;
        p[3] = a1*b5*v3;
        p[4] = a5*b5*v4;
        p[5] = 0;
        p[6] = 0;
        p[7] = 0;
        p[8] = 0;
    }
    else{
        printf("ERROR: derivative of expectation reaching region it shouldn't\n");
	printf("state = (%f, %f)\n", state[0], state[1]);
    }
    double out = 0;
    for(size_t i=0; i<9; i++){
        out += p[i];
    }
    // if (false){
    //     double v0 = predict(sprime0);
    //     double v1 = predict(sprime1);
    //     double v2 = predict(sprime2);
    //     double v3 = predict(sprime3);
    //     double v4 = predict(sprime4);
    //     double v5 = predict(sprime5);
    //     double v6 = predict(sprime6);
    //     double v7 = predict(sprime7);
    //     double v8 = predict(sprime8);
    //     p[0] = a1*b1*v0;
    //     p[1] = a2*b1*v1;
    //     p[2] = a3*b1*v2;
    //     p[3] = a1*b2*v3;
    //     p[4] = a2*b2*v4;
    //     p[5] = a3*b2*v5;
    //     p[6] = a1*b3*v6;
    //     p[7] = a2*b3*v7;
    //     p[8] = a3*b3*v8;
    //     printf("\n");
    // }
    return out;
}

inline double D_objective_d_p0(double *state, double p_other, const std::vector<double> &actions){
    double a1 = exp(2*(g(state[0]) - actions[0]));
    double a2 = exp(g(state[0]) - actions[0]);
    double a3 = 1 + a2 + exp(g(state[1]) - p_other);
    return -(M*(actions[0]-c)*(a1/(a3*a3)-a2/a3) + M*a2/a3);
}

inline double D_objective_d_x0(double *state, double policy_other, const std::vector<double> &actions){
    return -(-1.0 + beta*D_expectation_d_x0(state, policy_other, actions));
}

double objective(const std::vector<double> &actions, std::vector<double> &grad, void *data){
    objective_data *data_struct = (objective_data *) data;
    double *state_ptr = data_struct->state;
    double policy_other = data_struct->policy_other;
    double p_other = data_struct->p_other;
    if(!grad.empty()){
        grad[0] = D_objective_d_p0(state_ptr, p_other, actions);
        grad[1] = D_objective_d_x0(state_ptr, policy_other, actions);
    }
    double pi = profit(state_ptr, p_other, actions);
    double out = -(pi - actions[1] + beta*(expectation(state_ptr, policy_other, actions)));
    if(isnan(out)){
	printf("state = (%f, %f)\n", state_ptr[0], state_ptr[1]);
    }
    func_count++;
    return out;
}
