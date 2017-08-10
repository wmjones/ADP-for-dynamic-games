#include <iostream>
#include <math.h>
#include <stdio.h>
#include <nlopt.hpp>
#include <fstream>
#include <vector>
#include "fitting.hpp"
#include "parameters.hpp"
#include "model.hpp"
#include <omp.h>
#include <stdio.h>

using namespace std;
ofstream myfile;

void bellman(double **xy_knots, double *value, double *value_last, double *policy, double *policy_last, double *price, double *price_last, double *dvalue0, double *dvalue1, size_t nmax){
// #pragma omp parallel for
    for(size_t i=0; i<S; i++){
	// printf("Optimization[%zd]\n", i);
	// printf("state = (%f, %f)\n", xy_knots[i][0], xy_knots[i][1]);
	size_t index_other = i/num_of_knots + (i%num_of_knots)*num_of_knots;
	double p_other = price[index_other];
	double policy_other = policy[index_other];
        objective_data data_init = {xy_knots[i], p_other, policy_other};
//        if(false){
//            double h = 1e-5;
//            double state[2] = {xy_knots[i][0], xy_knots[i][1]};
//            double state_x[2] = {xy_knots[i][0] + h, xy_knots[i][1]};
//            double state_y[2] = {xy_knots[i][0], xy_knots[i][1] + h};
//            vector<double> grad(2);
//            vector<double> actions(2);
//            actions[0] = 6;
//            actions[1] = 0;
//            vector<double> actions_x(2);
//            actions_x[0] = 6+h;
//            actions_x[1] = 0;
//            vector<double> actions_y(2);
//            actions_y[0] = 6;
//            actions_y[1] = 0+h;
//            double v = objective(actions, grad, &data_init);
//            double v_x = objective(actions_x, grad, &data_init);
//            double v_y = objective(actions_y, grad, &data_init);
//            double dv_dx = D_objective_d_p0(state, actions);
//            double dv_dy = D_objective_d_x0(state, actions);
//            double approx_dv_dx = (v_x-v)/h;
//            double approx_dv_dy = (v_y-v)/h;
//            printf("approx grad\n");
//            printf("(%.8f, %.8f)\n", approx_dv_dx, approx_dv_dy);
//            printf("grad\n");
//            printf("(%.8f, %.8f)\n", dv_dx, dv_dy);
//            printf("\n");
//        }
        nlopt::opt opt(nlopt::LD_LBFGS, 2);
//        LD_LBFGS
//        LN_COBYLA
        opt.set_xtol_rel(1e-5);
	std::vector<double> x_l(2);
        x_l[0] = 0; x_l[1] = 0;
        opt.set_lower_bounds(x_l);
        opt.set_min_objective(objective, &data_init);
	std::vector<double> actions(2);
        actions[0] = price_last[i]+.001; actions[1] = policy_last[i]+.001;
	// actions[0] = 6.5; actions[1] = .05;
        double minf;
        nlopt::result result = opt.optimize(actions, minf);
        if(isnan(-minf) || isnan(actions[0]) || isnan(actions[1]) || -minf>1000 || actions[1]<0){
            printf("Error in Optimization\n");
        }
        value[i] = -minf;
        price[i] = actions[0];
        policy[i] = actions[1];
        if(result != 1 && result != 4 && result != 3){
            printf("NLOPT ERROR CODE %d v=%f p=%f x=%f\n", result, -minf, actions[0], actions[1]);
            value[i] = value_last[i];
            price[i] = price_last[i];
            policy[i] = policy_last[i];
        }
    }
}

void plotting(double **xy_knots, double *value, double *price, double *policy, size_t k){
    // printf("Plotting\n");
    double *x_test, *y_test, *v_test, *p_test, *inv_test;
    x_test = new double [num_of_test];
    y_test = new double [num_of_test];
    v_test = new double [S_test];
    p_test = new double [S_test];
    inv_test = new double [S_test];
    double **xy_test = new double*[S_test];
    for(size_t i=0; i<S_test; i++){
        xy_test[i] = new double[d];
    }
    for(size_t i=0; i<num_of_test; i++){
        x_test[i] = xmin[0] + i/double(num_of_test-1)*(xmax[0] - xmin[0]);
    }
    for(size_t i=0; i<num_of_test; i++){
        y_test[i] = xmin[1] + i/double(num_of_test-1)*(xmax[1]-xmin[1]);
    }
    for(size_t i=0; i<S_test; i++){
        xy_test[i][0] = x_test[i % num_of_test];
        xy_test[i][1] = y_test[i / num_of_test];
    }
    for(size_t i=0; i<S_test; i++){
	v_test[i] = predict(xy_test[i]);
	if(i<S){
	    p_test[i] = price[i];
	    inv_test[i] = policy[i];
	}
	else{
	    p_test[i] = 0;
	    inv_test[i] = 0;
	}
    }
    for(size_t i=0; i<S_test; i++){
        double x1 = 0;
        double y1 = 0;
        double v1 = 0;
        if(i<S){
            x1 = xy_knots[i][0];
            y1 = xy_knots[i][1];
            v1 = value[i];
        }
        myfile << xy_test[i][0] << "," << xy_test[i][1] << "," <<
        v_test[i] << "," << p_test[i] << "," << inv_test[i] << "," << x1 << "," << y1 << "," << v1 << "\n";
    }
    if(k+1==Nmax){
        myfile.close();
    }
    for(size_t i=0; i<S_test; i++){
        delete[] xy_test[i];
    }
    delete[] inv_test;
    delete[] xy_test;
    delete[] x_test;
    delete[] y_test;
    delete[] v_test;
    delete[] p_test;
}

int main(int argv, char* argc[]){
    auto start = std::chrono::high_resolution_clock::now();
    myfile.open("model_data.csv");
    myfile << "x,y,v,p,inv,x1,y1,v1\n";
    double *x_knots, *y_knots, *value, *policy, *price, *dvalue0, *dvalue1, *z;
    x_knots = new double [num_of_knots];
    y_knots = new double [num_of_knots];
    z = new double [num_of_knots];
    value = new double [S];
    policy = new double [S];
    price = new double [S];
    dvalue0 = new double [S];
    dvalue1 = new double [S];
    double **xy_knots = new double *[S];
    double **z_knots = new double *[S];
    for(size_t i=0; i<S; i++){
        xy_knots[i] = new double [d];
        z_knots[i] = new double [d];
	dvalue0[i] = 0;
	dvalue1[i] = 0;
	value[i] = 0;		// inital value function
        policy[i] = 0;		// inital investment guess
        price[i] = 6.5;		// initial price guess
    }
    for(size_t i=0; i<num_of_knots; i++){
        z[i] = -cos((2*(i+1)-1)*M_PI/(2*num_of_knots));
	if(approx_type==2){
	    x_knots[i] = (z[i]+1)*(xmax[0]-xmin[0])/2+xmin[0];
	    y_knots[i] = (z[i]+1)*(xmax[1]-xmin[1])/2+xmin[1];
	    // x_knots[i] = xmin[0] + i/double(num_of_knots-1)*(xmax[0]-xmin[0]);
	    // y_knots[i] = xmin[1] + i/double(num_of_knots-1)*(xmax[1]-xmin[1]);
	}
	else{
	    x_knots[i] = (z[i]+1)*(xmax[0]-xmin[0])/2+xmin[0];
	    y_knots[i] = (z[i]+1)*(xmax[1]-xmin[1])/2+xmin[1];
	}
	// printf("(%f, %f)\n", x_knots[i], y_knots[i]);
    }
    for(size_t i=0; i<S; i++){
        xy_knots[i][0] = x_knots[i % num_of_knots];
        xy_knots[i][1] = y_knots[i / num_of_knots];
        z_knots[i][0] = z[i % num_of_knots];
        z_knots[i][1] = z[i / num_of_knots];
    }
    double *value_last, *price_last, *policy_last;
    value_last = new double [S];
    price_last = new double [S];
    policy_last = new double [S];
    fitting(xy_knots, value, dvalue0, dvalue1, value_coef, z_knots);
    plotting(xy_knots, value, price, policy, -1);
    for(size_t k=0; k<Nmax; k++){
	for(size_t i=0; i<S; i++){
	    value_last[i] = value[i];
	    price_last[i] = price[i];
	    policy_last[i] = policy[i];
	}
        bellman(xy_knots, value, value_last, policy, policy_last, price, price_last, dvalue0, dvalue1, k);
        for(size_t i=0; i<S; i++){
            value[i] = value[i]*eta + (1-eta)*value_last[i];
            price[i] = price[i]*eta + (1-eta)*price_last[i];
            policy[i] = policy[i]*eta + (1-eta)*policy_last[i];
        }
	if(approx_type==2){
	    fitting(xy_knots, value, dvalue0, dvalue1, value_coef, z_knots);
	}
	else{
	    fitting(xy_knots, value, dvalue0, dvalue1, value_coef, z_knots);
	    fitting(xy_knots, price, dvalue0, dvalue1, price_coef, z_knots);
	    fitting(xy_knots, policy, dvalue0, dvalue1, policy_coef, z_knots);
	}
        plotting(xy_knots, value, price, policy, k);
        printf("%zd\n", k);
    }
    auto end = std::chrono::high_resolution_clock::now();
    printf("Finished Iterations\n");
    printf("Objective Evaluated %f times per state\n", func_count/((double) (Nmax * S)));
    cout << "Time " << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count()/1000.0 << endl;

    delete[] x_knots;
    delete[] y_knots;
    delete[] value;
    delete[] value_last;
    delete[] policy;
    delete[] policy_last;
    delete[] price;
    delete[] price_last;
    delete[] dvalue0;
    delete[] dvalue1;
    delete[] z;
    for(size_t i=0; i<S; i++){
        delete[] xy_knots[i];
	delete[] z_knots[i];
    }
    delete[] xy_knots;
    for(size_t i=0; i<num_of_coef; i++){
        delete[] alpha_coef[i];
    }
    delete[] alpha_coef;
    delete[] xmin;
    delete[] xmax;
    delete[] value_coef;
    delete[] policy_coef;

    return 0;
}
