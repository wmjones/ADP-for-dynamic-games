#include <math.h>
#include <stdio.h>
#include <fstream>
#include <eigen3/Eigen/Dense>
#include "parameters.hpp"
#include "fitting.hpp"

#include <dlib/dnn.h>
#include <dlib/data_io.h>
using namespace std;
using namespace dlib;
using net_type = loss_mean_squared<fc<1,
				      relu<fc<10,
				      relu<fc<10,
				      relu<fc<10,
				      relu<fc<10,
				      relu<fc<10,
				      // htan<fc<100,
				      input<matrix<float,0,1>>
					      >>>>>>>>>>>>;
net_type net;
// double v_min;
// double v_max;

// void vander_hermite(double **xy_data, double **X, int ** alpha, size_t num_of_alpha){
//     double z[d];
//     for(size_t i=0; i<d; i++){
//         z[i] = 0.0;
//     }
//     for(size_t i=0; i<S; i++){
//         for(size_t j=0; j<num_of_alpha; j++){
//             for(size_t k=0; k<d; k++){
//                 z[k] = (2*xy_data[i][k] - xmin[k] - xmax[k])/(xmax[k] - xmin[k]);
//             }
//             X[i][j] = T_alpha(z, alpha[j]);
//         }
//     }
//     for(size_t i=0; i<S; i++){
//         for(size_t j=0; j<num_of_alpha; j++){
//             for(size_t k=0; k<d; k++){
//                 z[k] = (2*xy_data[i][k] - xmin[k] - xmax[k])/(xmax[k] - xmin[k]);
//             }
//             X[S + i][j] = dT_alpha(z, 0, alpha[j]);
//         }
//     }
//     for(size_t i=0; i<S; i++){
//         for(size_t j=0; j<num_of_alpha; j++){
//             for(size_t k=0; k<d; k++){
//                 z[k] = (2*xy_data[i][k] - xmin[k] - xmax[k])/(xmax[k] - xmin[k]);
//             }
//             X[2*S + i][j] = dT_alpha(z, 1, alpha[j]);
//         }
//     }
// }

void vander_lagrange(double **xy_data, double **X, int **alpha, size_t num_of_alpha){
    double z[d];
    for(size_t i=0; i<d; i++){
        z[i] = 0.0;
    }
    for(size_t i=0; i<S; i++){
        for(size_t j=0; j<num_of_alpha; j++){
            for(size_t k=0; k<d; k++){
                z[k] = (2*xy_data[i][k] - xmin[k] - xmax[k])/(xmax[k] - xmin[k]);
            }
            X[i][j] = T_alpha(z, alpha[j]);
        }
    }
}


double predict(double *state){
    // if(approx_type == "cheb_hermite"){
    // 	return V_hat(state, value_coef, alpha, num_of_coef)
    // }
    if(approx_type == 0){
	return V_hat(state, value_coef, alpha_coef, num_of_coef);
    }
    else if(approx_type == 1){
	return V_hat(state, value_coef, alpha_coef, num_of_coef);
    }
    else if(approx_type == 2){
	std::vector<matrix<float, 0, 1>> sample(1);
	sample[0] = {(state[0] - xmin[0])/(xmax[0] - xmin[0])*2-1,
		     (state[1] - xmin[1])/(xmax[1] - xmin[1])*2-1};
	auto out = net(sample);
	if(isnan(out[0])){
	    printf("ann is nan at state = (%f, %f)\n", state[0], state[1]);
	}
	// return out[0]*(v_max-v_min)+v_min;
	return out[0];
	// printf("%f\n", out[0]);
	// std::cout << out[0] << std::endl;
	// return 0;
    }
    else{
	return 0;
    }
}


void fitting(double **xy_data, double *v_hat, double *dz_data0, double *dz_data1, double *b, double **z_knots){
    // if(approx_type == "cheb_hermite"){
    // 	// void fitting_hermite(double **xy_data, double *v_hat,
    // 	// double *dz_data0, double *dz_data1, double *b,  int **alpha, size_t num_of_coef)
    // 	Eigen::MatrixXd data(3*S, num_of_coef);
    // 	Eigen::VectorXd y = Eigen::VectorXd::Zero(3*S);
    // 	Eigen::VectorXd coef = Eigen::VectorXd::Zero(num_of_coef);

    // 	for(size_t i=0; i<S; i++){
    // 	    y[i] = v_hat[i];
    // 	}
    // 	for(size_t i=0; i<S; i++){
    // 	    y[S + i] = dz_data0[i];
    // 	}
    // 	for(size_t i=0; i<S; i++){
    // 	    y[2*S + i] = dz_data1[i];
    // 	}
    // 	double **X = new double*[3*S];
    // 	for(size_t i=0; i<3*S; i++){
    // 	    X[i] = new double[num_of_coef];
    // 	}
    // 	vander_hermite(xy_data, X, alpha, num_of_coef);
    // 	for(size_t i=0; i<3*S; i++){
    // 	    for(size_t j=0; j<num_of_coef; j++){
    // 		data(i, j) =  X[i][j];
    // 	    }
    // 	}
    // 	coef = (data.transpose() * data).ldlt().solve(data.transpose() * y);
    // 	//    coef = data.colPivHouseholderQr().solve(y);
    // 	for(size_t i=0; i<num_of_coef; i++){
    // 	    b[i] = coef[i];
    // 	}
    // 	for(size_t i=0; i<3*S; i++){
    // 	    delete[] X[i];
    // 	}
    // 	delete[] X;
    // }

    if(approx_type == 0){
	// void fitting_lagrange(double **xy_data, double *v, double *b,
        //               int ** alpha, size_t num_of_coef)
	Eigen::MatrixXd data(S, num_of_coef);
	Eigen::VectorXd y = Eigen::VectorXd::Zero(S);
	Eigen::VectorXd coef = Eigen::VectorXd::Zero(num_of_coef);
	for(size_t i=0; i<S; i++){
	    y(i) = v_hat[i];
	}
	double **X = new double*[S];
	for(size_t i=0; i<S; i++){
	    X[i] = new double[num_of_coef];
	}
	vander_lagrange(xy_data, X, alpha_coef, num_of_coef);
	for(size_t i=0; i<S; i++){
	    for(size_t j=0; j<num_of_coef; j++){
		data(i, j) =  X[i][j];
	    }
	}
    //    coef = (data.transpose() * data).ldlt().solve(data.transpose() * y);
	coef = data.colPivHouseholderQr().solve(y);
	for(size_t i=0; i<num_of_coef; i++){
	    b[i] = coef(i);
	}

	for(size_t i=0; i<S; i++){
	    delete[] X[i];
	}
	delete[] X;
    }

    else if(approx_type == 1){
	// void fitting_Cai(double **xy_data, double **z_knots, double *v_hat, double *b,
        //          int **alpha, size_t num_of_coef)
	for(size_t i=0; i<num_of_coef; i++){
	    double out = 0;
	    for(size_t j=0; j<S; j++){
		out+=v_hat[j]*T_alpha(z_knots[j], alpha_coef[i]);
	    }
	    double d_tilde = 0;
	    for(size_t j=0; j<d; j++){
		d_tilde+=(alpha_coef[i][j]>0)?1:0;
	    }
	    b[i] = out/((double) pow(num_of_knots, d))*pow(2, d_tilde);
	}
    }

    else if(approx_type == 2){
	net.clean();		// not sure if i should use this or not
	dlib::rand rnd(time(0));
	std::vector<matrix<float, 0, 1>> samples(S);
	std::vector<float> labels(S);

	// v_max = *std::max_element(v_hat, v_hat+S);
	// v_min = *std::min_element(v_hat, v_hat+S)-1e-5;
	for(size_t i=0; i<S; i++){
	    samples[i] = {(xy_data[i][0] - xmin[0])/(xmax[0] - xmin[0])*2-1,
			  (xy_data[i][1] - xmin[1])/(xmax[1] - xmin[1])*2-1};
	    // labels[i] = {(v_hat[i]-v_min)/(v_max-v_min)};
	    labels[i] = {v_hat[i]};
	    // std::cout << "i=" << i << "\tsample=(" << samples[i](0,0) << ", "
	    // 	      << samples[i](1,0) << ")\t\tlabel=" << v_hat[i] << std::endl;
	}
	// printf("(%f, %f)\n", v_max, v_min);
	// dnn_trainer<net_type> trainer(net);
	dnn_trainer<net_type, sgd> trainer(net, sgd(0.0005, 0.9));
	trainer.set_learning_rate(0.01);

	for(int i=0; i<10000; i++)
	    trainer.train_one_step(samples, labels);
	// trainer.train(samples, labels);
	trainer.get_net();
	// for(size_t i=0; i<S; i++){
        //     std::cout << "i=" << i << "\tsample=(" << samples[i](0,0) << ", "
	//     	      << samples[i](1,0) << ")\t\tlabel=" << labels[i]
	// 	      << "\tpredicted = " << net(samples[i]) << std::endl;
	// }
	// samples.clear();
        // labels.clear();
	// sample[0] = {1, 1};
	// auto out = net(sample);
	// return out[0];
	// printf("%f\n", out[0]);
	// std::cout << out[0] << std::endl;
    }
}
