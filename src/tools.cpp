#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
	VectorXd rmse(4);
	rmse << 0, 0, 0, 0;

	// check the validity of the following inputs:
	//  * the estimation vector size should not be zero
	//  * the estimation vector size should equal ground truth vector size
	if (estimations.size() != ground_truth.size()
		|| estimations.size() == 0) {
		cout << "Invalid estimation or ground_truth data" << endl;
		return rmse;
	}

	//accumulate squared residuals
	for (unsigned int i = 0; i < estimations.size(); ++i) {

		VectorXd residual = estimations[i] - ground_truth[i];

		//coefficient-wise multiplication
		residual = residual.array()*residual.array();
		rmse += residual;
	}

	//calculate the mean
	rmse = rmse / estimations.size();

	//calculate the squared root
	rmse = rmse.array().sqrt();

	//return the result
	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd &x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
	MatrixXd Hj(3, 4);
	Hj << 1, 1, 0, 0,
		1, 1, 0, 0,
		1, 1, 1, 1;
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	//pre-compute a set of terms to avoid repeated calculation
	float c1 = px * px + py * py;
	float c2 = sqrt(c1);
	float c3 = (c1*c2);

	//check division by zero
	if (fabs(c1) < 0.0001) {
		cout << "CalculateJacobian () - Error - Division by Zero" << endl;
		return Hj;
	}

	//compute the Jacobian matrix
	Hj << (px / c2), (py / c2), 0, 0,
		-(py / c1), (px / c1), 0, 0,
		py*(vx*py - vy * px) / c3, px*(px*vy - py * vx) / c3, px / c2, py / c2;

	return Hj;
}

VectorXd Tools::Calculateh(const VectorXd &x_state) {
	//Calculate h(x)
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	float rho = sqrt(px*px + py * py);
	float phi = atan2(py, px);

	float rho_dot;
	
	if (fabs(rho) < 0.0001) {
		rho_dot = 0;
	}
	else {
		rho_dot = (px*vx + py * vy) / rho;
	}
	VectorXd h = VectorXd(3);
	h << rho, phi, rho_dot;

	return h;
}

VectorXd Tools::ConvertPolarToCartesian(const VectorXd &raw_measurements_) {
	float rho = raw_measurements_[0]; // range
	float phi = raw_measurements_[1]; // bearing
	float rho_dot = raw_measurements_[2]; // velocity of rho
	// Coordinates convertion from polar to cartesian
	float x = rho * cos(phi);
	if (x < 0.0001) {
		x = 0.0001;
	}
	float y = rho * sin(phi);
	if (y < 0.0001) {
		y = 0.0001;
	}
	/*double vx = rho_dot * cos(phi);
	double vy = rho_dot * sin(phi);*/
	VectorXd x_ = VectorXd(4);
	x_ << x, y, 0, 0;

	return x_;
}
