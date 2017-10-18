#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd>& estimations,
							  const vector<VectorXd>& ground_truth) {
	VectorXd rmse(4);
	rmse << 0, 0, 0, 0;

	if (estimations.size() != ground_truth.size() || estimations.size() == 0) {
		cout << "ERROR: invalid data to calculate RMSE" << endl;
		return rmse;
	}

	//accumulate squared residuals
	for (unsigned int i = 0; i < estimations.size(); ++i) {
		rmse += (estimations[i] - ground_truth[i]).array().pow(2).matrix();
	}

	rmse /= (const double)estimations.size();
	rmse = rmse.array().sqrt();
	return rmse;
}

VectorXd Tools::PolarToCartesian(const double& rho, const double& phi, const double& rhodot) {
	double px = rho * cos(phi);
	double py = rho * sin(phi);
	double vx = rhodot * sin(phi);
	double vy = rhodot * cos(phi);
	double v = sqrt(vx * vx + vy * vy);

	VectorXd cartesian(5);
	cartesian << px, py, v, 0, 0;
	return cartesian;
}

VectorXd Tools::CartesianToPolar(const double& px, const double& py, const double& vx, const double& vy) {
	double rho = sqrt(px * px + py * py);
	double phi = atan2(py, px);
	double rho_dot = (abs(rho) < 0.00001) ? 0 : (px * vx + py * vy) / rho;

	VectorXd polar(3);
	polar << rho, phi, rho_dot;
	return polar;
}

double Tools::ConstrainRadian(double x) {
	if (x > M_PI)
		return fmod(x - M_PI, 2 * M_PI) - M_PI;
	else if (x < -M_PI)
		return fmod(x + M_PI, 2 * M_PI) + M_PI;
	return x;
}
