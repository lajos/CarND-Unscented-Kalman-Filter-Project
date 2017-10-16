#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

class Tools {
public:
	/**
	* Constructor.
	*/
	Tools();

	/**
	* Destructor.
	*/
	virtual ~Tools();

	/**
	* A helper method to calculate RMSE.
	*/
	VectorXd CalculateRMSE(const vector<VectorXd>& estimations, const vector<VectorXd>& ground_truth);

	/**
	* Convert polar to cartesian coordinates.
	*/
	static VectorXd PolarToCartesian(const double& rho, const double& phi, const double& rhodot);

	/**
	* Convert cartesian to polar coordinates.
	*/
	static VectorXd CartesianToPolar(const double& px, const double& py, const double& vx, const double& vy);

	/**
	* Constrain radian angle between -M_PI and M_PI
	*/
	static double ConstrainRadian(double x);

};

#endif /* TOOLS_H_ */