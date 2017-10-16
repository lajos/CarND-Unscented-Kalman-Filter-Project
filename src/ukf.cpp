#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include "tools.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
	// if this is false, laser measurements will be ignored (except during init)
	use_laser_ = true;

	// if this is false, radar measurements will be ignored (except during init)
	use_radar_ = false;

	// initial state vector
	x_ = VectorXd(5);
	x_.fill(0);

	// initial covariance matrix
	P_ = MatrixXd(5, 5);
	P_.fill(0);

	// Process noise standard deviation longitudinal acceleration in m/s^2
	std_a_ = 30;

	// Process noise standard deviation yaw acceleration in rad/s^2
	std_yawdd_ = 30;

	// Laser measurement noise standard deviation position1 in m
	std_laspx_ = 0.15;

	// Laser measurement noise standard deviation position2 in m
	std_laspy_ = 0.15;

	// Radar measurement noise standard deviation radius in m
	std_radr_ = 0.3;

	// Radar measurement noise standard deviation angle in rad
	std_radphi_ = 0.03;

	// Radar measurement noise standard deviation radius change in m/s
	std_radrd_ = 0.3;

	// State dimension
	n_x_ = 5;

	// Augmented state dimension
	//int n_aug_;

	// Sigma point spreading parameter
	double lambda_ = 3 - n_x_;

	/**
	TODO:

	Complete the initialization. See ukf.h for other member properties.

	Hint: one or more values initialized above might be wildly off...
	*/
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} measurement_pack The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(const MeasurementPackage& measurement_pack) {

	// ignore measurement if not used
	if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR && !use_radar_) return;
	if (measurement_pack.sensor_type_ == MeasurementPackage::LASER && !use_laser_) return;

	/*****************************************************************************
	*  Initialization
	****************************************************************************/
	if (!is_initialized_) {
		// first measurement
		cout << "EKF: " << endl;

		if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
			// convert radar from polar to cartesian coordinates and initialize state.
			x_ << Tools::PolarToCartesian(measurement_pack.raw_measurements_[0],   // rho
										  measurement_pack.raw_measurements_[1],   // phi
										  measurement_pack.raw_measurements_[2]);  // rhodot
		} else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
			x_ << measurement_pack.raw_measurements_[0],   // px
			measurement_pack.raw_measurements_[1],   // py
			0,                                       // vx (unknown)
			0;                                       // vy (unknown)
		}


		// debug print x vector
		cout << "x_ initialized to: " << x_ << endl;

		// initialize timestamp
		previous_timestamp_ = measurement_pack.timestamp_;

		// done initializing, no need to predict or update
		is_initialized_ = true;
		return;
	}

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
	/**
	TODO:

	Complete this function! Estimate the object's location. Modify the state
	vector, x_. Predict sigma points, the state, and the state covariance matrix.
	*/
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(const MeasurementPackage& measurement_pack) {
	/**
	TODO:

	Complete this function! Use lidar data to update the belief about the object's
	position. Modify the state vector, x_, and covariance, P_.

	You'll also need to calculate the lidar NIS.
	*/
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(const MeasurementPackage& measurement_pack) {
	/**
	TODO:

	Complete this function! Use radar data to update the belief about the object's
	position. Modify the state vector, x_, and covariance, P_.

	You'll also need to calculate the radar NIS.
	*/
}
