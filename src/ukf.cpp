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
	P_ <<
	   1, 0, 0, 0, 0,
	   0, 1, 0, 0, 0,
	   0, 0, 1, 0, 0,
	   0, 0, 0, 1, 0,
	   0, 0, 0, 0, 1;

	// Process noise standard deviation longitudinal acceleration in m/s^2
	std_a_ = 0.2;

	// Process noise standard deviation yaw acceleration in rad/s^2
	std_yawdd_ = 0.2;

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
	n_aug_ = 7;

	// Sigma point spreading parameter
	lambda_ = 3 - n_x_;

	// process noise covariance matrix
	Q_.fill(0);
	Q_(0, 0) = pow(std_a_, 2);
	Q_(1, 1) = pow(std_yawdd_, 2);

	// augmented mean vector
	VectorXd x_aug_ = VectorXd(7);
	x_aug_.fill(0);

	// augmented state covariance
	P_aug_ = MatrixXd(n_aug_, n_aug_);
	P_aug_.fill(0);

	// augmented sigma point matrix
	Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);
	Xsig_aug_.fill(0);

	// predicted sigma points
	Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
	Xsig_pred_.fill(0);


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
			measurement_pack.raw_measurements_[1],		   // py
			0,                                             // v (unknown)
			0,                                             // psi (yaw angle) (unknown)
			0;                                             // psi dot (unknown)
		}


		// debug print x vector
		cout << "x_ initialized to: " << x_ << endl;

		// initialize timestamp
		previous_timestamp_ = measurement_pack.timestamp_;

		// done initializing, no need to predict or update
		is_initialized_ = true;
		return;
	}

	// calculate elapsed time in seconds
	double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
	previous_timestamp_ = measurement_pack.timestamp_;

	Prediction(dt);
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
	x_aug_.block(0, 0, n_x_, 1) = x_;
	GenerateAugmentedSigmaPoints();
	PredictSigmaPoints(delta_t);
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

void UKF::GenerateAugmentedSigmaPoints() {
	// create augmented covariance matrix
	P_aug_.fill(0);								// TODO: is this needed? or does it carry from previous state?
	P_aug_.block(0, 0, n_x_, n_x_) = P_;
	P_aug_.block(n_x_, n_x_, 2, 2) = Q_;
	//cout << "P_aug" << endl << P_aug_ << endl;

	//create square root matrix
	MatrixXd A = P_aug_.llt().matrixL();

	//create augmented sigma points
	Xsig_aug_.col(0) = x_aug_;

	for (int i = 0; i < n_aug_; ++i) {
		Xsig_aug_.col(i + 1) = x_aug_ + sqrt(lambda_ + n_aug_) * A.col(i);
		Xsig_aug_.col(n_aug_ + i + 1) = x_aug_ - sqrt(lambda_ + n_aug_) * A.col(i);
	}
}

void UKF::PredictSigmaPoints(const double& delta_t) {
	int n_col = 2 * n_aug_ + 1;

	for (int i = 0; i < n_col; ++i) {
		VectorXd x_k = Xsig_aug_.col(i);
		//cout << "x_k =" << endl << x_k << endl;
		double v = x_k(2);
		double psi = x_k(3);
		double psi_dot = x_k(4);
		double nu_a = x_k(5);
		double nu_psi_dotdot = x_k(6);

		VectorXd a = VectorXd(n_x_);
		VectorXd b = VectorXd(n_x_);

		if (abs(psi_dot) < 0.000001) {
			a << v / psi_dot*(sin(psi + psi_dot * delta_t) - sin(psi)),
			v / psi_dot*(-cos(psi + psi_dot * delta_t) + cos(psi)),
			0,
			psi_dot *delta_t,
			0;
		} else {
			a << v *cos(psi) * delta_t,
			v *sin(psi) * delta_t,
			0,
			0,
			0;
		}

		b << 0.5 * pow(delta_t, 2)*cos(psi)*nu_a,
		0.5 * pow(delta_t, 2)*sin(psi)*nu_a,
		delta_t *nu_a,
		0.5 * pow(delta_t, 2)*nu_psi_dotdot,
		delta_t *nu_psi_dotdot;

		Xsig_pred_.col(i) = x_k.segment(0, n_x_) + a + b;

	}
}

