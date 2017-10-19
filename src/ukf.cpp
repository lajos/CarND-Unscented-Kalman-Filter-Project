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
	use_laser_ = false;

	// if this is false, radar measurements will be ignored (except during init)
	use_radar_ = true;

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
	Q_ = MatrixXd(2, 2);
	Q_.fill(0);
	Q_(0, 0) = pow(std_a_, 2);
	Q_(1, 1) = pow(std_yawdd_, 2);

	// augmented mean vector
	x_aug_ = VectorXd(7);
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

	// sigma point weights
	weights_ = VectorXd(2 * n_aug_ + 1);

	// set weights
	weights_(0) = 1.0 * lambda_ / (lambda_ + n_aug_);
	for (int i = 1; i < 2 * n_aug_ + 1; ++i) {
		weights_(i) = 1.0 / (2 * (n_aug_ + lambda_));
	}

	//set measurement dimension, radar can measure r, phi, and r_dot
	n_z_radar_ = 3;

	// radar predicted mean
	z_pred_radar_ = VectorXd(n_z_radar_);

	// radar predicted covariance
	S_radar_ = MatrixXd(n_z_radar_, n_z_radar_);

	// radar measurement noise covariance matrix
	R_radar_ = MatrixXd(n_z_radar_, n_z_radar_);
	R_radar_ << std_radr_ *std_radr_, 0, 0,
			 0, std_radphi_ *std_radphi_, 0,
			 0, 0, std_radrd_ *std_radrd_;

	// sigma points in radar measurement space
	Zsig_radar_ = MatrixXd(n_z_radar_, 2 * n_aug_ + 1);

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} measurement_pack The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(const MeasurementPackage& measurement_pack) {

	if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
		cout << "process RADAR measurement" << endl;
	}

	if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
		cout << "process LASER measurement" << endl;
	}

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

	if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
		UpdateRadar(measurement_pack);
	} else {
		UpdateLidar(measurement_pack);
	}

	cout << "measurrement processed" << endl;
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
	PredictMeanAndCovariance();
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
	PredictRadarMeasurement();
	VectorXd z = VectorXd(n_z_radar_);
	z << measurement_pack.raw_measurements_[0],   // r
	measurement_pack.raw_measurements_[1],    // phi
	measurement_pack.raw_measurements_[2];    // r_dot

	UpdateRadarState(z);
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

		if (abs(psi_dot) > 0.000001) {
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

void UKF::PredictMeanAndCovariance() {
	// vector for predicted state
	VectorXd x = VectorXd(n_x_);

	// covariance matrix for prediction
	MatrixXd P = MatrixXd(n_x_, n_x_);

	//predict state mean
	for (int i = 0; i < n_x_; ++i) {
		x(i) = (weights_.array() * Xsig_pred_.row(i).transpose().array()).sum();
	}

	//predicted state covariance matrix
	P.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x;

		//angle normalization
		x_diff(3) = Tools::ConstrainRadian(x_diff(3));
		//while (x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
		//while (x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;

		P = P + weights_(i) * x_diff * x_diff.transpose();
	}

	x_ = x;
	P_ = P;
}


void UKF::PredictRadarMeasurement() {

	//transform sigma points into measurement space
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		double p_x = Xsig_pred_(0, i);
		double p_y = Xsig_pred_(1, i);
		double v = Xsig_pred_(2, i);
		double yaw = Xsig_pred_(3, i);

		double v1 = cos(yaw) * v;
		double v2 = sin(yaw) * v;

		// TODO: move this conversion to tools?
		// measurement model
		Zsig_radar_(0, i) = sqrt(p_x * p_x + p_y * p_y);                           // r
		Zsig_radar_(1, i) = atan2(p_y, p_x);                                       // phi
		Zsig_radar_(2, i) = (p_x * v1 + p_y * v2) / sqrt(p_x * p_x + p_y * p_y);   // r_dot
	}

	// mean predicted measurement
	VectorXd z_pred = VectorXd(n_z_radar_);
	z_pred.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		z_pred = z_pred + weights_(i) * Zsig_radar_.col(i);
	}

	//measurement covariance matrix S
	MatrixXd S = MatrixXd(n_z_radar_, n_z_radar_);
	S.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 sigma points
		//residual
		VectorXd z_diff = Zsig_radar_.col(i) - z_pred;

		//angle normalization
		z_diff(1) = Tools::ConstrainRadian(z_diff(1));
		//while (z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
		//while (z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;

		S = S + weights_(i) * z_diff * z_diff.transpose();
	}

	//add measurement noise covariance matrix
	S = S + R_radar_;

	z_pred_radar_ = z_pred;
	S_radar_ = S;
}

void UKF::UpdateRadarState(const VectorXd& z) {
	//create matrix for cross correlation Tc
	MatrixXd Tc = MatrixXd(n_x_, n_z_radar_);

	//calculate cross correlation matrix
	Tc.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		VectorXd z_diff = Zsig_radar_.col(i) - z_pred_radar_;

		//angle normalization
		z_diff(1) = Tools::ConstrainRadian(z_diff(1));
		//while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
		//while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;

		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_;

		//angle normalization
		x_diff(3) = Tools::ConstrainRadian(x_diff(3));
		//while (x_diff(3)> M_PI) x_diff(3) -= 2.*M_PI;
		//while (x_diff(3)<-M_PI) x_diff(3) += 2.*M_PI;

		Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
	}

	//Kalman gain K;
	MatrixXd K = Tc * S_radar_.inverse();

	//residual
	VectorXd z_diff = z - z_pred_radar_;

	//angle normalization
	z_diff(1) = Tools::ConstrainRadian(z_diff(1));
	//while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
	//while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;

	//update state mean and covariance matrix
	x_ = x_ + K * z_diff;
	P_ = P_ - K * S_radar_ * K.transpose();

}
