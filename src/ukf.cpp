#include "ukf.h"
#include "Eigen/Dense"
#include <cmath>
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1.5;
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

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
  
  /**
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */

  // number of dimensions of state vector
  n_x_ = 5;

  // number of dimensions of augmented state vector
  n_aug_ = n_x_ + 2;

  // Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  // Number of sigma points for n_aug
  Xsig_dim = 2 * n_aug_ + 1;

  
  // Weights of sigma points
  weights_ = VectorXd(Xsig_dim);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (size_t i = 1; i < Xsig_dim; i++)
      weights_(i) = 0.5 / (lambda_ + n_aug_);

  // is_initialised should be initialised to false, until a measurement comes in
  is_initialized_ = false;

}

UKF::~UKF() {}

double UKF::LimitAngles_in_360 (double angle)
{
    while(angle > M_PI) { angle -= (2. * M_PI); }
    while(angle < -M_PI) { angle += (2. * M_PI); }
    return angle;
}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */ 
    if (!is_initialized_)
    {
        is_initialized_ = true;
        previous_timestamp = meas_package.timestamp_;
        double first_element{ meas_package.raw_measurements_(0) };
        double second_element{ meas_package.raw_measurements_(1) };

        if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
        {
            double third_element{ meas_package.raw_measurements_(2) };
            double velocity{ third_element };
            x_ << first_element * cos(second_element),
                first_element* sin(second_element),
                velocity,
                0.0,
                0.0;

            P_ << pow(std_radr_,2)*pow(std_radphi_,2), 0.0, 0.0, 0.0, 0.0,
                0, pow(std_radr_,2)*pow(std_radphi_,2), 0.0, 0.0, 0.0,
                0., 0., 1., 0., 0.,
                0., 0., 0., 1., 0.,
                0., 0., 0., 0., 1.;
        }
        else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
            x_ << first_element,
                second_element,
                0.,
                0.,
                0.;
            P_ << pow(std_laspx_, 2), 0., 0., 0., 0.,
                0., pow(std_laspy_, 2), 0., 0., 0.,
                0., 0., 2., 0., 0.,
                0., 0., 0., 1., 0.,
                0., 0., 0., 0., 1.;
        }
        return;
    }

    double dt = (meas_package.timestamp_ - previous_timestamp)/1000000.0;
    previous_timestamp = meas_package.timestamp_;
    Prediction(dt);

    if (meas_package.sensor_type_ == MeasurementPackage::LASER)
        UpdateLidar(meas_package);

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
        UpdateRadar(meas_package);

}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estim^12345678�0p+
   <-.  
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */
    // To include noise, x_aug is used for generating sigma points
    VectorXd x_aug_ = VectorXd(n_aug_);
    x_aug_.head(5) = x_;
    x_aug_(5) = 0.0;
    x_aug_(6) = 0.0;

    // Including the variance of the noises in Covariance matrix
    MatrixXd P_aug_ = MatrixXd(n_aug_, n_aug_);
    P_aug_.fill(0.0);
    P_aug_.topLeftCorner(5, 5) = P_;
    P_aug_.bottomRightCorner(2, 2) << pow(std_a_, 2), 0,
                                      0, pow(std_yawdd_, 2);

    // Creating a matrix for Sigma points
    MatrixXd Xsig = MatrixXd(n_aug_, Xsig_dim);

    // Sigma point generation
    Xsig.fill(0.0);
    Xsig.col(0) = x_aug_;
    MatrixXd A = P_aug_.llt().matrixL();

    for (size_t i = 0; i < n_aug_; i++)
    {
        Xsig.col(i+1) = x_aug_ + sqrt(lambda_+ n_aug_) * A.col(i);
        Xsig.col(i+1+n_aug_) = x_aug_ - sqrt(lambda_+ n_aug_) * A.col(i);
    }

    // Plug-in the sigma points in the process model to get the predicted sigma points

    // Matrix to hold the predicted sigma points
    Xsig_pred_ = MatrixXd(n_x_, Xsig_dim);
    for (size_t i = 0; i < Xsig_dim; i++)
    {
        double v = Xsig(2, i);
        double si = Xsig(3, i);
        double si_dot = Xsig(4, i);
        double nu_a = Xsig(5,i);
        double nu_yawdd = Xsig(6,i);
        VectorXd diff_term = VectorXd(n_x_);
        if (fabs(si_dot) > 0.001)
        {
            diff_term << (v / si_dot) * (sin(si + si_dot * delta_t) - sin(si)),
                (v / si_dot)* (-cos(si + si_dot * delta_t) + cos(si)),
                0.0,
                si_dot* delta_t,
                0.0;
        }
        else {
            diff_term << v * cos(si) * delta_t,
                v* sin(si)* delta_t,
                0.0,
                si_dot* delta_t,
                0.0;
        }
        VectorXd noise_influence = VectorXd(n_x_);
        noise_influence << 0.5 * pow(delta_t, 2) * cos(si) * nu_a,
            0.5 * pow(delta_t, 2) * sin(si) * nu_a,
            delta_t* nu_a,
            0.5 * pow(delta_t, 2) * nu_yawdd,
            delta_t* nu_yawdd;

        Xsig_pred_.col(i) = Xsig.col(i).head(5) + diff_term + noise_influence;
    }

    x_.fill(0.0);
    for (size_t i = 0; i < Xsig_dim; i++)
    {
        x_ = x_ + weights_(i) * Xsig_pred_.col(i);
        
    }

    P_.fill(0.0);
    for (size_t i = 0; i < Xsig_dim; i++)
    {
        VectorXd diff_term = Xsig_pred_.col(i) - x_;
        diff_term(3) = LimitAngles_in_360(diff_term(3));
        P_ = P_ + weights_(i) * diff_term * diff_term.transpose();
    }
}



void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
    MatrixXd H = MatrixXd(2, n_x_);
    H << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0;

    MatrixXd R = MatrixXd(2, 2);
    R << pow(std_laspx_, 2), 0,
        0,pow(std_laspy_, 2);

    VectorXd z_pred = VectorXd(2, 1);
    z_pred = H * x_;

    VectorXd y = VectorXd(2, 1);
    y = meas_package.raw_measurements_ - z_pred;

    MatrixXd S = MatrixXd(2, 2);
    S = H * P_ * H.transpose() + R;

    MatrixXd K_gain = MatrixXd(n_x_, 2);
    K_gain = P_ * H.transpose() * S.inverse();

    VectorXd dx_ = K_gain * y;
    x_ = x_ + dx_;

    MatrixXd I = MatrixXd::Identity(5, 5);
    P_ = (I - K_gain * H) * P_;
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
    
    MatrixXd Zsig = MatrixXd(3, Xsig_dim);
    for (size_t i = 0; i < Xsig_dim; i++)
    {
        //std::cout << Xsig_pred_ << std::endl;
        double px = Xsig_pred_(0,i);
        double py = Xsig_pred_(1,i);
        double v = Xsig_pred_(2,i);
        double si = Xsig_pred_(3,i);

        double dist = sqrt(px * px + py * py);
        Zsig.col(i) << dist,
            atan2(py, px),
            (px * cos(si) * v + py * sin(si) * v) / dist;  
    }
    VectorXd z_pred = VectorXd(3);
    z_pred.fill(0.0);
    for (size_t i = 0; i < Xsig_dim; i++)
    {
        z_pred = z_pred + weights_(i) * Zsig.col(i);
        
    }

    MatrixXd S = MatrixXd(3, 3);
    S.fill(0.0);
    for (size_t i = 0; i < Xsig_dim; i++)
    {
        VectorXd diff_term = Zsig.col(i) - z_pred;
        diff_term(1) = LimitAngles_in_360(diff_term(1));
        S = S + weights_(i) * diff_term * diff_term.transpose();
    }
    S(0, 0) += pow(std_radr_, 2);
    S(1, 1) += pow(std_radphi_, 2);
    S(2, 2) += pow(std_radrd_, 2);

    MatrixXd T_ = MatrixXd(n_x_, 3);
    T_.fill(0.0);
    for (size_t i = 0; i < Xsig_dim; i++)
    {
        VectorXd diff_X = Xsig_pred_.col(i) - x_;
        diff_X(3) = LimitAngles_in_360(diff_X(3));
        VectorXd diff_Z = Zsig.col(i) - z_pred;
        diff_Z(1) = LimitAngles_in_360(diff_Z(1));
        
        T_ = T_ + weights_(i) * diff_X * diff_Z.transpose();
    }

    MatrixXd K_gain = MatrixXd(n_x_, 3);
    K_gain = T_ * S.inverse();

    VectorXd dz_ = (meas_package.raw_measurements_ - z_pred);
    dz_(1) = LimitAngles_in_360(dz_(1));

    x_ = x_ + K_gain * dz_;

    P_ = P_ - K_gain * S * K_gain.transpose();
}