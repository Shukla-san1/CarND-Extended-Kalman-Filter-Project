#include "kalman_filter.h"
#include<iostream>


using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
     predict the state
  */
	x_ = F_ * x_;
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /*
     update the state by using Kalman Filter equations
  */
	VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred;
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}

//normalize the angle

/**

  * This method is to normalized the angle, so that it should remain between pi and -pi
  * Input: take angle as input
  * output: return the normalized angle
  *
*/

float KalmanFilter::Normalize_angle(float angle){
	if(fabs(angle) > PI){
		angle -= round(angle / (2. * PI)) * (2.* PI);
	}

	return angle;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {


  /**

    * update the state by using Extended Kalman Filter equations
  */
	float px = x_(0);
    float py = x_(1);
    float vx = x_(2);
    float vy = x_(3);

    float rho_dot;
    float rho = sqrtf(powf(px, 2) + powf(py, 2));
    float phi;

    //Division by Zero check

    if(fabs(px) < 0.0001){
        cout << "Error while converting vector x_ to polar coordinates: Division by Zero" << endl;
      }
      else {
        phi = atan2(py, px);
      }




    //Division by Zero check

    if (fabs(rho) < 0.0001) {

        cout << "Error while converting vector x_ to polar coordinates: Division by Zero" << endl;
      }
      else {

         rho_dot = (px*vx + py*vy)/rho;

      }

    VectorXd hx(3);
    hx << rho, phi, rho_dot;


	VectorXd y = z - hx;

//Normalizing the y(1) value.

	y(1) = Normalize_angle(y(1));

	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;
	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;

}


