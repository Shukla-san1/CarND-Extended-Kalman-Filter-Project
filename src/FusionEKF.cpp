#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.F_ = MatrixXd(4, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;
  H_laser_ << 1,0,0,0,
		     0,1,0,0;


   // Init object covariance matrix.
  ekf_.P_ << 1, 0, 0,    0,
              0, 1, 0,    0,
              0, 0, 1000, 0,
              0, 0, 0,    1000;

   // Init state transition matrix. Indices (0, 2) and (1, 3) will be replaced
   // with dt at each measurement step.
   ekf_.F_ << 1, 0, 0, 0,
              0, 1, 0, 0,
              0, 0, 1, 0,
              0, 0, 0, 1;

	//set the acceleration noise components
	noise_ax = 9;
	noise_ay = 9;

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
      * This is function which predict and update based on sensor data. Below are different thing to be done.
      	  * Initialize the state ekf_.x_ with the first measurement.
      	  * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */

    // first measurement
       cout << "EKF: " << endl;
       ekf_.x_ = VectorXd(4);
       ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {



      //Convert radar data from polar to cartesian coordinates and initialize state.



    	float r = measurement_pack.raw_measurements_[0];
    	float theta = measurement_pack.raw_measurements_[1];
    	theta = ekf_.Normalize_angle(theta);
    	float rho_dot = measurement_pack.raw_measurements_[2];
    	float vx = rho_dot * cos(theta);
    	float vy = rho_dot * sin(theta);
    	float px = r*cos(theta);
    	float py = r*sin(theta);
    	if(fabs(px) < 0.0001){
    		px = 0.0001;
    	}
    	else if(fabs(py) < 0.0001){
    		py = 0.0001;
    	}

    	ekf_.x_ << px, py, vx, vy;
    	previous_timestamp_ = measurement_pack.timestamp_;


    }

    if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {

        // Initialize state.
    	float px = measurement_pack.raw_measurements_[0];
    	float py = measurement_pack.raw_measurements_[1];
      	if(fabs(px) < .0001){
        		px = .0001;
        	}

      	else if(fabs(py) < .0001){
        		py = .0001;
        	}

    	ekf_.x_ << px, py, 0, 0;
    	previous_timestamp_ = measurement_pack.timestamp_;



    }

    // done initializing, no need to predict or update

    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /*
   Below are the different thing to be done:

   	     * Update the state transition matrix F according to the new elapsed time.
      	 * Update the process noise covariance matrix.
     	 * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
	float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	 //dt - expressed in seconds
	previous_timestamp_ = measurement_pack.timestamp_;

	// discard predict and update in case if dt is very less

	if (dt <0.001){
		return;
	}

	// Updating the state transition matrix F according to the new elapsed time.

	ekf_.F_ << 1, 0, dt, 0,
			  0, 1, 0, dt,
			  0, 0, 1, 0,
			  0, 0, 0, 1;

	// Set the process covariance matrix Q

	ekf_.Q_ = MatrixXd(4, 4);
	ekf_.Q_ << (pow(dt,4)/4)*noise_ax,0, (pow(dt,3)/2)*noise_ax,0,
	          0,(pow(dt,4)/4)*noise_ay,0,(pow(dt,3)/2)*noise_ay,
	          (pow(dt,3)/2)*noise_ax,0,pow(dt,2)*noise_ax,0,
	          0,(pow(dt,3)/2)*noise_ay,0,pow(dt,2)*noise_ay;


  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /*

     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */
 // kf_.Hj_ = tools.CalculateJacobian(kf_.x_);

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
	  ekf_.R_ = R_radar_;
	  ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
	  ekf_.UpdateEKF(measurement_pack.raw_measurements_);



  }

  if (measurement_pack.sensor_type_ == MeasurementPackage::LASER){
    // Laser updates
	  ekf_.R_ = R_laser_;
	  ekf_.H_ = H_laser_;
	  ekf_.Update(measurement_pack.raw_measurements_);

  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
