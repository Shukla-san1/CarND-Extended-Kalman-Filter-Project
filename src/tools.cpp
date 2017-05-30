#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {

}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**

    * This method is to calculate the RMSE here.
  */
	VectorXd rmse(4);
	rmse << 0,0,0,0;
	if(estimations.size() != ground_truth.size() || estimations.size() == 0  )
	{
	    cout <<"Error: invalid input" <<endl;
	    return rmse;
	}

	//accumulate squared error
	for(int i=0; i < estimations.size(); ++i){

        VectorXd e = estimations[i] - ground_truth[i];
        e = e.array()*e.array();
        rmse +=e;

	}

	//calculate the mean

	rmse = rmse/estimations.size();

	//calculate the squared root

	rmse = rmse.array().sqrt();

	//return the result
	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**

    * This method is to calculate a Jacobian
  */
	MatrixXd Hj(3,4);

	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);


	//pre-compute a set of terms to avoid repeated calculation

	float c1 = px*px+py*py;
	float c2 = sqrt(c1);
	float c3 = (c1*c2);

	//checking division by zero

	if(fabs(c1) < 0.0001){
		cout << "CalculateJacobian () - Error - Division by Zero" << endl;
		return Hj;
	}
	//compute the Jacobian matrix

	Hj << (px/c2), (py/c2), 0, 0,
		  -(py/c1), (px/c1), 0, 0,
		  py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;

	return Hj;

}
