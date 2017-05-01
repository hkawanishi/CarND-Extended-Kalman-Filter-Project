#include <iostream>
#include "tools.h"

using namespace std;

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
	rmse << 0,0,0,0;
	VectorXd diff;
  // check to make sure:
  // - the estimation vector size should not be zero
  // - the estimation vector size should equal to ground truth.
  if(estimations.size() != ground_truth.size() || estimations.size() == 0){
  	cout << "CalculateRMSE() Error - Check the estimation or ground truth data" << endl;
  	return rmse;
  }
  // calculate the squared residuals.
  for(unsigned int i=0; i < estimations.size(); ++i){
  	diff = estimations[i] - ground_truth[i];
  	diff = diff.array()*diff.array();
  	rmse += diff;
  }
  // calculate the mean
  rmse = rmse/estimations.size();

  // calculate the squared root.
  rmse = rmse.array().sqrt();

  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  MatrixXd Hj(3,4);
  // state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  // px_py_2 is px^2 + py^2
  float px_py_2;
  px_py_2 = pow(px,2) + pow(py,2);

  // if the denominator is zero or very small, give an error message
  if (fabs(px_py_2) < 0.0001){
  	cout << "CalculateJacobian() - Error - Divided by zero" << endl;
    return Hj;
  } else {
		Hj << px/(sqrt(px_py_2)),py/(sqrt(px_py_2)), 0, 0,
	    -py/(px_py_2), px/(px_py_2), 0, 0,
	    py*(vx*py-vy*px)/pow(px_py_2, 1.5), px*(vy*px-vx*py)/pow(px_py_2, 1.5), px/sqrt(px_py_2), py/sqrt(px_py_2);
		return Hj;
  }
}
