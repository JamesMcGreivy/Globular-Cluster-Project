#include <eigen3/Eigen/Dense>
#include <cmath>

Eigen::Vector3d E(const Eigen::Vector3d &pos, const Eigen::Vector3d &vel, double t)
{
	Eigen::Vector3d E_field(0,1,0);
	return E_field;
}


Eigen::Vector3d B(const Eigen::Vector3d &pos, const Eigen::Vector3d &vel, double t)
{
	Eigen::Vector3d B_field(0,0,1);
	return B_field;
}