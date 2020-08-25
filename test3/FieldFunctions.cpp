#include <eigen3/Eigen/Dense>
#include <cmath>
double E_0 = 0.1;
Eigen::Vector3d E(const Eigen::Vector3d &pos, const Eigen::Vector3d &vel, double t)
{
	Eigen::Vector3d E_field(-E_0*sin(pos[0]),0,0);
	return E_field;
}


Eigen::Vector3d B(const Eigen::Vector3d &pos, const Eigen::Vector3d &vel, double t)
{
	Eigen::Vector3d B_field(0, 0, pow( (1 + pos[1]*pos[1]) , 0.5 ));
	return B_field;
}