#include <eigen3/Eigen/Dense>
#define _USE_MATH_DEFINES
#include <cmath>

Eigen::Vector3d E(const Eigen::Vector3d &pos, const Eigen::Vector3d &vel, double t)
{

    Eigen::Vector3d E_field(0,0,0);
    return E_field;

}

Eigen::Vector3d B(const Eigen::Vector3d &pos, const Eigen::Vector3d &vel, double t)
{

    Eigen::Vector3d B_field(0,0,1);
    return B_field;
}

Eigen::Vector3d posInit(1, 0, 0);
Eigen::Vector3d velInit(0, -1, 0);

double tInit = 0;
double tFinal = 628.318530718;

bool use_adaptive_timestep = true;