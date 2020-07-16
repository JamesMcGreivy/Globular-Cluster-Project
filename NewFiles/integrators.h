#include <eigen3/Eigen/Dense>
#include <cmath>
#include <iostream>


namespace ODE_Integrator
{

	template <typename v>
	double magnitude(const v vector)
	{
		return std::sqrt( vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2] );
	}
		
	
	//The main variables that are Set at the start of each integration
	Eigen::MatrixXd positions;
	Eigen::MatrixXd velocities;
	double t;
	double dt;
	int num_particles;
	Eigen::MatrixXd (*ODE_function)(const Eigen::MatrixXd &positions, const Eigen::MatrixXd &velocities, double t);
	
		
	//Used to return the position and velocity arrays from the main integrator function
	struct ArrayPair 
	{
			
		Eigen::MatrixXd array1;
		Eigen::MatrixXd array2;
		double t;

		ArrayPair(Eigen::MatrixXd &arr1, Eigen::MatrixXd &arr2) 
		{
			array1 = arr1;
			array2 = arr2;
			t = -1000;
		}

		ArrayPair(Eigen::MatrixXd &arr1, Eigen::MatrixXd &arr2, double t_) 
		{
			array1 = arr1;
			array2 = arr2;
			t = t_;
		}
			

		ArrayPair()
		{

		}

		Eigen::MatrixXd &operator[](int i)
		{
			assert(i == 0 or i == 1);
			if (i == 0)
			{
				return array1;
			}
			if (i == 1)
			{
				return array2;
			}
		}
		
	};	

	void leap_frog(Eigen::MatrixXd &positions, Eigen::MatrixXd &velocities, double dt);
	ODE_Integrator::ArrayPair Integrate(Eigen::MatrixXd positions, Eigen::MatrixXd velocities, double dt, Eigen::MatrixXd (*calc_acc)(const Eigen::MatrixXd &positions, const Eigen::MatrixXd &velocities, double t), std::string method);	

}



