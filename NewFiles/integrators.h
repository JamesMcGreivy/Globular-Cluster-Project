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
		
	//Used to return the position and velocity arrays from the main integrator function
	struct Data 
	{
			
		Eigen::MatrixXd position;
		Eigen::MatrixXd velocity;
		double time;

		Data(Eigen::MatrixXd &pos, Eigen::MatrixXd &vel, double t) 
		{
			position = pos;
			velocity = vel;
			time = t;
		}
			

		Data()
		{

		}

		Eigen::MatrixXd &operator[](int i)
		{
			assert(i == 0 or i == 1);
			if (i == 0)
			{
				return position;
			}
			if (i == 1)
			{
				return velocity;
			}
		}
		
	};	
	

}



