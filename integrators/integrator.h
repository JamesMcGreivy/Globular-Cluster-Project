#include <eigen3/Eigen/Dense>
#include <cmath>
#include <iostream>

	//Used to return the position and velocity arrays from the main integrator function
	struct Data 
	{
			
		Eigen::Vector3d position;
		Eigen::Vector3d velocity;
		double time;

		Data(Eigen::Vector3d &pos, Eigen::Vector3d &vel, double t) 
		{
			position = pos;
			velocity = vel;
			time = t;
		}
			

		Data()
		{

		}

		Eigen::Vector3d &operator[](int i)
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
