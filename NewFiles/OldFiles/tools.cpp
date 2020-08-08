#include <eigen3/Eigen/Dense>
#include <cmath>
#include <iostream>
#include <vector>

namespace ODE_Integrator
{

	namespace Tools
	{
		//Only works for gravity
		
		double G = 1;
		
		double potential_energy(Eigen::Vector3d pos1, Eigen::Vector3d pos2, double mass1, double mass2)
		{

			double distance = magnitude ( (pos1 - pos2) );
			return ( -G * mass1 * mass2 / distance );

		}

		double total_PE(const Eigen::MatrixXd &positions, Eigen::VectorXd &masses)
		{
			int num_particles = positions.rows();
			double PE = 0;

			for (int p1 = 0; p1 < num_particles; p1 ++)
			{
				
				for (int p2 = 0; p2 < p1; p2++)
				{

					PE += potential_energy( positions.row(p1), positions.row(p2), masses[p1], masses[p2]);

				}
			
			}

			return PE;
		
		}


		double total_KE(const Eigen::MatrixXd &velocities, const Eigen::VectorXd &masses)
		{
			int num_particles = velocities.rows();
			double KE = 0;

			for (int particle = 0; particle < num_particles; particle ++)
			{
				KE += (0.5) * masses[particle] * pow( magnitude( velocities.row(particle) ) , 2 );
			}

			return KE;
		
		}




	}
}