#include <eigen3/Eigen/Dense>
#include <cmath>
#include <iostream>
#include <vector>
#include "NewFiles/integrators.cpp"
#include "NewFiles/tools.cpp"

//Example acceleration function for gravity
// ( G * m2 / |r|^3 ) * r
double G = 1;
	
Eigen::VectorXd masses;

Eigen::VectorXd calc_gravity(int index_of_1, const Eigen::VectorXd pos1, const Eigen::VectorXd vel1, int index_of_2, const Eigen::VectorXd pos2, const Eigen::VectorXd vel2)
{
	Eigen::VectorXd r = pos2 - pos1;
	double distance = 0;
	for (int comp = 0; comp < pos1.rows(); comp++)
	{
		distance += pow(pos2[comp]-pos1[comp],2);
	}

	distance = pow(distance, 0.5); 

	if (distance < 0.05)
	{
		return 0 * r;
	}
				
	return ( (G * masses[index_of_2] / pow(distance,3) ) * r );
}

Eigen::MatrixXd gravity(const Eigen::MatrixXd &positions, const Eigen::MatrixXd &velocities, double t)
{
	
	assert (positions.cols() == velocities.cols() and positions.rows() == velocities.rows());
	const int num_particles = positions.rows();
	const int num_components = positions.cols();
			
	Eigen::MatrixXd accelerations(num_particles,num_components);

	//Loops through each particle and calculates the net force on that particle
	for (int particle = 0; particle < num_particles; particle ++)
	{
				
	//Initializes the acceleration matrix with all zeroe
	for (int component = 0; component < num_components; component++)
		{
			accelerations(particle,component) = 0; 
		}

		//Loops through all particles besides the current one
		for (int other_particle = 0; other_particle < num_particles; other_particle ++)
		{
			if (particle != other_particle)
			{	

			accelerations.row(particle) += calc_gravity( particle, positions.row(particle), velocities.row(particle), other_particle, positions.row(other_particle), velocities.row(other_particle) );
						
			}
		}
	}

	return accelerations;

}

int main()
{

	Eigen::Matrix<double,2,3> positions;
	Eigen::Matrix<double,2,3> velocities;

	positions <<	-3,0,0,
					3,0,0;


	velocities <<	0,0.3,0,
					0,-.3,0;


	masses = Eigen::Matrix<double,2,1>();
	masses << 1, 1;

	double num_particles = positions.rows();
	double t_final = 120;
	double dt = 20;


	//Creates t_eval with evaluations at each time-step from t_initial to t_final
	Eigen::VectorXd t_eval(int(t_final/dt));
	for (int i = 0; i < (t_final/(dt)); i ++)
	{
		t_eval[i] = (i+1)*(dt);
	}

	//Integrates the function
	std::vector<ODE_Integrator::ArrayPair> data = ODE_Integrator::Integrate(gravity, 0, t_final, positions, velocities, dt, t_eval, "RK");

	//Outputs data to a csv file
	std::ofstream position_data;
	std::ofstream velocity_data;
	position_data.open("NewFiles/Output/positions.csv");
	velocity_data.open("NewFiles/Output/velocities.csv");

	for (int i = 0; i < data.size(); i += 1)
	{	
		for (int particle = 0; particle < num_particles; particle ++)
		{
			for (int component = 0; component < 3; component ++)
			{
				position_data << data[i][0](particle,component) << ",";
				velocity_data << data[i][1](particle,component) << ",";
			}
		}

		position_data << data[i].t << std::endl;
		velocity_data << data[i].t << std::endl;

	}


	//Outputs energy data to a csv file
	std::ofstream energy_data;
	energy_data.open("NewFiles/Output/energy.csv");
	for (int i = 0; i < data.size(); i += 1)
	{
		double total_KE = ODE_Integrator::Tools::total_KE(data[i][1], masses);
		double total_PE = ODE_Integrator::Tools::total_PE(data[i][0], masses);
		double total_energy = total_KE + total_PE;
		double initial_energy;
		if (i == 0)
		{
			initial_energy = total_energy;
		}

		energy_data << "KE: , " 	<< total_KE << " , "
					<< "PE: , " 	<< total_PE << " , "
					<< "total: , " 	<< total_energy << " , "
					<< "delta E: , "<< (total_energy - initial_energy) / initial_energy << " , "
					<< "t: , " 		<< t_eval[i] << std::endl;

	}

	return 0;

}