#include <eigen3/Eigen/Dense>
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <vector>
#include "NewFiles/integrators.cpp"
#include <chrono>

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

	if (distance < 0.8)
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
				
		//Initializes the acceleration matrix with all zeros
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

//Example (p,q) function representing simple harmonic motion
Eigen::MatrixXd p_dot_SHM(const Eigen::MatrixXd &p, const Eigen::MatrixXd &q, double t)
{
	const int num_particles = p.rows();
	Eigen::MatrixXd p_dot_values(num_particles,3);

	for (int particle = 0; particle < num_particles; particle ++)
	{
		p_dot_values.row(particle) = -1*q.row(particle);
	}

	return p_dot_values;
}

Eigen::MatrixXd q_dot_SHM(const Eigen::MatrixXd &p, const Eigen::MatrixXd &q, double t)
{
	const int num_particles = p.rows();
	Eigen::MatrixXd q_dot_values(num_particles,3);

	for (int particle = 0; particle < num_particles; particle ++)
	{
		for (int component = 0; component < 3; component ++)
		{
			q_dot_values.row(particle) = p.row(particle);
		}
	}

	return q_dot_values;
}


//ODE Representing non-relativistic (E + v x B) force
Eigen::Vector3d E(const Eigen::Vector3d &pos, const Eigen::Vector3d &vel, double t)
{

	Eigen::Vector3d E_field(1-(0.000001)*pos[0],0,0);
	return E_field;

}

Eigen::Vector3d B(const Eigen::Vector3d &pos, const Eigen::Vector3d &vel, double t)
{

	Eigen::Vector3d B_field(0,0,1);
	return B_field;

}



void output_data(std::vector<Data> data)
{

	double num_particles = data[0][0].rows();

	std::ofstream position_data;
	position_data.open("NewFiles/Output/positions.csv");
	
	std::ofstream velocity_data;
	velocity_data.open("NewFiles/Output/velocities.csv");

	for (int i = 0; i < data.size(); i += 1)
	{	
		
		for (int particle = 0; particle < num_particles; particle ++)
		{
			for (int component = 0; component < 3; component ++)
			{
				position_data << std::setprecision(15) << data[i][0](particle,component) << ",";
				velocity_data << std::setprecision(15) << data[i][1](particle,component) << ",";
			}
		}

		position_data << data[i].time << std::endl;
		velocity_data << data[i].time << std::endl;
	}

}


int main(int argc, char *argv[])
{

	Eigen::Matrix<double,1,3> pos;
	Eigen::Matrix<double,1,3> vel;

	pos << 	0,0,0;
	
	vel <<	0,0,0;

	double t_init = 0;
	double dt = (2*M_PI)/10;
	double t_final = 200*dt;

	auto start = std::chrono::high_resolution_clock::now();
	std::vector<Data> data1 = LorentzForceIntegrator::Integrate(E, B, t_init, t_final, pos, vel, dt);
	auto end = std::chrono::high_resolution_clock::now();

	auto duration = duration_cast<std::chrono::milliseconds>(end - start);

	std::cout << "time to run: " << duration.count()/1000.0 << "s" << std::endl;

	output_data(data1);


	return 0;

}