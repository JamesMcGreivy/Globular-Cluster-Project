
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