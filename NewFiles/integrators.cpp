#include <eigen3/Eigen/Dense>
#include <cmath>
#include <iostream>
#include <vector>
#include "integrators.h"
#include <fstream>

namespace ODE_Integrator
{

	namespace LeapFrog
	{
		
		void LF_Step(double &dt)
		{
			positions += (dt/2) * velocities;
			velocities += dt * ODE_function(positions, velocities, t);
			positions += (dt/2) * velocities;
			t += dt;
		}

	}

	namespace AdaptiveStep
	{

		double dt_min;
	
		void kick(double dt, Eigen::VectorXi &to_be_kicked)
		{

			Eigen::MatrixXd accelerations = ODE_function(positions, velocities, t);

			for (int i = 0; i < num_particles; i++)
			{
				if (to_be_kicked[i])
				{
					velocities.row(i) += dt * accelerations.row(i); 
				}
			}

		}

		void kick(double dt)
		{
			Eigen::MatrixXd accelerations = ODE_function(positions, velocities, t);

			for (int i = 0; i < num_particles; i++)
			{
				velocities.row(i) += dt * accelerations.row(i); 
			}
		}

		void drift(double dt, Eigen::VectorXi &to_be_drifted)
		{
			for (int i = 0; i < num_particles; i++)
			{
				if (to_be_drifted[i])
				{
					positions.row(i) += dt * velocities.row(i); 
				}
			}				
		}

		void drift(double dt)
		{
			for (int i = 0; i < num_particles; i++)
			{
				positions.row(i) += dt * velocities.row(i); 
			}				
		}

		Eigen::VectorXi select(double dt, Eigen::VectorXi &already_on_timestep)
		{
				
			Eigen::VectorXi on_timestep;
			on_timestep.resize(num_particles,1);

			bool is_on_timestep = false;
			double average_velocity = 0;
			double ideal_dt;
			dt_min = dt;
	
			for (int i = 0; i < num_particles; i ++)
			{
				average_velocity += magnitude( velocities.row(i) );
			}
			average_velocity = average_velocity / num_particles;

			for (int i = 0; i < num_particles; i++)
			{

				ideal_dt = average_velocity / ( 1 + magnitude( velocities.row(i) ) );
					
				is_on_timestep = false;

				if (ideal_dt >= dt)
				{
					is_on_timestep = true;
				}
				else
				{
					dt_min = ideal_dt;
				}

				if ( (is_on_timestep) and (not already_on_timestep[i]) )
				{
					already_on_timestep[i] = 1;
					on_timestep[i] = 1;
				}
				else
				{	
					on_timestep[i] = 0;
				}
			}

			return on_timestep;

		}

		void adaptive_step_recurse(double dt, Eigen::VectorXi already_on_timestep)
		{

			drift(dt/2);

			Eigen::VectorXi on_this_timestep = select(dt, already_on_timestep);

			std::cout << "dt: " << dt << " dt_min: " << dt_min << std::endl;
 				
 			if (dt <= dt_min)
			{
				kick(dt, on_this_timestep);
				drift(dt/2);
			}

			else
			{

				drift(-dt/2);
					
				adaptive_step_recurse(dt/2, already_on_timestep);
					
				kick(dt,on_this_timestep);

				adaptive_step_recurse(dt/2, already_on_timestep);

			}
		}

		void Adaptive_Step(double &dt)
		{
			Eigen::VectorXi already_on_timestep;
			already_on_timestep.resize(num_particles,1);

			for (int i = 0; i < num_particles; i ++)
			{
				already_on_timestep[i] = 0;		
			}

			AdaptiveStep::adaptive_step_recurse(dt,already_on_timestep);
			t += dt;

		}
	}


	namespace RungaKutta
	{

		std::vector<Eigen::Matrix<double,8,3> > b_values;
		std::vector<Eigen::Matrix<double,8,3> > g_values;
			
		Eigen::MatrixXd init_pos;
		Eigen::MatrixXd init_vel;
		Eigen::MatrixXd init_acc;
			
		Eigen::Matrix<double,8,8> c_m;

		static const double h[9]  = { 0.0, 0.0562625605369221464656521910318, 0.180240691736892364987579942780, 0.352624717113169637373907769648, 0.547153626330555383001448554766, 0.734210177215410531523210605558, 0.885320946839095768090359771030, 0.977520613561287501891174488626};
		void get_c_m() {

			c_m <<	1.00000000000000000000, -0.0562625605369221500, 0.01014080283006363000, -0.0035758977292516170, 0.00195656540994722100, -0.0014365302363708915, 0.00127179030902686780, -0.0012432012432012432, 
					0.00000000000000000000, 1.00000000000000000000, -0.2365032522738145200, 0.09353769525946207000, -0.0547553868890686900, 0.04215852772126870600, -0.0387603579159067700, 0.03916083916083916000, 
					0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000, -0.5891279693869842000, 0.41588120008230690000, -0.3600995965020568000, 0.36096224345284600000, -0.3916083916083916000, 
					0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000, -1.1362815957175396000, 1.25015071184069100000, -1.4668842084004270000, 1.79487179487179470000, 
					0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000, -1.8704917729329500000, 2.90613625930842900000, -4.3076923076923075000, 
					0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000, -2.7558127197720457000, 5.60000000000000000000, 
					0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000, -3.7333333333333334000, 
					0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000;

		}


		void substep(double h, double dt) 
		{

			Eigen::Vector3d pos(0,0,0);
			Eigen::Vector3d vel(0,0,0);

			for (int i = 0; i < num_particles; i++) 
			{

				for (int j = 0; j < 3; j++) 
				{
					pos[j] = init_pos(i,j)
							+ ((h*dt) * init_vel(i,j))
							+ ((h*h*dt*dt/2)*init_acc(i,j)) 
							+ ((h*h*h*dt*dt/6)*b_values[i](0,j)) 
							+ ((h*h*h*h*dt*dt/12)*b_values[i](1,j)) 
							+ ((h*h*h*h*h*dt*dt/20)*b_values[i](2,j))
							+ ((h*h*h*h*h*h*dt*dt/30)*b_values[i](3,j))
							+ ((h*h*h*h*h*h*h*dt*dt/42)*b_values[i](4,j))
							+ ((h*h*h*h*h*h*h*h*dt*dt/56)*b_values[i](5,j))
							+ ((h*h*h*h*h*h*h*h*h*dt*dt/72)*b_values[i](6,j))
							+ ((h*h*h*h*h*h*h*h*h*h*dt*dt/90)*b_values[i](7,j));

					vel[j] = init_vel(i,j)
							+ ((h*dt)*init_acc(i,j)) 
							+ ((h*h*dt/2)*b_values[i](0,j)) 
							+ ((h*h*h*dt/3)*b_values[i](1,j)) 
							+ ((h*h*h*h*dt/4)*b_values[i](2,j))
							+ ((h*h*h*h*h*dt/5)*b_values[i](3,j))
							+ ((h*h*h*h*h*h*dt/6)*b_values[i](4,j))
							+ ((h*h*h*h*h*h*h*dt/7)*b_values[i](5,j))
							+ ((h*h*h*h*h*h*h*h*dt/8)*b_values[i](6,j))
							+ ((h*h*h*h*h*h*h*h*h*dt/9)*b_values[i](7,j));
				}

				positions.row(i) = pos;
				velocities.row(i) = vel;

			}
		}

		//Reassigns the b_values based upon the current g_values
		void convert_g_to_b() 
		{		
			
			Eigen::Matrix<double,8,3> new_b;
				
			for(int i = 0; i < num_particles; i ++) {

				for(int j = 0; j < 3; j ++) {
						
					new_b.col(j) = c_m * g_values[i].col(j);
						
				}

				b_values[i] = new_b;
			}
		}

		//Finds the g_values based on the initial conditions and the current b_values
		void find_g_values(double dt) {
			
			Eigen::MatrixXd next_acc;

			//Steps system forward to dt = h[1];
			substep(h[1], dt);
			next_acc = ODE_function(positions, velocities, t);
			//Calculates all g1 values
			for (int i = 0; i < num_particles; i ++) {

				for (int j = 0; j < 3; j ++) {

					g_values[i](0,j) = ( (next_acc(i,j) - init_acc(i,j)) / h[1] );

				}
			}

			//Steps system forward to dt = h[2];
			substep(h[2], dt);

			next_acc = ODE_function(positions, velocities, t);

			//Calculates all g2 values
			for (int i = 0; i < num_particles; i ++) {

				for (int j = 0; j < 3; j ++) {

					g_values[i](1,j) = ( (next_acc(i,j) - init_acc(i,j) - g_values[i](0,j)*h[2] ) / (h[2]*(h[2]-h[1])) );

				}
			}

			//Steps system forward to dt = h[3];
			substep(h[3], dt);

			next_acc = ODE_function(positions, velocities, t);

			//Calculates all g3 values
			for (int i = 0; i < num_particles; i ++) {

				for (int j = 0; j < 3; j ++) {

					g_values[i](2,j) = ( (next_acc(i,j) - init_acc(i,j) - g_values[i](0,j)*h[3] - g_values[i](1,j)*h[3]*(h[3]-h[1]) ) / (h[3]*(h[3]-h[1])*(h[3]-h[2])) );

				}
			}

			//Steps system forward to dt = h[4]
			substep(h[4], dt);

			next_acc = ODE_function(positions, velocities, t);

			//Calculates all g4 values
			for (int i = 0; i < num_particles; i ++) {

				for (int j = 0; j < 3; j ++) {

					g_values[i](3,j) = ( (next_acc(i,j) - init_acc(i,j) - g_values[i](0,j)*h[4] - g_values[i](1,j)*h[4]*(h[4]-h[1]) - g_values[i](2,j)*h[4]*(h[4]-h[1])*(h[4]-h[2]) ) / (h[4]*(h[4]-h[1])*(h[4]-h[2]) * (h[4]-h[3])) );

				}
			}

			//Steps system forward to dt = h[5]
			substep(h[5], dt);

			next_acc = ODE_function(positions, velocities, t);

			//Calculates all g5 values
			for (int i = 0; i < num_particles; i ++) {

				for (int j = 0; j < 3; j ++) {

					g_values[i](4,j) = ( (next_acc(i,j) - init_acc(i,j) - g_values[i](0,j)*h[5] - g_values[i](1,j)*h[5]*(h[5]-h[1]) - g_values[i](2,j)*h[5]*(h[5]-h[1])*(h[5]-h[2]) - g_values[i](3,j)*h[5]*(h[5]-h[1])*(h[5]-h[2])*(h[5]-h[3]) ) / (h[5]*(h[5]-h[1])*(h[5]-h[2])*(h[5]-h[3])*(h[5]-h[4])) );

				}
			}

			//Steps system forward to dt = h[6]
			substep(h[6], dt);

			next_acc = ODE_function(positions, velocities, t);

			//Calculates all g6 values
			for (int i = 0; i < num_particles; i ++) {

				for (int j = 0; j < 3; j ++) {

					g_values[i](5,j) = ( (next_acc(i,j) - init_acc(i,j) - g_values[i](0,j)*h[6] - g_values[i](1,j)*h[6]*(h[6]-h[1]) - g_values[i](2,j)*h[6]*(h[6]-h[1])*(h[6]-h[2]) - g_values[i](3,j)*h[6]*(h[6]-h[1])*(h[6]-h[2])*(h[6]-h[3]) - g_values[i](4,j)*h[6]*(h[6]-h[1])*(h[6]-h[2])*(h[6]-h[3])*(h[6]-h[4]) ) / (h[6]*(h[6]-h[1])*(h[6]-h[2])*(h[6]-h[3])*(h[6]-h[4])*(h[6]-h[5])) );

				}

			}

			//Steps system forward to dt=h[7]
			substep(h[7], dt);

			next_acc = ODE_function(positions, velocities, t);

			//Calculates all g7 values
			for (int i = 0; i < num_particles; i ++) {

				for (int j = 0; j < 3; j ++) {

					g_values[i](6,j) = ( (next_acc(i,j) - init_acc(i,j) - g_values[i](0,j)*h[7] - g_values[i](1,j)*h[7]*(h[7]-h[1]) - g_values[i](2,j)*h[7]*(h[7]-h[1])*(h[7]-h[2]) - g_values[i](3,j)*h[7]*(h[7]-h[1])*(h[7]-h[2])*(h[7]-h[3]) - g_values[i](4,j)*h[7]*(h[7]-h[1])*(h[7]-h[2])*(h[7]-h[3])*(h[7]-h[4]) - g_values[i](5,j)*h[7]*(h[7]-h[1])*(h[7]-h[2])*(h[7]-h[3])*(h[7]-h[4])*(h[7]-h[5]) ) / (h[7]*(h[7]-h[1])*(h[7]-h[2])*(h[7]-h[3])*(h[7]-h[4])*(h[7]-h[5])*(h[7]-h[6])) );

				}

			}

			//Steps system forward to dt = dt
			substep(1, dt);

			next_acc = ODE_function(positions, velocities, t);

			//Calculates all g8 values
			for (int i = 0; i < num_particles; i ++) {

				for (int j = 0; j < 3; j ++) {

					g_values[i](7,j) = ( (next_acc(i,j) - init_acc(i,j) - g_values[i](0,j) - g_values[i](1,j)*(1-h[1]) - g_values[i](2,j)*(1-h[1])*(1-h[2]) - g_values[i](3,j)*(1-h[1])*(1-h[2])*(1-h[3]) - g_values[i](4,j)*(1-h[1])*(1-h[2])*(1-h[3])*(1-h[4]) - g_values[i](5,j)*(1-h[1])*(1-h[2])*(1-h[3])*(1-h[4])*(1-h[5]) - g_values[i](6,j)*(1-h[1])*(1-h[2])*(1-h[3])*(1-h[4])*(1-h[5])*(1-h[6]) ) / ((1-h[1])*(1-h[2])*(1-h[3])*(1-h[4])*(1-h[5])*(1-h[6])*(1-h[7])) );

				}

			}

			//Returns system to dt = 0 and converts the g values to b values
			substep(0, dt);
			convert_g_to_b();

		}

		bool initialized = false;
		void initialize_RK_step(double dt) 
		{

			init_pos = positions;
			init_vel = velocities;
			init_acc = ODE_function(positions,velocities,t);
				
			if (not initialized)
			{

				Eigen::Matrix<double,8,3> coeff;

				coeff <<	0, 0, 0, 
							0, 0, 0, 
							0, 0, 0, 
							0, 0, 0, 
							0, 0, 0, 
							0, 0, 0, 
							0, 0, 0, 
							0, 0, 0;

				for (int i = 0; i < num_particles; i ++)
				{
					b_values.push_back(coeff);
					g_values.push_back(coeff);
				}

				get_c_m();
				initialized = true;

			}

			find_g_values(dt);

		}

		std::vector<Eigen::Matrix<double,8,3> > old_b_values;
		
		int count = 0;
		double epsilon = 0.028;
		double previous_error = 0;
		void RK_Step(double &dt)
		{
			if (count == 0)
			{
				initialize_RK_step(dt);
			}

			old_b_values = b_values;
			find_g_values(dt);

			//Determines whether the g_values have converged to machine precision yet
			double max_del_b8 = 0;
			double max_y_pp = 0;

			for (int i = 0; i < num_particles; i++)
			{
				for (int j = 0; j < 3; j ++)
				{
					max_del_b8 = std::max( abs(max_del_b8) , abs(old_b_values[i](7,j) - b_values[i](7,j) ) );
					max_y_pp = std::max( abs( max_y_pp) , abs( init_acc(i,j) ) );
				}
			}

			double global_error = max_del_b8 / max_y_pp ;


			if ( global_error <= pow(10,-16) or (count > 2 and previous_error <= global_error) or count >= 12 )
			{

				//Calculates if the current step-size is acceptable
				double max_b8 = 0;
				for (int i = 0; i < num_particles; i++)
				{
					for (int j = 0; j < 3; j ++)
					{
						max_b8 = std::max( abs(max_b8) , abs(b_values[i](7,j) ) );
					}
				}
				max_b8 = max_b8 / max_y_pp;

				double dt_req = dt * pow( ( (epsilon * max_y_pp) / max_b8) , 0.14285714285 );
				if (not isnormal(dt_req) or global_error == 0)
				{
					dt_req = 1.5*dt;
				}

				if (dt > dt_req)
				{
					count = 0;
					RK_Step(dt_req);
				}
				
				else
				{
					substep(1, dt);
					std::cout << "converged with count = " << count << " and error = " << global_error << std::endl; 
					std::cout << "dt = " << dt << " and dt_Req = " << dt_req << std::endl;
					count = 0;
					t += dt;
					dt = dt_req;
					
				}
			}

			else
			{
				count += 1;
				previous_error = global_error;
				RK_Step(dt);
			}
		}	
	}

	//fun, t_span, y0, method='RK45', t_eval=None, dense_output=False, events=None, vectorized=False, args=None, **options
	std::vector<ArrayPair> Integrate(Eigen::MatrixXd (*func)(const Eigen::MatrixXd &positions, const Eigen::MatrixXd &velocities, double t), double t_initial, double t_final, const Eigen::MatrixXd &init_positions, const Eigen::MatrixXd &init_velocities, double dt, const Eigen::VectorXd &t_eval = Eigen::VectorXd(), std::string method = "RK")
	{

		//Initializes all of the values
		t = t_initial;
		positions = init_positions;
		velocities = init_velocities;
		std::vector<ArrayPair> data;
		ODE_function = func;
		num_particles = positions.rows();

		//To be used if t_eval is provided by user
		int max_t_eval_index = t_eval.rows();
		bool use_t_eval = max_t_eval_index;
		int t_eval_index = 0;
		double dt_to_t_eval;
		
		//Chooses the integration method used
		void (*integrate_step)(double &dt);
		
		if (method == "RK")
		{
			integrate_step = &RungaKutta::RK_Step;
		}

		if (method == "LF")
		{
			integrate_step = &LeapFrog::LF_Step;
		}

		if (method == "AT")
		{
			integrate_step = &AdaptiveStep::Adaptive_Step;
		}

		while (t < t_final)
		{
			if (use_t_eval and t + dt > t_eval[t_eval_index])
			{	
				dt = t_eval[t_eval_index] - t;
			}
			if (t + dt > t_final)
			{
				dt = t_final - t;
			}

			integrate_step(dt);
			std::cout << t << std::endl;

			if (t == t_eval[t_eval_index])
			{
				data.push_back( ArrayPair(positions, velocities, t));
				t_eval_index += 1;
			}
			if(t_eval_index == max_t_eval_index)
			{
				use_t_eval = false;
			}

		}
		
		return data;

	}
	
}


