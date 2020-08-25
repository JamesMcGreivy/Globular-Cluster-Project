#define _USE_MATH_DEFINES
#include <eigen3/Eigen/Dense>
#include <cmath>
#include <iostream>
#include <vector>
#include "integrator.h"
#include <fstream>

namespace LorentzForceIntegrator
{

	template <typename v>
	double magnitude(const v vector)
	{
		return std::sqrt( vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2] );
	}

	class RungaKuttaIntegrator
	{
	public:
		bool use_predicted_b;
		bool use_adaptive_timestep;
		double adaptive_factor;

		double t = 0;
		double dt;
		int num_particles;
		Eigen::MatrixXd positions;
		Eigen::MatrixXd velocities;

		std::vector<Eigen::Matrix<double,7,3> > b_values;
		std::vector<Eigen::Matrix<double,7,3> > g_values;
		Eigen::MatrixXd init_pos;
		Eigen::MatrixXd init_vel;
		Eigen::MatrixXd init_acc;
		double init_t;

		int count = 0;
		
		Eigen::Matrix<double,7,7> c_m;
		const double h[8] = { 0.0, 0.0562625605369221464656521910318, 0.180240691736892364987579942780, 0.352624717113169637373907769648, 0.547153626330555383001448554766, 0.734210177215410531523210605558, 0.885320946839095768090359771030, 0.977520613561287501891174488626};


		Eigen::Vector3d (*E)(const Eigen::Vector3d &pos, const Eigen::Vector3d &vel, double t);
		Eigen::Vector3d (*B)(const Eigen::Vector3d &pos, const Eigen::Vector3d &vel, double t);

		RungaKuttaIntegrator(const Eigen::MatrixXd &init_pos, const Eigen::MatrixXd &init_vel, Eigen::Vector3d (*E_function)(const Eigen::Vector3d &pos, const Eigen::Vector3d &vel, double t), Eigen::Vector3d (*B_function)(const Eigen::Vector3d &pos, const Eigen::Vector3d &vel, double t), double init_t, double init_dt, bool predict_b = true, bool adaptive_timestep = true, double adaptive_period_factor = 0.125)
		{ 
			c_m <<	1.00000000000000000000, -0.0562625605369221500, 0.01014080283006363000, -0.0035758977292516170, 0.00195656540994722100, -0.0014365302363708915, 0.00127179030902686780, 
					0.00000000000000000000, 1.00000000000000000000, -0.2365032522738145200, 0.09353769525946207000, -0.0547553868890686900, 0.04215852772126870600, -0.0387603579159067700, 
					0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000, -0.5891279693869842000, 0.41588120008230690000, -0.3600995965020568000, 0.36096224345284600000, 
					0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000, -1.1362815957175396000, 1.25015071184069100000, -1.4668842084004270000, 
					0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000, -1.8704917729329500000, 2.90613625930842900000, 
					0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000, -2.7558127197720457000, 
					0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000;
			
			use_predicted_b = predict_b;
			use_adaptive_timestep = adaptive_timestep;
			adaptive_factor = adaptive_period_factor;

			E = E_function;
			B = B_function;

			positions = init_pos;
			velocities = init_vel;
			num_particles = positions.rows();
			t = init_t;
			dt = init_dt;

		}

	private:

		Eigen::MatrixXd ODE_function(const Eigen::MatrixXd &pos, const Eigen::MatrixXd &vel, double t)
		{
			Eigen::MatrixXd acceleration(num_particles,3);
			
			Eigen::Vector3d acc(0,0,0);
			Eigen::Vector3d v;
			
			for (int p = 0; p < num_particles; p++)
			{
				v = vel.row(p);
				
				acc = E( pos.row(p), v, t) + v.cross( B(pos.row(p),v,t) );
				acceleration.row(p) = acc;

			}

			return acceleration;

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
						+ (h*h*dt*dt/2)		*	(init_acc(i,j)
						+ (h/3)				*	(b_values[i](0,j) 
						+ (h/2)				*	(b_values[i](1,j) 
						+ (3*h/5)			*	(b_values[i](2,j)
						+ (4*h/6)			*	(b_values[i](3,j)
						+ (5*h/7)			*	(b_values[i](4,j)
						+ (6*h/8)			*	(b_values[i](5,j)
						+ (7*h/9)			*	(b_values[i](6,j)))))))));

					vel[j] = init_vel(i,j)
						+ (h*dt)	*	(init_acc(i,j)
						+ (h/2)		*	(b_values[i](0,j) 
						+ (2*h/3)	*	(b_values[i](1,j) 
						+ (3*h/4)	*	(b_values[i](2,j)
						+ (4*h/5)	*	(b_values[i](3,j)
						+ (5*h/6)	*	(b_values[i](4,j)
						+ (6*h/7)	*	(b_values[i](5,j)
						+ (7*h/8)	*	(b_values[i](6,j)))))))));
				}

				positions.row(i) = pos;
				velocities.row(i) = vel;

			}

			t = init_t + dt*h;

		}

		//Reassigns the b_values based upon the current g_values
		void convert_g_to_b() 
		{		
			Eigen::Matrix<double,7,3> new_b;

			for(int i = 0; i < num_particles; i ++) 
			{
				for(int j = 0; j < 3; j ++) 
				{					
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

					g_values[i](1,j) = ( (next_acc(i,j) - init_acc(i,j) - h[2]*g_values[i](0,j)) / (h[2]*(h[2]-h[1])) );

				}
			}

			//Steps system forward to dt = h[3];
			substep(h[3], dt);

			next_acc = ODE_function(positions, velocities, t);

			//Calculates all g3 values
			for (int i = 0; i < num_particles; i ++) {

				for (int j = 0; j < 3; j ++) {

					g_values[i](2,j) = ( (next_acc(i,j) - init_acc(i,j) - h[3]*(g_values[i](0,j) + (h[3]-h[1])*g_values[i](1,j))) / (h[3]*(h[3]-h[1])*(h[3]-h[2])) );

				}
			}

			//Steps system forward to dt = h[4]
			substep(h[4], dt);

			next_acc = ODE_function(positions, velocities, t);

			//Calculates all g4 values
			for (int i = 0; i < num_particles; i ++) {

				for (int j = 0; j < 3; j ++) {

					g_values[i](3,j) = ( (next_acc(i,j) - init_acc(i,j) - h[4]*(g_values[i](0,j) + (h[4]-h[1])*(g_values[i](1,j) + (h[4]-h[2])*(g_values[i](2,j))))) / (h[4]*(h[4]-h[1])*(h[4]-h[2])*(h[4]-h[3])) );

				}
			}

			//Steps system forward to dt = h[5]
			substep(h[5], dt);

			next_acc = ODE_function(positions, velocities, t);

			//Calculates all g5 values
			for (int i = 0; i < num_particles; i ++) {

				for (int j = 0; j < 3; j ++) {

					g_values[i](4,j) = ( (next_acc(i,j) - init_acc(i,j) - h[5]*(g_values[i](0,j) + (h[5]-h[1])*(g_values[i](1,j) + (h[5]-h[2])*(g_values[i](2,j) + (h[5]-h[3])*(g_values[i](3,j)))))) / (h[5]*(h[5]-h[1])*(h[5]-h[2])*(h[5]-h[3])*(h[5]-h[4])) );

				}
			}

			//Steps system forward to dt = h[6]
			substep(h[6], dt);

			next_acc = ODE_function(positions, velocities, t);

			//Calculates all g6 values
			for (int i = 0; i < num_particles; i ++) {

				for (int j = 0; j < 3; j ++) {

					g_values[i](5,j) = ( (next_acc(i,j) - init_acc(i,j) - h[6]*(g_values[i](0,j) + (h[6]-h[1])*(g_values[i](1,j) + (h[6]-h[2])*(g_values[i](2,j) + (h[6]-h[3])*(g_values[i](3,j) + (h[6]-h[4])*(g_values[i](4,j))))))) / (h[6]*(h[6]-h[1])*(h[6]-h[2])*(h[6]-h[3])*(h[6]-h[4])*(h[6]-h[5])) );

				}

			}

			//Steps system forward to dt=h[7]
			substep(h[7], dt);

			next_acc = ODE_function(positions, velocities, t);

			//Calculates all g7 values
			for (int i = 0; i < num_particles; i ++) {

				for (int j = 0; j < 3; j ++) {

					g_values[i](6,j) = ( (next_acc(i,j) - init_acc(i,j) - h[7]*(g_values[i](0,j) + (h[7]-h[1])*(g_values[i](1,j) + (h[7]-h[2])*(g_values[i](2,j) + (h[7]-h[3])*(g_values[i](3,j) + (h[7]-h[4])*(g_values[i](4,j) + (h[7]-h[5])*(g_values[i](5,j)))))))) / (h[7]*(h[7]-h[1])*(h[7]-h[2])*(h[7]-h[3])*(h[7]-h[4])*(h[7]-h[5])*(h[7]-h[6])) );

				}

			}

			//Returns system to dt = 0 and converts the g values to b values
			substep(0, dt);
			convert_g_to_b();

		}


		void init_substep(double h)
		{

			Eigen::Vector3d pos(0,0,0);
			Eigen::Vector3d vel(0,0,0);

			Eigen::Vector3d C1(0,0,0);
			Eigen::Vector3d C2(0,0,0);
			
			Eigen::Vector3d V_par(0,0,0);
			Eigen::Vector3d E_par(0,0,0);
			
			Eigen::Vector3d const_E(0,0,0);
			Eigen::Vector3d const_B(0,0,0);
			
			Eigen::Vector3d V_0(0,0,0);
			Eigen::Vector3d X_0(0,0,0);
			double Omega;

			for (int i = 0; i < num_particles; i++) 
			{
				V_0 = init_vel.row(i);
				X_0 = init_pos.row(i);

				const_E = E(X_0, V_0, t);
				const_B = B(X_0, V_0, t);
				Omega = magnitude(const_B);

				if (Omega == 0)
				{
					pos = (const_E*(dt*h)*(dt*h)/2) + V_0*(dt*h) + X_0;
					vel = const_E*(dt*h) + V_0;
				}
				else
				{	

					E_par = (const_E.dot(const_B)/(Omega*Omega))*const_B;
					V_par = (V_0.dot(const_B)/(Omega*Omega))*const_B;

					C1 = ( V_0 - ( ((V_0.dot(const_B) * const_B) + const_E.cross(const_B)) / (Omega*Omega)) );
					C2 = (1/Omega) * (const_E - ( (const_E.dot(const_B)/(Omega*Omega))*const_B ) + V_0.cross(const_B));

					pos = 	E_par*(dt*h)*(dt*h)/2 + V_par*(dt*h)
							+ (C1/Omega)*std::sin(Omega*dt*h) - (C2/Omega)*std::cos(Omega*dt*h) 
							+ ((const_E.cross(const_B)/(Omega*Omega))*dt*h) + (C2/Omega) + X_0;

					vel = 	E_par*dt*h + V_par 
							+ C1*std::cos(Omega*dt*h) + C2*std::sin(Omega*dt*h) 
							+ (const_E.cross(const_B)/(Omega*Omega));
				}

				positions.row(i) = pos;
				velocities.row(i) = vel;

			}

			t = init_t + dt*h;

		}

		//Finds the g_values based on the initial conditions and the current b_values
		void find_initial_g_values(double dt) {

			Eigen::MatrixXd next_acc;

			//Steps system forward to dt = h[1];
			init_substep(h[1]);
			
			next_acc = ODE_function(positions, velocities, t);
			//Calculates all g1 values
			for (int i = 0; i < num_particles; i ++) {

				for (int j = 0; j < 3; j ++) {

					g_values[i](0,j) = ( (next_acc(i,j) - init_acc(i,j)) / h[1] );
		
				}
			}

			//Steps system forward to dt = h[2];
			init_substep(h[2]);

			next_acc = ODE_function(positions, velocities, t);

			//Calculates all g2 values
			for (int i = 0; i < num_particles; i ++) {

				for (int j = 0; j < 3; j ++) {

					g_values[i](1,j) = ( (next_acc(i,j) - init_acc(i,j) - h[2]*g_values[i](0,j)) / (h[2]*(h[2]-h[1])) );

				}
			}

			//Steps system forward to dt = h[3];
			init_substep(h[3]);

			next_acc = ODE_function(positions, velocities, t);

			//Calculates all g3 values
			for (int i = 0; i < num_particles; i ++) {

				for (int j = 0; j < 3; j ++) {

					g_values[i](2,j) = ( (next_acc(i,j) - init_acc(i,j) - h[3]*(g_values[i](0,j) + (h[3]-h[1])*g_values[i](1,j))) / (h[3]*(h[3]-h[1])*(h[3]-h[2])) );

				}
			}

			//Steps system forward to dt = h[4]
			init_substep(h[4]);

			next_acc = ODE_function(positions, velocities, t);

			//Calculates all g4 values
			for (int i = 0; i < num_particles; i ++) {

				for (int j = 0; j < 3; j ++) {

					g_values[i](3,j) = ( (next_acc(i,j) - init_acc(i,j) - h[4]*(g_values[i](0,j) + (h[4]-h[1])*(g_values[i](1,j) + (h[4]-h[2])*(g_values[i](2,j))))) / (h[4]*(h[4]-h[1])*(h[4]-h[2])*(h[4]-h[3])) );

				}
			}

			//Steps system forward to dt = h[5]
			init_substep(h[5]);

			next_acc = ODE_function(positions, velocities, t);

			//Calculates all g5 values
			for (int i = 0; i < num_particles; i ++) {

				for (int j = 0; j < 3; j ++) {

					g_values[i](4,j) = ( (next_acc(i,j) - init_acc(i,j) - h[5]*(g_values[i](0,j) + (h[5]-h[1])*(g_values[i](1,j) + (h[5]-h[2])*(g_values[i](2,j) + (h[5]-h[3])*(g_values[i](3,j)))))) / (h[5]*(h[5]-h[1])*(h[5]-h[2])*(h[5]-h[3])*(h[5]-h[4])) );

				}
			}

			//Steps system forward to dt = h[6]
			init_substep(h[6]);

			next_acc = ODE_function(positions, velocities, t);

			//Calculates all g6 values
			for (int i = 0; i < num_particles; i ++) {

				for (int j = 0; j < 3; j ++) {

					g_values[i](5,j) = ( (next_acc(i,j) - init_acc(i,j) - h[6]*(g_values[i](0,j) + (h[6]-h[1])*(g_values[i](1,j) + (h[6]-h[2])*(g_values[i](2,j) + (h[6]-h[3])*(g_values[i](3,j) + (h[6]-h[4])*(g_values[i](4,j))))))) / (h[6]*(h[6]-h[1])*(h[6]-h[2])*(h[6]-h[3])*(h[6]-h[4])*(h[6]-h[5])) );

				}

			}

			//Steps system forward to dt=h[7]
			init_substep(h[7]);

			next_acc = ODE_function(positions, velocities, t);

			//Calculates all g7 values
			for (int i = 0; i < num_particles; i ++) {

				for (int j = 0; j < 3; j ++) {

					g_values[i](6,j) = ( (next_acc(i,j) - init_acc(i,j) - h[7]*(g_values[i](0,j) + (h[7]-h[1])*(g_values[i](1,j) + (h[7]-h[2])*(g_values[i](2,j) + (h[7]-h[3])*(g_values[i](3,j) + (h[7]-h[4])*(g_values[i](4,j) + (h[7]-h[5])*(g_values[i](5,j)))))))) / (h[7]*(h[7]-h[1])*(h[7]-h[2])*(h[7]-h[3])*(h[7]-h[4])*(h[7]-h[5])*(h[7]-h[6])) );

				}

			}

			//Returns system to dt = 0 and converts the g values to b values
			substep(0, dt);

			convert_g_to_b();
		
		}

		Eigen::Matrix<double,7,3> init_b_values;

		bool initialized = false;
		void initialize_step() 
		{

			init_pos = positions;
			init_vel = velocities;
			init_acc = ODE_function(positions,velocities,t);
			init_t = t;
				
			if (not initialized)
			{

				Eigen::Matrix<double,7,3> coeff;

				coeff <<	0, 0, 0, 
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

				initialized = true;

			}

			if (use_adaptive_timestep)
			{
				init_substep(1);
				double Omega = std::max( 
								magnitude(B(init_pos.row(0),init_vel.row(0),init_t)),
								magnitude(B(positions.row(0),velocities.row(0),t))
								);
				if (Omega != 0)
				{
					dt = (1/Omega) * (2*M_PI) * (adaptive_factor);
				}
			}

			if (use_predicted_b)
			{
				find_initial_g_values(dt);
			}
			else
			{
				find_g_values(dt);
			}

		}

	public:
		
		std::vector<Eigen::Matrix<double,7,3> > old_b_values;
		
		double max_del_b8;
		double max_y_pp;
		double global_error;
		double total_count = 0;

		void step()
		{
			if (count == 0)
			{
				initialize_step();
			}
			
			old_b_values = b_values;
			find_g_values(dt);

			//Determines whether the g_values have converged to machine precision yet
			max_del_b8 = 0;
			max_y_pp = 0;

			for (int i = 0; i < num_particles; i++)
			{
				for (int j = 0; j < 3; j ++)
				{
					max_del_b8 = std::max( abs(max_del_b8) , abs(old_b_values[i](6,j) - b_values[i](6,j) ) );
					max_y_pp = std::max( abs( max_y_pp) , abs( init_acc(i,j) ) );
				}
			}
			
			global_error = max_del_b8 / max_y_pp;
			if (count > 4 and (global_error <= 1e-15 or count >= 20))
			{	
				substep(1, dt);
				std::cout << t << "  " << count << std::endl;
				total_count += count;
				count = 0;

			}
			else
			{
				count += 1;
				step();
			}
		}
	};

	std::vector<Data> Integrate(Eigen::Vector3d (*E_function)(const Eigen::Vector3d &pos, const Eigen::Vector3d &vel, double t), Eigen::Vector3d (*B_function)(const Eigen::Vector3d &pos, const Eigen::Vector3d &vel, double t), double t_initial, double t_final, const Eigen::MatrixXd &init_pos, const Eigen::MatrixXd &init_vel, double dt, bool use_predicted_b = true, bool use_adaptive_timestep = true, double adaptive_factor = 0.125)
	{

		std::vector<Data> return_data;
		RungaKuttaIntegrator integrator = RungaKuttaIntegrator(init_pos, init_vel, E_function, B_function, t_initial, dt, use_predicted_b, use_adaptive_timestep, adaptive_factor);
		return_data.push_back( Data(integrator.positions, integrator.velocities, integrator.t) );

		while (integrator.t < t_final)
		{
			
			integrator.step();


			if (integrator.t <= t_final)
			{
				return_data.push_back( Data(integrator.positions, integrator.velocities, integrator.t) );	
			}
		
		}

		std::cout << integrator.total_count << std::endl;

		return return_data;

	}	

}

namespace BorisIntegrator
{

	template <typename v>
	double magnitude(const v vector)
	{
		return std::sqrt( vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2] );
	}

	class BorisIntegrator
	{
	public:
		double t = 0;
		double dt;
		
		int num_particles;
		Eigen::MatrixXd positions;
		Eigen::MatrixXd velocities;

		bool use_adaptive_timestep;
		double adaptive_factor;

		Eigen::Vector3d (*E)(const Eigen::Vector3d &pos, const Eigen::Vector3d &vel, double t);
		Eigen::Vector3d (*B)(const Eigen::Vector3d &pos, const Eigen::Vector3d &vel, double t);

		BorisIntegrator(const Eigen::MatrixXd &init_pos, const Eigen::MatrixXd &init_vel, Eigen::Vector3d (*E_function)(const Eigen::Vector3d &pos, const Eigen::Vector3d &vel, double t), Eigen::Vector3d (*B_function)(const Eigen::Vector3d &pos, const Eigen::Vector3d &vel, double t), double init_t, double init_dt, bool adaptive_timestep = true, double adaptive_period_factor = 0.125)
		{

			positions = init_pos;
			velocities = init_vel;

			E = E_function;
			B = B_function;

			t = init_t;
			dt = init_dt;

			num_particles = positions.rows();

			use_adaptive_timestep = adaptive_timestep;
			adaptive_factor = adaptive_period_factor;

		}

		void step()
		{
			if (use_adaptive_timestep)
			{
				double Omega = magnitude(B(positions.row(0),velocities.row(0),t));
				if (Omega != 0)
				{
					dt = (1/Omega) * (2*M_PI) * (adaptive_factor);
				}
			} 

			double qmdt = 0.5*dt;
			double dt2 = 0.5*dt;

			Eigen::Vector3d pos(0,0,0);
			Eigen::Vector3d vel(0,0,0);
			
			Eigen::Vector3d E_field(0,0,0);
			Eigen::Vector3d B_field(0,0,0);
			
			Eigen::Vector3d vm(0,0,0);
			Eigen::Vector3d t_(0,0,0);
			Eigen::Vector3d s(0,0,0);
			Eigen::Vector3d cp(0,0,0);
			Eigen::Vector3d vprime(0,0,0);
			Eigen::Vector3d vp(0,0,0);

			

			for (int i = 0; i < num_particles; i++)
			{

				pos = positions.row(i);
				vel = velocities.row(i);

				double qmdt = 0.5*dt;
				double dt2 = 0.5*dt;

			    pos[0] = pos[0] + dt2*vel[0];
			    pos[1] = pos[1] + dt2*vel[1];
			    pos[2] = pos[2] + dt2*vel[2];  

			   	E_field = E(pos,vel,t);
			   	double E_0 = E_field[0];
			   	double E_1 = E_field[1];
				double E_2 = E_field[2];
			    
			    B_field = B(pos,vel,t);
			   	double B_0 = B_field[0];
			   	double B_1 = B_field[1];
				double B_2 = B_field[2];

			    double vm_0 = vel[0] + qmdt*E_0;
			    double vm_1 = vel[1] + qmdt*E_1;
			    double vm_2 = vel[2] + qmdt*E_2;

			    double t_0 = qmdt*B_0;
			    double t_1 = qmdt*B_1;
			    double t_2 = qmdt*B_2;

			    double tNorm = t_0*t_0 + t_1*t_1 + t_2*t_2;
			    double s_0 = 2*t_0/(1+tNorm);
			    double s_1 = 2*t_1/(1+tNorm);
			    double s_2 = 2*t_2/(1+tNorm);

			    double cp_0 =  vm_1*t_2 - vm_2*t_1;
			    double cp_1 =  vm_2*t_0 - vm_0*t_2;
			    double cp_2 =  vm_0*t_1 - vm_1*t_0;

			    double vprime_0 = vm_0 + cp_0;
			    double vprime_1 = vm_1 + cp_1;
			    double vprime_2 = vm_2 + cp_2;
			      
			    cp_0 =  vprime_1*s_2 - vprime_2*s_1;
			    cp_1 =  vprime_2*s_0 - vprime_0*s_2;
			    cp_2 =  vprime_0*s_1 - vprime_1*s_0;

			    double vp_0 = vm_0 + cp_0;
			    double vp_1 = vm_1 + cp_1;
			    double vp_2 = vm_2 + cp_2;

			    vel[0] = vp_0 + qmdt*E_0;
			    vel[1] = vp_1 + qmdt*E_1;
			    vel[2] = vp_2 + qmdt*E_2;

			    pos[0] = pos[0] + dt2*vel[0];
			    pos[1] = pos[1] + dt2*vel[1];
			    pos[2] = pos[2] + dt2*vel[2];

			    positions.row(i) = pos;
			    velocities.row(i) = vel;

			}

			t += dt;

		}

	};

	std::vector<Data> Integrate(Eigen::Vector3d (*E_function)(const Eigen::Vector3d &pos, const Eigen::Vector3d &vel, double t), Eigen::Vector3d (*B_function)(const Eigen::Vector3d &pos, const Eigen::Vector3d &vel, double t), double t_initial, double t_final, const Eigen::MatrixXd &init_pos, const Eigen::MatrixXd &init_vel, double dt, bool use_predicted_b = true, bool use_adaptive_timestep = true, double adaptive_factor = 0.125)
	{

		std::vector<Data> return_data;
		BorisIntegrator integrator = BorisIntegrator(init_pos, init_vel, E_function, B_function, t_initial, dt, use_adaptive_timestep, adaptive_factor);
		
		return_data.push_back( Data(integrator.positions, integrator.velocities, integrator.t) );

		while (integrator.t < t_final)
		{
			
			integrator.step();


			if (integrator.t <= t_final)
			{
				return_data.push_back( Data(integrator.positions, integrator.velocities, integrator.t) );	
			}
		
		}

		return return_data;

	}	
}


