#define _USE_MATH_DEFINES
#include <eigen3/Eigen/Dense>
#include <cmath>
#include <iostream>
#include <vector>
#include "integrator.h"
#include <fstream>
#include <chrono>

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
		bool use_adaptive_timestep;
		double adaptive_factor;

		double t = 0;
		double dt;
		Eigen::Vector3d position;
		Eigen::Vector3d velocity;

		std::vector<Eigen::Matrix<double, 7, 1> > b_values;
		std::vector<Eigen::Matrix<double, 7, 1> > g_values;
		Eigen::Vector3d init_pos;
		Eigen::Vector3d init_vel;
		Eigen::Vector3d init_acc;
		double init_t;

		int count = 0;
		
		Eigen::Matrix<double,7,7> c_m;
		const double h[8] = { 0.0, 0.0562625605369221464656521910318, 0.180240691736892364987579942780, 0.352624717113169637373907769648, 0.547153626330555383001448554766, 0.734210177215410531523210605558, 0.885320946839095768090359771030, 0.977520613561287501891174488626};

		Eigen::Vector3d (*E)(const Eigen::Vector3d &pos, const Eigen::Vector3d &vel, double t);
		Eigen::Vector3d (*B)(const Eigen::Vector3d &pos, const Eigen::Vector3d &vel, double t);

		RungaKuttaIntegrator(const Eigen::Vector3d &init_pos, const Eigen::Vector3d &init_vel, Eigen::Vector3d (*E_function)(const Eigen::Vector3d &pos, const Eigen::Vector3d &vel, double t), Eigen::Vector3d (*B_function)(const Eigen::Vector3d &pos, const Eigen::Vector3d &vel, double t), double init_t, double init_dt, bool predict_b = true, bool adaptive_timestep = true, double adaptive_period_factor = 0.125)
		{

			c_m <<	1.00000000000000000000, -0.0562625605369221500, 0.01014080283006363000, -0.0035758977292516170, 0.00195656540994722100, -0.0014365302363708915, 0.00127179030902686780, 
					0.00000000000000000000, 1.00000000000000000000, -0.2365032522738145200, 0.09353769525946207000, -0.0547553868890686900, 0.04215852772126870600, -0.0387603579159067700, 
					0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000, -0.5891279693869842000, 0.41588120008230690000, -0.3600995965020568000, 0.36096224345284600000, 
					0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000, -1.1362815957175396000, 1.25015071184069100000, -1.4668842084004270000, 
					0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000, -1.8704917729329500000, 2.90613625930842900000, 
					0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000, -2.7558127197720457000, 
					0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000;

			use_adaptive_timestep = adaptive_timestep;
			adaptive_factor = adaptive_period_factor;

			E = E_function;
			B = B_function;

			position = init_pos;
			velocity = init_vel;
			t = init_t;
			dt = init_dt;

		}

	private:
		//Takes in the position vector and velocity vector of a particle and returns
		//the acceleration vector of that particle
		Eigen::Vector3d ODE_function(const Eigen::Vector3d &pos, const Eigen::Vector3d &vel, double t)
		{
			Eigen::Vector3d acceleration = E(pos, vel, t) + vel.cross( B(pos, vel, t) );
			return acceleration;
		}

		//Steps all particles in the system forward to time init_t + h*dt
		void substep(double h, double dt) 
		{

			double omega = magnitude(B(init_pos, init_vel, init_t));
			Eigen::Vector3d A_ = omega * init_vel.cross(B(init_pos, init_vel, init_t)) + E(init_pos, init_vel, init_t);
			Eigen::Vector3d B_ = - init_vel * omega - E(init_pos, init_vel, init_t).cross(B(init_pos, init_vel, init_t));


			for (int j = 0; j < 2; j++)
			{
				position[j] = init_pos[j]
					+ (h*dt) * init_vel[j]
					+ (h*h*dt*dt/2)		*	(init_acc[j]
					+ (h/3)				*	(b_values[j](0) 
					+ (h/2)				*	(b_values[j](1)
					+ (3*h/5)			*	(b_values[j](2)
					+ (4*h/6)			*	(b_values[j](3)
					+ (5*h/7)			*	(b_values[j](4)
					+ (6*h/8)			*	(b_values[j](5)
					+ (7*h/9)			*	(b_values[j](6)))))))));

				velocity[j] = init_vel[j]
					+ (h*dt)	*	(init_acc[j]
					+ (h/2)		*	(b_values[j](0) 
					+ (2*h/3)	*	(b_values[j](1) 
					+ (3*h/4)	*	(b_values[j](2)
					+ (4*h/5)	*	(b_values[j](3)
					+ (5*h/6)	*	(b_values[j](4)
					+ (6*h/7)	*	(b_values[j](5)
					+ (7*h/8)	*	(b_values[j](6)))))))));
			}


			position[j] = init_pos[j]
				+ (h*dt) * init_vel[j]
				+ (h*h*dt*dt/2)		*	(init_acc[j]
				+ (h/3)				*	(b_values[j](0) 
				+ (h/2)				*	(b_values[j](1)
				+ (3*h/5)			*	(b_values[j](2)
				+ (4*h/6)			*	(b_values[j](3)
				+ (5*h/7)			*	(b_values[j](4)
				+ (6*h/8)			*	(b_values[j](5)
				+ (7*h/9)			*	(b_values[j](6)))))))));

			velocity[j] = init_vel[j]
				+ (h*dt)	*	(init_acc[j]
				+ (h/2)		*	(b_values[j](0) 
				+ (2*h/3)	*	(b_values[j](1) 
				+ (3*h/4)	*	(b_values[j](2)
				+ (4*h/5)	*	(b_values[j](3)
				+ (5*h/6)	*	(b_values[j](4)
				+ (6*h/7)	*	(b_values[j](5)
				+ (7*h/8)	*	(b_values[j](6)))))))));

			t = init_t + dt*h;

		}

		Eigen::Vector3d acceleration(double h) {
			Eigen::Vector3d acc = init_acc;

			for (int j = 0; j < 3; j ++) {
				acc[j] += h * (b_values[j][0] + h * (b_values[j][1] + h * (b_values[j][2] + h * (b_values[j][3] + h * (b_values[j][4] + h * (b_values[j][5] + h * (b_values[j][6])))))));
			}
			return acc;
		}

		//Reassigns the b_values based upon the current g_values
		void convert_g_to_b() 
		{		

			for(int j = 0; j < 3; j ++) {					
				b_values[j] = c_m * g_values[j];						
			}

		}

		//Finds the g_values based on the initial conditions and the current b_values
		double max_diff;
		void find_g_values(double dt) {

			Eigen::Vector3d next_acc;
			max_diff = 0;

			//Steps system forward to dt = h[1];
			substep(h[1], dt);
			
			next_acc = ODE_function(position, velocity, t);

			max_diff = std::max(magnitude(next_acc - acceleration(h[1])) / magnitude(next_acc), max_diff);
			
			//Calculates all g1 values
			for (int j = 0; j < 3; j ++) {

				g_values[j](0) = ( (next_acc[j] - init_acc[j]) / h[1] );

			}

			//Steps system forward to dt = h[2];
			substep(h[2], dt);

			next_acc = ODE_function(position, velocity, t);

			//Calculates all g2 values
			for (int j = 0; j < 3; j ++) {

				g_values[j](1) = ( (next_acc[j] - init_acc[j] - h[2]*g_values[j](0)) / (h[2]*(h[2]-h[1])) );

			}

			//Steps system forward to dt = h[3];
			substep(h[3], dt);

			next_acc = ODE_function(position, velocity, t);

			//Calculates all g3 values
			for (int j = 0; j < 3; j ++) {

				g_values[j](2) = ( (next_acc[j] - init_acc[j] - h[3]*(g_values[j](0) + (h[3]-h[1])*g_values[j](1))) / (h[3]*(h[3]-h[1])*(h[3]-h[2])) );

			}

			//Steps system forward to dt = h[4]
			substep(h[4], dt);

			next_acc = ODE_function(position, velocity, t);

			//Calculates all g4 values
			for (int j = 0; j < 3; j ++) {

				g_values[j](3) = ( (next_acc[j] - init_acc[j] - h[4]*(g_values[j](0) + (h[4]-h[1])*(g_values[j](1) + (h[4]-h[2])*(g_values[j](2))))) / (h[4]*(h[4]-h[1])*(h[4]-h[2])*(h[4]-h[3])) );

			}

			//Steps system forward to dt = h[5]
			substep(h[5], dt);

			next_acc = ODE_function(position, velocity, t);

			//Calculates all g5 values
			for (int j = 0; j < 3; j ++) {

				g_values[j](4) = ( (next_acc[j] - init_acc[j] - h[5]*(g_values[j](0) + (h[5]-h[1])*(g_values[j](1) + (h[5]-h[2])*(g_values[j](2) + (h[5]-h[3])*(g_values[j](3)))))) / (h[5]*(h[5]-h[1])*(h[5]-h[2])*(h[5]-h[3])*(h[5]-h[4])) );

			}

			//Steps system forward to dt = h[6]
			substep(h[6], dt);

			next_acc = ODE_function(position, velocity, t);

			//Calculates all g6 values
			for (int j = 0; j < 3; j ++) {

				g_values[j](5) = ( (next_acc[j] - init_acc[j] - h[6]*(g_values[j](0) + (h[6]-h[1])*(g_values[j](1) + (h[6]-h[2])*(g_values[j](2) + (h[6]-h[3])*(g_values[j](3) + (h[6]-h[4])*(g_values[j](4))))))) / (h[6]*(h[6]-h[1])*(h[6]-h[2])*(h[6]-h[3])*(h[6]-h[4])*(h[6]-h[5])) );

			}


			//Steps system forward to dt=h[7]
			substep(h[7], dt);

			next_acc = ODE_function(position, velocity, t);

			max_diff = std::max(magnitude(next_acc - acceleration(h[7])) / magnitude(next_acc), max_diff);

			//Calculates all g7 values
			for (int j = 0; j < 3; j ++) {

				g_values[j](6) = ( (next_acc[j] - init_acc[j] - h[7]*(g_values[j](0) + (h[7]-h[1])*(g_values[j](1) + (h[7]-h[2])*(g_values[j](2) + (h[7]-h[3])*(g_values[j](3) + (h[7]-h[4])*(g_values[j](4) + (h[7]-h[5])*(g_values[j](5)))))))) / (h[7]*(h[7]-h[1])*(h[7]-h[2])*(h[7]-h[3])*(h[7]-h[4])*(h[7]-h[5])*(h[7]-h[6])) );

			}

			//Returns system to dt = 0 and converts the g values to b values
			substep(0, dt);
			
			convert_g_to_b();
		}


		Eigen::Matrix<double, 7, 1> b_x;
		Eigen::Matrix<double, 7, 1> b_y;
		Eigen::Matrix<double, 7, 1> b_z;

		bool initialized = false;
		//Initializes the arrays on the first time-step.
		//Calculates the step-size and initial b-values for each time-step
		void initialize_step() 
		{

			init_pos = position;
			init_vel = velocity;
			init_acc = ODE_function(position,velocity,t);
			init_t = t;


			//Calculates the size of the time-step that should be taken by
			//averaging the magnetic field at the start and end of a theoretical
			//time-step
			if (use_adaptive_timestep)
			{
				double Omega = magnitude(B(init_pos,init_vel,init_t));

				if (Omega != 0)
				{
					dt = (1/Omega) * (adaptive_factor);
				}
			}

			initialize_b_values();
				
			if (not initialized)
			{

				Eigen::Matrix<double, 7, 1> coeff;

				coeff << 0, 0, 0, 0, 0, 0, 0;

				for (int j = 0; j < 3; j ++) {
					b_values.push_back(coeff);
					g_values.push_back(coeff);
				}

				initialized = true;
			
			}

			b_values[0] = b_x;
			b_values[1] = b_y;
			b_values[2] = b_z;
			//find_g_values(dt);
		}

	public:
		
		double total_count = 0;
		double time_steps_taken = 0;
		
		std::vector<Eigen::Matrix<double, 7, 1> > old_b_values;
		//The main function, that converges a predictor-corrector
		//loop and advances the stystem forward one time-step 
		void step()
		{
			
			if (count == 0)
			{
				initialize_step();
			}
			
			old_b_values = b_values;
			find_g_values(dt);

			//"If the predictor-corrector loop has converged"
			if (max_diff < 1e-14 or count > 20)
			{	
				std::cout << count << std::endl;
				substep(1, dt);
				total_count += count;
				time_steps_taken ++;
				count = 0;

			}
			else
			{
				count += 1;
				step();
			}
		}
	};

	std::vector<Data> Integrate(Eigen::Vector3d (*E_function)(const Eigen::Vector3d &pos, const Eigen::Vector3d &vel, double t), Eigen::Vector3d (*B_function)(const Eigen::Vector3d &pos, const Eigen::Vector3d &vel, double t), double t_initial, double t_final, const Eigen::Vector3d &init_pos, const Eigen::Vector3d &init_vel, double dt, bool use_predicted_b = true, bool use_adaptive_timestep = true, double adaptive_factor = 0.125)
	{
		std::vector<Data> return_data;
		RungaKuttaIntegrator integrator = RungaKuttaIntegrator(init_pos, init_vel, E_function, B_function, t_initial, dt, use_predicted_b, use_adaptive_timestep, adaptive_factor);

		return_data.push_back( Data(integrator.position, integrator.velocity, integrator.t) );

		while (integrator.t < t_final)
		{
			
			integrator.step();

			if (integrator.t <= t_final)
			{
				return_data.push_back( Data(integrator.position, integrator.velocity, integrator.t) );	
			}
		
		}

		std::cout << integrator.total_count / (integrator.time_steps_taken) << std::endl;

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
		
		Eigen::Vector3d position;
		Eigen::Vector3d velocity;

		bool use_adaptive_timestep;
		double adaptive_factor;

		Eigen::Vector3d (*E)(const Eigen::Vector3d &pos, const Eigen::Vector3d &vel, double t);
		Eigen::Vector3d (*B)(const Eigen::Vector3d &pos, const Eigen::Vector3d &vel, double t);

		BorisIntegrator(const Eigen::Vector3d &init_pos, const Eigen::Vector3d &init_vel, Eigen::Vector3d (*E_function)(const Eigen::Vector3d &pos, const Eigen::Vector3d &vel, double t), Eigen::Vector3d (*B_function)(const Eigen::Vector3d &pos, const Eigen::Vector3d &vel, double t), double init_t, double init_dt, bool adaptive_timestep = true, double adaptive_period_factor = 0.125)
		{

			position = init_pos;
			velocity = init_vel;

			E = E_function;
			B = B_function;

			t = init_t;
			dt = init_dt;

			use_adaptive_timestep = adaptive_timestep;
			adaptive_factor = adaptive_period_factor;

		}

		void step()
		{
			if (use_adaptive_timestep)
			{
				double Omega = magnitude(B(position, velocity, t));
				if (Omega != 0)
				{
					dt = (1/Omega) * (adaptive_factor);
					std::cout << dt << std::endl;
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

			

			
			{

				pos = position;
				vel = velocity;

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

			    position = pos;
			    velocity = vel;

			}

			t += dt;

		}

	};

	std::vector<Data> Integrate(Eigen::Vector3d (*E_function)(const Eigen::Vector3d &pos, const Eigen::Vector3d &vel, double t), Eigen::Vector3d (*B_function)(const Eigen::Vector3d &pos, const Eigen::Vector3d &vel, double t), double t_initial, double t_final, const Eigen::MatrixXd &init_pos, const Eigen::MatrixXd &init_vel, double dt, bool use_predicted_b = true, bool use_adaptive_timestep = true, double adaptive_factor = 0.125)
	{

		std::vector<Data> return_data;
		BorisIntegrator integrator = BorisIntegrator(init_pos, init_vel, E_function, B_function, t_initial, dt, use_adaptive_timestep, adaptive_factor);
		
		return_data.push_back( Data(integrator.position, integrator.velocity, integrator.t) );

		while (integrator.t < t_final)
		{
			
			integrator.step();


			if (integrator.t <= t_final)
			{
				return_data.push_back( Data(integrator.position, integrator.velocity, integrator.t) );	
			}
		
		}

		return return_data;

	}	
}

