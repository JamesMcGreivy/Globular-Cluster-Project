#define _USE_MATH_DEFINES
#include <eigen3/Eigen/Dense>
#include <cmath>
#include <iostream>
#include <vector>
#include "integrator.h"
#include <fstream>
#include <chrono>

namespace Integrator
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

		std::vector<Eigen::Matrix<double, 8, 1> > b_values;
		std::vector<Eigen::Matrix<double, 8, 1> > g_values;
		Eigen::Vector3d init_pos;
		Eigen::Vector3d init_vel;
		Eigen::Vector3d init_acc;
		double init_t;

		double omega;
		Eigen::Vector3d a;
		Eigen::Vector3d b;

		int count = 0;
		
		Eigen::Matrix<double,8,8> c_m;
		const double h[8] = { 0.0, 0.0562625605369221464656521910318, 0.180240691736892364987579942780, 0.352624717113169637373907769648, 0.547153626330555383001448554766, 0.734210177215410531523210605558, 0.885320946839095768090359771030, 0.977520613561287501891174488626};

		Eigen::Vector3d (*E)(const Eigen::Vector3d &pos, const Eigen::Vector3d &vel, double t);
		Eigen::Vector3d (*B)(const Eigen::Vector3d &pos, const Eigen::Vector3d &vel, double t);

		RungaKuttaIntegrator(const Eigen::Vector3d &init_pos, const Eigen::Vector3d &init_vel, Eigen::Vector3d (*E_function)(const Eigen::Vector3d &pos, const Eigen::Vector3d &vel, double t), Eigen::Vector3d (*B_function)(const Eigen::Vector3d &pos, const Eigen::Vector3d &vel, double t), double init_t, double init_dt, bool adaptive_timestep = true, double adaptive_period_factor = 0.125)
		{

			c_m << 	1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
					0.0, 1, -0.0562625605369221464656521910318, 0.010140802830063629986481804785974182308866271229767, -0.0035758977292516175949344588994050706108720804043654, 0.0019565654099472210769005670603189885571956471712273, -0.0014365302363708915424459552998590393061684596650879, 0.0012717903090268677492943116148344573439495152714550,
					0.0, 0.0, 1, -0.2365032522738145114532321338118, 0.093537695259462065895748461145361180200417881960741, -0.054755386889068686440808429439418673464885107321435, 0.042158527721268707707297347035605828637126355974817, -0.038760357915906770369904624820589821903042636522055,
					0.0, 0.0, 0.0, 1, -0.5891279693869841488271399034598, 0.41588120008230686168862191119100434500510510288145, -0.36009959650205681228976646105758538951242389978524, 0.36096224345284598322533980803401415582035565375754,
					0.0, 0.0, 0.0, 0.0, 1, -1.1362815957175395318285884582258, 1.2501507118406910258505440751095981077084355810659, -1.4668842084004269643701552583099963247634348625813,
					0.0, 0.0, 0.0, 0.0, 0.0, 1, -1.8704917729329500633517990637838, 2.9061362593084293014237913072901328720574767857664,
					0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1, -2.7558127197720458314421588348138,
					0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1;
			use_adaptive_timestep = adaptive_timestep;
			adaptive_factor = adaptive_period_factor;

			E = E_function;
			B = B_function;

			position = init_pos;
			velocity = init_vel;
			t = init_t;
			dt = init_dt;

			omega = getOmega();
			a = getA();
			b = getB();

		}

	private:
		//Takes in the position vector and velocity vector of a particle and returns
		//the acceleration vector of that particle
		Eigen::Vector3d ODE_function(const Eigen::Vector3d &pos, const Eigen::Vector3d &vel, double t)
		{
			Eigen::Vector3d acceleration = E(pos, vel, t) + vel.cross( B(pos, vel, t) );
			return acceleration;
		}

		double getOmega() {
			return magnitude(B(init_pos, init_vel, t));
		}

		Eigen::Vector3d getA() {
			return (getOmega() * init_vel).cross( B(init_pos, init_vel, t)) / magnitude(B(init_pos, init_vel, t)) + E(init_pos, init_vel, t);
		}
		
		Eigen::Vector3d getB() {
			return -(getOmega() * init_vel ) + (E(init_pos, init_vel, t).cross(B(init_pos, init_vel, t)) / magnitude(B(init_pos, init_vel, t)));
		}

		double sinVal;
		double cosVal;

		//Steps all particles in the system forward to time init_t + h*dt
		void substep(double h, double dt) 
		{

			cosVal = cos(omega * h * dt);
			sinVal = sin(omega * h * dt);
			for (int j = 0; j < 3; j++)
			{
				position[j] = ( - a[j] / pow(omega , 2) * cosVal - b[j] / pow(omega , 2) * sinVal + a[j] / pow(omega , 2) + init_pos[j] + (b[j] / omega + init_vel[j]) * (h * dt)
			+ (pow((h * dt), 2) / 2.0) * (b_values[j][0] + (h / 3.0) * (b_values[j][1] + (2.0 * h / 4.0) * (b_values[j][2] + (3.0 * h / 5.0) * (b_values[j][3] + (4.0 * h / 6.0) * (b_values[j][4] + (5.0 * h / 7.0) * (b_values[j][5] + (6.0 * h / 8.0) * (b_values[j][6] + (7.0 * h / 9.0) * (b_values[j][7]))))))))
				);

				velocity[j] = ( ( a[j] / omega * sinVal ) - ( b[j] / omega * cosVal ) + ( b[j] / omega ) + init_vel[j]
			+ dt * h * (b_values[j][0] + (h / 2.0) * (b_values[j][1] + (2.0 * h / 3.0) * (b_values[j][2] + (3.0 * h / 4.0) * (b_values[j][3] + (4.0 * h / 5.0) * (b_values[j][4] + (5.0 * h / 6.0) * (b_values[j][5] + (6.0 * h / 7.0) * (b_values[j][6] + (7.0 * h / 8.0) * (b_values[j][7]))))))))
				);
			}

			t = init_t + dt*h;

		}

		Eigen::Vector3d acceleration(double h) {
			Eigen::Vector3d acc = init_acc;

			cosVal = cos(omega * h * dt);
			sinVal = sin(omega * h * dt);
			for (int j = 0; j < 3; j ++) {
				acc[j] = ((a[j] * cosVal + b[j] * sinVal)
			+ b_values[j][0] + (h) * (b_values[j][1] + (h) * (b_values[j][2] + (h) * (b_values[j][3] + (h) * (b_values[j][4] + (h) * (b_values[j][5] + (h) * (b_values[j][6] + (h) * (b_values[j][7])))))))
				);
			}
			return acc;
		}

		//Finds the g_values based on the initial conditions and the current b_values
		void find_g_values(double dt) {

			Eigen::Vector3d next_acc;
			double nextDenom;

			//Steps system forward to dt = h[1];
			substep(h[1], dt);

			next_acc = ODE_function(position, velocity, t);

			nextDenom = ((h[1] - h[0]));
			//Calculates all g1 values
			cosVal = cos(omega * h[1] * dt);
			sinVal = sin(omega * h[1] * dt);
			for (int j = 0; j < 3; j ++) {
				g_values[j][1] = ( next_acc[j] / nextDenom ) - ((a[j] * cosVal + b[j] * sinVal) / nextDenom ) - ((g_values[j][0] ) / nextDenom );
			}

			//Steps system forward to dt = h[2];
			substep(h[2], dt);

			next_acc = ODE_function(position, velocity, t);

			nextDenom = ((h[2] - h[0]) * (h[2] - h[1]));
			//Calculates all g2 values
			cosVal = cos(omega * h[2] * dt);
			sinVal = sin(omega * h[2] * dt);
			for (int j = 0; j < 3; j ++) {
				g_values[j][2] = ( next_acc[j] / nextDenom ) - ((a[j] * cosVal + b[j] * sinVal) / nextDenom ) - ((g_values[j][0] + ((h[2] - h[0])) * (g_values[j][1])) / nextDenom );
			}

			//Steps system forward to dt = h[3]
			substep(h[3], dt);

			next_acc = ODE_function(position, velocity, t);

			nextDenom = ((h[3] - h[0]) * (h[3] - h[1]) * (h[3] - h[2]));
			//Calculates all g3 values
			cosVal = cos(omega * h[3] * dt);
			sinVal = sin(omega * h[3] * dt);
			for (int j = 0; j < 3; j ++) {
				g_values[j][3] = ( next_acc[j] / nextDenom ) - ((a[j] * cosVal + b[j] * sinVal) / nextDenom ) - ((g_values[j][0] + ((h[3] - h[0])) * (g_values[j][1] + ((h[3] - h[1])) * (g_values[j][2]))) / nextDenom );
			}

			//Steps system forward to dt = h[4]
			substep(h[4], dt);

			next_acc = ODE_function(position, velocity, t);

			nextDenom = ((h[4] - h[0]) * (h[4] - h[1]) * (h[4] - h[2]) * (h[4] - h[3]));
			//Calculates all g4 values
			cosVal = cos(omega * h[4] * dt);
			sinVal = sin(omega * h[4] * dt);
			for (int j = 0; j < 3; j ++) {
				g_values[j][4] = ( next_acc[j] / nextDenom ) - ((a[j] * cosVal + b[j] * sinVal) / nextDenom ) - ((g_values[j][0] + ((h[4] - h[0])) * (g_values[j][1] + ((h[4] - h[1])) * (g_values[j][2] + ((h[4] - h[2])) * (g_values[j][3])))) / nextDenom );
			}

			//Steps system forward to dt = h[5]
			substep(h[5], dt);

			next_acc = ODE_function(position, velocity, t);

			nextDenom = ((h[5] - h[0]) * (h[5] - h[1]) * (h[5] - h[2]) * (h[5] - h[3]) * (h[5] - h[4]));
			//Calculates all g5 values
			cosVal = cos(omega * h[5] * dt);
			sinVal = sin(omega * h[5] * dt);
			for (int j = 0; j < 3; j ++) {
				g_values[j][5] = ( next_acc[j] / nextDenom ) - ((a[j] * cosVal + b[j] * sinVal) / nextDenom ) - ((g_values[j][0] + ((h[5] - h[0])) * (g_values[j][1] + ((h[5] - h[1])) * (g_values[j][2] + ((h[5] - h[2])) * (g_values[j][3] + ((h[5] - h[3])) * (g_values[j][4]))))) / nextDenom );
			}


			//Steps system forward to dt=h[6]
			substep(h[6], dt);

			next_acc = ODE_function(position, velocity, t);

			nextDenom = ((h[6] - h[0]) * (h[6] - h[1]) * (h[6] - h[2]) * (h[6] - h[3]) * (h[6] - h[4]) * (h[6] - h[5]));
			//Calculates all g6 values
			cosVal = cos(omega * h[6] * dt);
			sinVal = sin(omega * h[6] * dt);
			for (int j = 0; j < 3; j ++) {
				g_values[j][6] = ( next_acc[j] / nextDenom ) - ((a[j] * cosVal + b[j] * sinVal) / nextDenom ) - ((g_values[j][0] + ((h[6] - h[0])) * (g_values[j][1] + ((h[6] - h[1])) * (g_values[j][2] + ((h[6] - h[2])) * (g_values[j][3] + ((h[6] - h[3])) * (g_values[j][4] + ((h[6] - h[4])) * (g_values[j][5])))))) / nextDenom );
			}

			//Steps system forward to dt = h[7];
			substep(h[7], dt);

			next_acc = ODE_function(position, velocity, t);

			nextDenom = ((h[7] - h[0]) * (h[7] - h[1]) * (h[7] - h[2]) * (h[7] - h[3]) * (h[7] - h[4]) * (h[7] - h[5]) * (h[7] - h[6]));
			//Calculates all g7 values
			cosVal = cos(omega * h[7] * dt);
			sinVal = sin(omega * h[7] * dt);
			for (int j = 0; j < 3; j ++) {
				g_values[j][7] = ( next_acc[j] / nextDenom ) - ((a[j] * cosVal + b[j] * sinVal) / nextDenom ) - ((g_values[j][0] + ((h[7] - h[0])) * (g_values[j][1] + ((h[7] - h[1])) * (g_values[j][2] + ((h[7] - h[2])) * (g_values[j][3] + ((h[7] - h[3])) * (g_values[j][4] + ((h[7] - h[4])) * (g_values[j][5] + ((h[7] - h[5])) * (g_values[j][6]))))))) / nextDenom );
			}

			//Returns system to dt = 0 and converts the g values to b values
			substep(0, dt);
			
			for (int j = 0; j < 3; j++) {
				b_values[j] = c_m * g_values[j];
			}
		}


		bool initialized = false;
		//Initializes the arrays on the first time-step.
		void initialize_step() 
		{

			init_pos = position;
			init_vel = velocity;
			init_acc = ODE_function(position,velocity,t);
			init_t = t;

			omega = getOmega();
			a = getA();
			b = getB();


			//Calculates the size of the time-step that should be taken by
			//averaging the magnetic field at the start and end of a theoretical
			//time-step
				
			if (not initialized)
			{

				Eigen::Matrix<double, 8, 1> coeff;

				coeff << 0, 0, 0, 0, 0, 0, 0, 0;

				for (int j = 0; j < 3; j ++) {
					b_values.push_back(coeff);
					g_values.push_back(coeff);
				}

				initialized = true;
			
			}

			if (use_adaptive_timestep)
			{
				double Omega1 = magnitude(B(init_pos,init_vel,init_t));
				substep(1, dt);
				double Omega2 = magnitude(B(position,velocity,t));
				substep(0, dt);

				if (Omega1 != 0 && Omega2 != 0)
				{
					dt = std::min(adaptive_factor * 2 * M_PI / Omega1, adaptive_factor * 2 * M_PI / Omega2);
				}
			}
		}

	public:
		
		double total_count = 0;
		double time_steps_taken = 0;

		Eigen::Vector3d next_acc;
		double converged;
		//The main function, that converges a predictor-corrector
		//loop and advances the stystem forward one time-step 
		void step()
		{
			
			if (count == 0)
			{
				initialize_step();
			}
			
			find_g_values(dt);
			count += 1;
			
			substep(h[7], dt);
			next_acc = ODE_function(position, velocity, t);
			converged = magnitude(next_acc - acceleration(h[7]));

			//"If the predictor-corrector loop has converged"
			if (converged <= 1e-14 or count >= 20)
			{	
				substep(1, dt);
				total_count += count;
				time_steps_taken ++;
				count = 0;

			}
			else
			{
				step();
			}
		}
	};

	std::vector<Data> Integrate(Eigen::Vector3d (*E_function)(const Eigen::Vector3d &pos, const Eigen::Vector3d &vel, double t), Eigen::Vector3d (*B_function)(const Eigen::Vector3d &pos, const Eigen::Vector3d &vel, double t), double t_initial, double t_final, const Eigen::MatrixXd &init_pos, const Eigen::MatrixXd &init_vel, double dt, bool use_adaptive_timestep = true, double adaptive_factor = 0.125)
	{
		std::vector<Data> return_data;
		RungaKuttaIntegrator integrator = RungaKuttaIntegrator(init_pos, init_vel, E_function, B_function, t_initial, dt, use_adaptive_timestep, adaptive_factor);

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
