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

	std::vector<Data> Integrate(Eigen::Vector3d (*E_function)(const Eigen::Vector3d &pos, const Eigen::Vector3d &vel, double t), Eigen::Vector3d (*B_function)(const Eigen::Vector3d &pos, const Eigen::Vector3d &vel, double t), double t_initial, double t_final, const Eigen::MatrixXd &init_pos, const Eigen::MatrixXd &init_vel, double dt, bool use_adaptive_timestep = true, double adaptive_factor = 0.125)
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

