
namespace ODE_Integrator
{

	class Integrator
	{
	public:
		double t = 0;
		double dt;
		int num_particles;
		Eigen::MatrixXd positions;
		Eigen::MatrixXd velocities;
		Eigen::MatrixXd (*ODE_function)(const Eigen::MatrixXd &positions, const Eigen::MatrixXd &velocities, double t);

		Integrator(const Eigen::MatrixXd &init_pos, const Eigen::MatrixXd &init_vel, double init_t, double init_dt, Eigen::MatrixXd (*y_pp_function)(const Eigen::MatrixXd &pos, const Eigen::MatrixXd &vel, double t))
		{
			positions = init_pos;
			velocities = init_vel;
			num_particles = positions.rows();
			ODE_function = y_pp_function;
			t = init_t;
			dt = init_dt;
		}

		virtual void step() {std::cout << "ERROR" << std::endl; t += dt;}


	};


	class RungaKuttaIntegrator : public Integrator
	{
	private:
		std::vector<Eigen::Matrix<double,8,3> > b_values;
		std::vector<Eigen::Matrix<double,8,3> > g_values;
		Eigen::MatrixXd init_pos;
		Eigen::MatrixXd init_vel;
		Eigen::MatrixXd init_acc;

		std::vector<Eigen::Matrix<double,8,3> > old_b_values;
		int count = 0;
		double epsilon = 0.028;
		double previous_error = 0;
		
		Eigen::Matrix<double,8,8> c_m;
		const double h[8] = { 0.0, 0.0562625605369221464656521910318, 0.180240691736892364987579942780, 0.352624717113169637373907769648, 0.547153626330555383001448554766, 0.734210177215410531523210605558, 0.885320946839095768090359771030, 0.977520613561287501891174488626};



	public:
		bool adaptive_timestepping;

		RungaKuttaIntegrator(const Eigen::MatrixXd &init_pos, const Eigen::MatrixXd &init_vel, double init_t, double init_dt, Eigen::MatrixXd (*y_pp_function)(const Eigen::MatrixXd &pos, const Eigen::MatrixXd &vel, double t), bool at = false)
		:Integrator(init_pos, init_vel, init_t, init_dt, y_pp_function) 
		{ 
			c_m <<	1.00000000000000000000, -0.0562625605369221500, 0.01014080283006363000, -0.0035758977292516170, 0.00195656540994722100, -0.0014365302363708915, 0.00127179030902686780, -0.0012432012432012432, 
					0.00000000000000000000, 1.00000000000000000000, -0.2365032522738145200, 0.09353769525946207000, -0.0547553868890686900, 0.04215852772126870600, -0.0387603579159067700, 0.03916083916083916000, 
					0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000, -0.5891279693869842000, 0.41588120008230690000, -0.3600995965020568000, 0.36096224345284600000, -0.3916083916083916000, 
					0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000, -1.1362815957175396000, 1.25015071184069100000, -1.4668842084004270000, 1.79487179487179470000, 
					0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000, -1.8704917729329500000, 2.90613625930842900000, -4.3076923076923075000, 
					0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000, -2.7558127197720457000, 5.60000000000000000000, 
					0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000, -3.7333333333333334000, 
					0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000;
			
			adaptive_timestepping = at;
		
		}

	private:
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
						+ (7*h/9)			*	(b_values[i](6,j)
						+ (8*h/10)			*	(b_values[i](7,j))))))))));

					vel[j] = init_vel(i,j)
						+ (h*dt)	*	(init_acc(i,j)
						+ (h/2)		*	(b_values[i](0,j) 
						+ (2*h/3)	*	(b_values[i](1,j) 
						+ (3*h/4)	*	(b_values[i](2,j)
						+ (4*h/5)	*	(b_values[i](3,j)
						+ (5*h/6)	*	(b_values[i](4,j)
						+ (6*h/7)	*	(b_values[i](5,j)
						+ (7*h/8)	*	(b_values[i](6,j)
						+ (8*h/9)	*	(b_values[i](7,j))))))))));
				}

				positions.row(i) = pos;
				velocities.row(i) = vel;

			}
		}

		//Reassigns the b_values based upon the current g_values
		void convert_g_to_b() 
		{		
			
			Eigen::Matrix<double,8,3> new_b;
				
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

			//Steps system forward to dt = dt
			substep(h[8], dt);
			next_acc = ODE_function(positions, velocities, t);

			//Calculates all g8 values
			for (int i = 0; i < num_particles; i ++) {

				for (int j = 0; j < 3; j ++) {

					g_values[i](7,j) = ( (next_acc(i,j) - init_acc(i,j) - h[8]*(g_values[i](0,j) + (h[8]-h[1])*(g_values[i](1,j) + (h[8]-h[2])*(g_values[i](2,j) + (h[8]-h[3])*(g_values[i](3,j) + (h[8]-h[4])*(g_values[i](4,j) + (h[8]-h[5])*(g_values[i](5,j) + (h[8]-h[6])*(g_values[i](6,j))))))))) / (h[8]*(h[8]-h[1])*(h[8]-h[2])*(h[8]-h[3])*(h[8]-h[4])*(h[8]-h[5])*(h[8]-h[6])*(h[8]-h[7])) );

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

				initialized = true;

			}

			find_g_values(dt);

		}

	public:
		void step()
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

			double global_error = max_del_b8 / max_y_pp;
			if (count > 60)
			{

				if (adaptive_timestepping)
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
					if (not isnormal(dt_req))
					{
						dt_req = 1.5*dt;
					}

					if (dt > dt_req)
					{
						count = 0;
						dt = dt_req;
						step();
					}
					
					else
					{
						substep(1, dt);
						std::cout << "converged" << dt << std::endl;
						count = 0;
						t += dt;
						dt = dt_req;
						
					}
				}
				else
				{
					//std::cout << "count: " << count << std::endl;
					substep(1, dt);
					count = 0;
					t += dt;

				}

			}

			else
			{
				count += 1;
				previous_error = global_error;
				step();
			}
		}

	};	

	//fun, t_span, y0, method='RK45', t_eval=None, dense_output=False, events=None, vectorized=False, args=None, **options
	std::vector<Data> Integrate(Eigen::MatrixXd (*y_pp_function)(const Eigen::MatrixXd &positions, const Eigen::MatrixXd &velocities, double t), double t_initial, double t_final, const Eigen::MatrixXd &init_positions, const Eigen::MatrixXd &init_velocities, double dt, const Eigen::VectorXd &t_eval = Eigen::VectorXd(), std::string method = "RK")
	{

		std::vector<Data> return_data;
		RungaKuttaIntegrator *integrator;

		//To be used if t_eval is provided by user
		int max_t_eval_index = t_eval.rows();
		bool use_t_eval = max_t_eval_index;
		int t_eval_index = 0;
		double dt_to_t_eval;
		
		//Chooses the integration method used
		integrator = new RungaKuttaIntegrator(init_positions, init_velocities, t_initial, dt, y_pp_function);

		return_data.push_back( Data(integrator->positions, integrator->velocities, integrator->t) );

		while (integrator->t < t_final)
		{

			if (not integrator->adaptive_timestepping)
			{
				integrator->dt = dt;
			}
			if (use_t_eval and integrator->t + integrator->dt > t_eval[t_eval_index])
			{
				integrator->dt = (t_eval[t_eval_index] - integrator->t);
			}
			if (integrator->t + integrator->dt > t_final)
			{
				integrator->dt = (t_final - integrator->t);
			}

			integrator->step();
			
			if (use_t_eval and integrator->t == t_eval[t_eval_index])
			{
				return_data.push_back( Data(integrator->positions, integrator->velocities, integrator->t) );
				t_eval_index += 1;
			}

			if (integrator->t == t_final)
			{
				return_data.push_back( Data(integrator->positions, integrator->velocities, integrator->t) );	
			}
			
			if (t_eval_index == max_t_eval_index)
			{
				use_t_eval = false;
			}

		
		}
		return return_data;

	}
	
}

namespace CoupledDifEqIntegrator
{

	class Integrator
	{
	public:
		double t = 0;
		double dt;
		int num_particles;
		Eigen::MatrixXd q_values;
		Eigen::MatrixXd p_values;

		Eigen::MatrixXd (*p_dot)(const Eigen::MatrixXd &q, const Eigen::MatrixXd &p, double t);
		Eigen::MatrixXd (*q_dot)(const Eigen::MatrixXd &q, const Eigen::MatrixXd &p, double t);

		Integrator(const Eigen::MatrixXd &init_q_values, const Eigen::MatrixXd &init_p_values, double init_t, double init_dt, Eigen::MatrixXd (*q_dot_function)(const Eigen::MatrixXd &p, const Eigen::MatrixXd &q, double t), Eigen::MatrixXd (*p_dot_function)(const Eigen::MatrixXd &p, const Eigen::MatrixXd &q, double t))
		{
			q_values = init_q_values;
			p_values = init_p_values;

			
			q_dot = q_dot_function;
			p_dot = p_dot_function;
			
			num_particles = q_values.rows();
			t = init_t;
			dt = init_dt;
		}

		virtual void step() {std::cout << "ERROR" << std::endl; t += dt;}

	};


	class RK_Integrator : public Integrator
	{
	public:
		std::vector<Eigen::Matrix<double,8,3> > q_b_values;
		std::vector<Eigen::Matrix<double,8,3> > q_g_values;

		std::vector<Eigen::Matrix<double,8,3> > p_b_values;
		std::vector<Eigen::Matrix<double,8,3> > p_g_values;
		
		Eigen::MatrixXd init_p;
		Eigen::MatrixXd init_p_dot;
		Eigen::MatrixXd init_q;
		Eigen::MatrixXd init_q_dot;

		Eigen::Matrix<double,8,8> c_m;
		const double h[8] = { 0.0, 0.0562625605369221464656521910318, 0.180240691736892364987579942780, 0.352624717113169637373907769648, 0.547153626330555383001448554766, 0.734210177215410531523210605558, 0.885320946839095768090359771030, 0.977520613561287501891174488626};

		RK_Integrator(const Eigen::MatrixXd &init_q_values, const Eigen::MatrixXd &init_p_values, double init_t, double init_dt, Eigen::MatrixXd (*q_dot_function)(const Eigen::MatrixXd &p, const Eigen::MatrixXd &q, double t), Eigen::MatrixXd (*p_dot_function)(const Eigen::MatrixXd &p, const Eigen::MatrixXd &q, double t))
		:Integrator(init_q_values, init_p_values, init_t, init_dt, q_dot_function, p_dot_function) 
		{ 
			c_m <<	1.00000000000000000000, -0.0562625605369221500, 0.01014080283006363000, -0.0035758977292516170, 0.00195656540994722100, -0.0014365302363708915, 0.00127179030902686780, -0.0012432012432012432, 
					0.00000000000000000000, 1.00000000000000000000, -0.2365032522738145200, 0.09353769525946207000, -0.0547553868890686900, 0.04215852772126870600, -0.0387603579159067700, 0.03916083916083916000, 
					0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000, -0.5891279693869842000, 0.41588120008230690000, -0.3600995965020568000, 0.36096224345284600000, -0.3916083916083916000, 
					0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000, -1.1362815957175396000, 1.25015071184069100000, -1.4668842084004270000, 1.79487179487179470000, 
					0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000, -1.8704917729329500000, 2.90613625930842900000, -4.3076923076923075000, 
					0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000, -2.7558127197720457000, 5.60000000000000000000, 
					0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000, -3.7333333333333334000, 
					0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000;
		}

	private:
		void substep(double h, double dt) 
		{

			Eigen::Vector3d q(0,0,0);
			Eigen::Vector3d p(0,0,0);

			for (int i = 0; i < num_particles; i++) 
			{

				for (int j = 0; j < 3; j++) 
				{
					q[j] = init_q(i,j)
						+ (h*dt)	*	(init_q_dot(i,j)
						+ (h/2)		*	(q_b_values[i](0,j) 
						+ (2*h/3)	*	(q_b_values[i](1,j) 
						+ (3*h/4)	*	(q_b_values[i](2,j)
						+ (4*h/5)	*	(q_b_values[i](3,j)
						+ (5*h/6)	*	(q_b_values[i](4,j)
						+ (6*h/7)	*	(q_b_values[i](5,j)
						+ (7*h/8)	*	(q_b_values[i](6,j)
						+ (8*h/9)	*	(q_b_values[i](7,j))))))))));

					p[j] = init_p(i,j)
						+ (h*dt)	*	(init_p_dot(i,j)
						+ (h/2)		*	(p_b_values[i](0,j) 
						+ (2*h/3)	*	(p_b_values[i](1,j) 
						+ (3*h/4)	*	(p_b_values[i](2,j)
						+ (4*h/5)	*	(p_b_values[i](3,j)
						+ (5*h/6)	*	(p_b_values[i](4,j)
						+ (6*h/7)	*	(p_b_values[i](5,j)
						+ (7*h/8)	*	(p_b_values[i](6,j)
						+ (8*h/9)	*	(p_b_values[i](7,j))))))))));

				}

				q_values.row(i) = q;
				p_values.row(i) = p;

			}

		}

		//Reassigns the b_values based upon the current g_values
		void convert_q_g_to_b() 
		{		
			
			Eigen::Matrix<double,8,3> new_q_b;
				
			for(int i = 0; i < num_particles; i ++) 
			{
				for(int j = 0; j < 3; j ++) 
				{	
					new_q_b.col(j) = c_m * q_g_values[i].col(j);
						
				}
				
				q_b_values[i] = new_q_b;
			
			}
		}

		void convert_p_g_to_b() 
		{		
			
			Eigen::Matrix<double,8,3> new_p_b;
				
			for(int i = 0; i < num_particles; i ++) 
			{
				for(int j = 0; j < 3; j ++) 
				{	
					new_p_b.col(j) = c_m * p_g_values[i].col(j);
						
				}

				p_b_values[i] = new_p_b;


			
			}
		}

		//Finds the g_values based on the initial conditions and the current b_values
		void find_g_values(double dt) {
			
			Eigen::MatrixXd next_p_dot;
			Eigen::MatrixXd next_q_dot;

			//Steps system forward to dt = h[1];
			substep(h[1], dt);
			next_p_dot = p_dot(q_values,p_values,t);
			
			//Calculates all g1 values
			for (int i = 0; i < num_particles; i ++) {

				for (int j = 0; j < 3; j ++) {

					p_g_values[i](0,j) = ( (next_p_dot(i,j) - init_p_dot(i,j)) / h[1] );

				}
			}

			//Steps system forward to dt = h[2];
			substep(h[2], dt);
			next_p_dot = p_dot(q_values,p_values,t);


			//Calculates all g2 values
			for (int i = 0; i < num_particles; i ++) {

				for (int j = 0; j < 3; j ++) {

					p_g_values[i](1,j) = ( (next_p_dot(i,j) - init_p_dot(i,j) - h[2]*p_g_values[i](0,j)) / (h[2]*(h[2]-h[1])) );

				}
			}

			//Steps system forward to dt = h[3];
			substep(h[3], dt);
			next_p_dot = p_dot(q_values,p_values,t);

			//Calculates all g3 values
			for (int i = 0; i < num_particles; i ++) {

				for (int j = 0; j < 3; j ++) {

					p_g_values[i](2,j) = ( (next_p_dot(i,j) - init_p_dot(i,j) - h[3]*(p_g_values[i](0,j) + (h[3]-h[1])*p_g_values[i](1,j))) / (h[3]*(h[3]-h[1])*(h[3]-h[2])) );

				}
			}

			//Steps system forward to dt = h[4]
			substep(h[4], dt);
			next_p_dot = p_dot(q_values,p_values,t);

			//Calculates all g4 values
			for (int i = 0; i < num_particles; i ++) {

				for (int j = 0; j < 3; j ++) {

					p_g_values[i](3,j) = ( (next_p_dot(i,j) - init_p_dot(i,j) - h[4]*(p_g_values[i](0,j) + (h[4]-h[1])*(p_g_values[i](1,j) + (h[4]-h[2])*(p_g_values[i](2,j))))) / (h[4]*(h[4]-h[1])*(h[4]-h[2])*(h[4]-h[3])) );

				}
			}

			//Steps system forward to dt = h[5]
			substep(h[5], dt);
			next_p_dot = p_dot(q_values,p_values,t);

			//Calculates all g5 values
			for (int i = 0; i < num_particles; i ++) {

				for (int j = 0; j < 3; j ++) {

					p_g_values[i](4,j) = ( (next_p_dot(i,j) - init_p_dot(i,j) - h[5]*(p_g_values[i](0,j) + (h[5]-h[1])*(p_g_values[i](1,j) + (h[5]-h[2])*(p_g_values[i](2,j) + (h[5]-h[3])*(p_g_values[i](3,j)))))) / (h[5]*(h[5]-h[1])*(h[5]-h[2])*(h[5]-h[3])*(h[5]-h[4])) );

				}
			}

			//Steps system forward to dt = h[6]
			substep(h[6], dt);
			next_p_dot = p_dot(q_values,p_values,t);
			
			//Calculates all g6 values
			for (int i = 0; i < num_particles; i ++) {

				for (int j = 0; j < 3; j ++) {

					p_g_values[i](5,j) = ( (next_p_dot(i,j) - init_p_dot(i,j) - h[6]*(p_g_values[i](0,j) + (h[6]-h[1])*(p_g_values[i](1,j) + (h[6]-h[2])*(p_g_values[i](2,j) + (h[6]-h[3])*(p_g_values[i](3,j) + (h[6]-h[4])*(p_g_values[i](4,j))))))) / (h[6]*(h[6]-h[1])*(h[6]-h[2])*(h[6]-h[3])*(h[6]-h[4])*(h[6]-h[5])) );

				}

			}

			//Steps system forward to dt=h[7]
			substep(h[7], dt);
			next_p_dot = p_dot(q_values,p_values,t);

			//Calculates all g7 values
			for (int i = 0; i < num_particles; i ++) {

				for (int j = 0; j < 3; j ++) {

					p_g_values[i](6,j) = ( (next_p_dot(i,j) - init_p_dot(i,j) - h[7]*(p_g_values[i](0,j) + (h[7]-h[1])*(p_g_values[i](1,j) + (h[7]-h[2])*(p_g_values[i](2,j) + (h[7]-h[3])*(p_g_values[i](3,j) + (h[7]-h[4])*(p_g_values[i](4,j) + (h[7]-h[5])*(p_g_values[i](5,j)))))))) / (h[7]*(h[7]-h[1])*(h[7]-h[2])*(h[7]-h[3])*(h[7]-h[4])*(h[7]-h[5])*(h[7]-h[6])) );

				}

			}

			//Steps system forward to dt = dt
			substep(1, dt);
			next_p_dot = p_dot(q_values,p_values,t);

			//Calculates all g8 values
			for (int i = 0; i < num_particles; i ++) {

				for (int j = 0; j < 3; j ++) {

					p_g_values[i](7,j) = ( (next_p_dot(i,j) - init_p_dot(i,j) - h[8]*(p_g_values[i](0,j) + (h[8]-h[1])*(p_g_values[i](1,j) + (h[8]-h[2])*(p_g_values[i](2,j) + (h[8]-h[3])*(p_g_values[i](3,j) + (h[8]-h[4])*(p_g_values[i](4,j) + (h[8]-h[5])*(p_g_values[i](5,j) + (h[8]-h[6])*(p_g_values[i](6,j))))))))) / (h[8]*(h[8]-h[1])*(h[8]-h[2])*(h[8]-h[3])*(h[8]-h[4])*(h[8]-h[5])*(h[8]-h[6])*(h[8]-h[7])) );

				}

			}

			//Returns system to dt = 0 and converts the g values to b values
			substep(0, dt);
			convert_p_g_to_b();


			//Steps system forward to dt = h[1];
			substep(h[1], dt);
			next_q_dot = q_dot(q_values,p_values,t);
			
			//Calculates all g1 values
			for (int i = 0; i < num_particles; i ++) {

				for (int j = 0; j < 3; j ++) {
	
					q_g_values[i](0,j) = ( (next_q_dot(i,j) - init_q_dot(i,j)) / h[1] );

				}
			}

			//Steps system forward to dt = h[2];
			substep(h[2], dt);
			next_q_dot = q_dot(q_values,p_values,t);

			//Calculates all g2 values
			for (int i = 0; i < num_particles; i ++) {

				for (int j = 0; j < 3; j ++) {

					q_g_values[i](1,j) = ( (next_q_dot(i,j) - init_q_dot(i,j) - h[2]*q_g_values[i](0,j)) / (h[2]*(h[2]-h[1])) );

				}
			}

			//Steps system forward to dt = h[3];
			substep(h[3], dt);
			next_q_dot = q_dot(q_values,p_values,t);

			//Calculates all g3 values
			for (int i = 0; i < num_particles; i ++) {

				for (int j = 0; j < 3; j ++) {

					q_g_values[i](2,j) = ( (next_q_dot(i,j) - init_q_dot(i,j) - h[3]*(q_g_values[i](0,j) + (h[3]-h[1])*q_g_values[i](1,j))) / (h[3]*(h[3]-h[1])*(h[3]-h[2])) );

				}
			}

			//Steps system forward to dt = h[4]
			substep(h[4], dt);
			next_q_dot = q_dot(q_values,p_values,t);

			//Calculates all g4 values
			for (int i = 0; i < num_particles; i ++) {

				for (int j = 0; j < 3; j ++) {

					q_g_values[i](3,j) = ( (next_q_dot(i,j) - init_q_dot(i,j) - h[4]*(q_g_values[i](0,j) + (h[4]-h[1])*(q_g_values[i](1,j) + (h[4]-h[2])*(q_g_values[i](2,j))))) / (h[4]*(h[4]-h[1])*(h[4]-h[2])*(h[4]-h[3])) );

				}
			}

			//Steps system forward to dt = h[5]
			substep(h[5], dt);
			next_q_dot = q_dot(q_values,p_values,t);

			//Calculates all g5 values
			for (int i = 0; i < num_particles; i ++) {

				for (int j = 0; j < 3; j ++) {

					q_g_values[i](4,j) = ( (next_q_dot(i,j) - init_q_dot(i,j) - h[5]*(q_g_values[i](0,j) + (h[5]-h[1])*(q_g_values[i](1,j) + (h[5]-h[2])*(q_g_values[i](2,j) + (h[5]-h[3])*(q_g_values[i](3,j)))))) / (h[5]*(h[5]-h[1])*(h[5]-h[2])*(h[5]-h[3])*(h[5]-h[4])) );

				}
			}

			//Steps system forward to dt = h[6]
			substep(h[6], dt);
			next_q_dot = q_dot(q_values,p_values,t);

			//Calculates all g6 values
			for (int i = 0; i < num_particles; i ++) {

				for (int j = 0; j < 3; j ++) {

					q_g_values[i](5,j) = ( (next_q_dot(i,j) - init_q_dot(i,j) - h[6]*(q_g_values[i](0,j) + (h[6]-h[1])*(q_g_values[i](1,j) + (h[6]-h[2])*(q_g_values[i](2,j) + (h[6]-h[3])*(q_g_values[i](3,j) + (h[6]-h[4])*(q_g_values[i](4,j))))))) / (h[6]*(h[6]-h[1])*(h[6]-h[2])*(h[6]-h[3])*(h[6]-h[4])*(h[6]-h[5])) );

				}

			}

			//Steps system forward to dt=h[7]
			substep(h[7], dt);
			next_q_dot = q_dot(q_values,p_values,t);

			//Calculates all g7 values
			for (int i = 0; i < num_particles; i ++) {

				for (int j = 0; j < 3; j ++) {

					q_g_values[i](6,j) = ( (next_q_dot(i,j) - init_q_dot(i,j) - h[7]*(q_g_values[i](0,j) + (h[7]-h[1])*(q_g_values[i](1,j) + (h[7]-h[2])*(q_g_values[i](2,j) + (h[7]-h[3])*(q_g_values[i](3,j) + (h[7]-h[4])*(q_g_values[i](4,j) + (h[7]-h[5])*(q_g_values[i](5,j)))))))) / (h[7]*(h[7]-h[1])*(h[7]-h[2])*(h[7]-h[3])*(h[7]-h[4])*(h[7]-h[5])*(h[7]-h[6])) );

				}

			}

			//Steps system forward to dt = dt
			substep(1, dt);
			next_q_dot = q_dot(q_values,p_values,t);

			//Calculates all g8 values
			for (int i = 0; i < num_particles; i ++) {

				for (int j = 0; j < 3; j ++) {

					p_g_values[i](7,j) = ( (next_p_dot(i,j) - init_p_dot(i,j) - h[8]*(p_g_values[i](0,j) + (h[8]-h[1])*(p_g_values[i](1,j) + (h[8]-h[2])*(p_g_values[i](2,j) + (h[8]-h[3])*(p_g_values[i](3,j) + (h[8]-h[4])*(p_g_values[i](4,j) + (h[8]-h[5])*(p_g_values[i](5,j) + (h[8]-h[6])*(p_g_values[i](6,j))))))))) / (h[8]*(h[8]-h[1])*(h[8]-h[2])*(h[8]-h[3])*(h[8]-h[4])*(h[8]-h[5])*(h[8]-h[6])*(h[8]-h[7])) );

				}

			}

			//Returns system to dt = 0 and converts the g values to b values
			substep(0, dt);
			convert_q_g_to_b();

		}

		bool initialized = false;
		void initialize_RK_step(double dt) 
		{

			init_q = q_values;
			init_p = p_values;
			init_q_dot = q_dot(q_values,p_values,t);
			init_p_dot = p_dot(q_values,p_values,t);

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
					q_b_values.push_back(coeff);
					q_g_values.push_back(coeff);
					p_b_values.push_back(coeff);
					p_g_values.push_back(coeff);
				}

				initialized = true;

			}

			find_g_values(dt);

		}


	public:
		int count = 0;
		void step()
		{

			if (count == 0)
			{
				initialize_RK_step(dt);
			}

			find_g_values(dt);



			if (count >= 5)
			{
				substep(1, dt);
				count = 0;
				t += dt;
				std::cout << p_b_values[0] << std::endl;

			}

			else
			{
				count += 1;
				step();
			}
		}

	};	


	std::vector<Data> Integrate(Eigen::MatrixXd (*q_dot_function)(const Eigen::MatrixXd &q, const Eigen::MatrixXd &p, double t), Eigen::MatrixXd (*p_dot_function)(const Eigen::MatrixXd &q, const Eigen::MatrixXd &p, double t), double t_initial, double t_final, const Eigen::MatrixXd &init_p, const Eigen::MatrixXd &init_q, double dt)
	{

		std::vector<Data> return_data;
		Integrator *integrator;
		
		//Chooses the integration method used
		integrator = new RK_Integrator(init_q, init_p, t_initial, dt, q_dot_function, p_dot_function);
		return_data.push_back( Data(integrator->q_values, integrator->p_values, integrator->t) );

		while (integrator->t < t_final)
		{

			integrator->step();

			if (integrator->t <= t_final)
			{
				return_data.push_back( Data(integrator->q_values, integrator->p_values, integrator->t) );	
			}
		
		}
		return return_data;

	}

}
