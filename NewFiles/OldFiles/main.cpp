#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
#include <random>
using namespace std;

double t = 0;
double dt = 0;
int force_eval_counter = 0;

//Overloading the Vector class and declaring a Particle and System class

//Overload the + , - , * operators and create a to_string function for <double> vectors
vector<double> operator+(const vector<double> &v1,const vector<double> &v2) {

	int size = v1.size();
	vector<double> v3(size);

	for (int i = 0; i < size; i ++) {
		v3[i] = v1[i] + v2[i];
	}

	return v3;
}

vector<double> operator-(const vector<double> &v1,const vector<double> &v2) {

	int size = v1.size();
	vector<double> v3(size);

	for (int i = 0; i < size; i ++) {
		v3[i] = v1[i] - v2[i];
	}

	return v3;	
}

vector<double> operator*(const double num, const vector<double> &v1) {

	int size = v1.size();
	vector<double> v3(size);

	for (int i = 0; i < size; i ++) {
		v3[i] = v1[i] * num;
	}

	return v3;
}

vector<vector<double> > operator*(const double num, const vector<vector<double> > &v1) {

	int size = v1.size();
	vector<vector<double> > v3(size);

	for (int i = 0; i < size; i ++) {
		v3[i] = num* v1[i];
	}

	return v3;
}

double magnitude(const vector<double> &v) {

	int size = v.size();
	double sum = 0;


	for (int i = 0; i < size; i ++) {

		sum += pow( v[i] , 2 );

	}

	return sqrt( sum );

}


// Makes it easier to print vectors as comma separated values  
string to_string(const vector<double> &v) {

	int size = v.size(); 
	string s = "";

	s = s + to_string(v[0]);

	for (int i = 1; i < size; i ++) {

		s = s + "," + to_string(v[i]);		

	}

	return s;
}

// Single particle. Keeps track of position, momentum, mass, and any other attributes 
struct Particle {

	vector<double> position;
	vector<double> momentum;

	double mass;

	// For use in the Adaptive Stepping Algorithm's "Select" operation
	double prev_force_on;
	bool already_has_timestep = false;


	Particle(double x,double y,double z, double vx, double vy, double vz, double mass_) {

		position.push_back(x);
		position.push_back(y);
		position.push_back(z);


		momentum.push_back(mass_*vx);
		momentum.push_back(mass_*vy);
		momentum.push_back(mass_*vz);

		mass = mass_;

	}

	void change_pos(const vector<double> &delta_pos) {
	
		position = position + delta_pos;
	
	}

	void change_mom(const vector<double> &delta_mom) {
		
		momentum = momentum+delta_mom;

	}

	string get_pos() {

		string s = "";

		s = s + to_string(position);

		return s;

	}

};

struct System {

	vector<Particle> sys;
	double total_mass;

	System() {
		total_mass = 0;
	}

	void add_particle(Particle &p) {

		sys.push_back(p);
		total_mass += p.mass;

	}

	Particle& operator[](int index) {

		return sys[index];

	}

	int size() {
		return sys.size();
	}

	string get_data() {

		string s = "";

		for (int i = 0; i < sys.size(); i++) {

			s = s + sys[i].get_pos() + " , " + to_string(sys[i].mass) + " , ";

		}

		return s;

	}

};



//Define the force function and functions for computing energy

// Force functions can be made more general, however I assume it is a force due to particle interactions
// Force ON p1 DUE TO p2, represented as a vector

//Gravitational Constant

double G = 1;

vector<double> force(const Particle &p1, const Particle &p2) {

	vector<double> f(3); 

	// This is specifically gravity 
	// ( G * m * m / |r|^3 ) * r

	f = ( G*p1.mass*p2.mass / ( pow( magnitude(p2.position - p1.position) , 3 ) ) )  *  (p2.position - p1.position);
	

	return f;

}


// Returns a vector, holding vector<double>, representing the force on each particle in the system.

vector<double> net_force_on_i(System &sys, int i) {
	force_eval_counter += 1;
	int size = sys.size();
	vector<double> net_f(3);

	for (int j = 0; j < size; j ++) {
		if (j != i) {
			net_f = net_f + force(sys[i],sys[j]);
		}
	}

	return net_f;
}

vector<vector<double> > net_force(System &sys) { 
	
	int size = sys.size();

	vector<vector<double> > net_f(size);
	
	for(int i = 0; i < size; i ++) {

		net_f[i] = net_force_on_i(sys,i);

	}

	return net_f;

}

// Functions for computing potential energy of the system, need to rewrite depending on the force

double potential_energy(const Particle &p1, const Particle &p2) {

	double distance = magnitude( (p1.position - p2.position) );
	return ( -G * p1.mass * p2.mass / distance );

}

double total_potential(System &sys) {

	int size = sys.size();
	double total_potential = 0;

	for (int i = 0; i < size; i ++) {

		for (int j = 0; j < i; j ++) {
			
			total_potential += potential_energy(sys[i],sys[j]);
	
		}
	}

	return total_potential;
}

// Functions for computing kinetic energy of the system

double kinetic_energy(Particle &p) {
	
	return pow(magnitude(p.momentum),2)/(2*p.mass);

}

double total_kinetic(System &sys) {

	int size = sys.size();
	double total_kinetic = 0;

	for (int i = 0; i < size; i++) {

		total_kinetic += kinetic_energy(sys[i]);
	
	}

	return total_kinetic;

}

// Angular momentum L = p x r
// Finds the angular momentum of the entire system

vector<double> cross_product(vector<double> &v1, vector<double> &v2) {

	vector<double> v3(3);

	v3[0] = v1[1]*v2[2] - v1[2]*v2[1];
	v3[1] = v1[2]*v2[0] - v1[0]*v2[2];
	v3[2] = v1[0]*v2[1] - v1[1]*v2[0];

	return v3;

}

vector<double> angular_momentum(System &sys) {

	int size = sys.size();

	vector<double> total_angular_momentum(3);

	for (int i = 0; i < size; i ++) {

		total_angular_momentum = total_angular_momentum + cross_product(sys[i].momentum, sys[i].position);

	}

	return total_angular_momentum;
}






// Functions defining different integration methods 



// The driver for the leap-frogging algorithm, where position is updated twice (DKD)
void leap_frog1(System &sys, double dt) {

	vector<vector<double> > f;
	
	int size = sys.size();

	for (int i = 0; i < size; i ++) {

		sys[i].change_pos( ((dt/2)/sys[i].mass) * sys[i].momentum );



	}

	f = net_force(sys);

	for (int i = 0; i < size; i ++) {
		
		sys[i].change_mom( dt * f[i] );
		sys[i].change_pos( ((dt/2)/sys[i].mass) * sys[i].momentum );

	}

	t += dt;

}


// The driver for the leap-frogging algorithm, where momentum is updated twice (KDK)
void leap_frog2(System &sys, double dt) {

	vector<vector<double> > f;
	
	f = net_force(sys);
	
	int size = sys.size();

	for (int i = 0; i < size; i ++) {

		sys[i].change_mom( (dt/2) * f[i] );
		sys[i].change_pos( (dt/sys[i].mass) * sys[i].momentum );

	}

	f = net_force(sys);

	for (int i = 0; i < size; i ++) {

		sys[i].change_mom( (dt/2) * f[i] );

	}

	t += dt;

}



// The driver for my first attempt at an adaptive time-stepping method
// However, this is not true adaptive stepping as it requires the whole
// system to be evolved simultaenously along the same time-step 

void kick1(System &sys, double dt) {

	int size = sys.size();

	vector<vector<double> > f = net_force(sys);

	for(int i = 0; i < size; i ++) {

		sys[i].change_mom( (dt / 2) * f[i] );
	
	}

}

void drift1(System &sys, double dt) {

	int size = sys.size();

	for(int i = 0; i < size; i++) {

		sys[i].change_pos( (dt/sys[i].mass) * sys[i].momentum );

	}

}

double pi = 3.14159;
double n = 0.001;

double select1(System &sys) {

	int size = sys.size();
	double max_distance = 0;
	double rho = 0;
	
	for (int i = 0; i < size; i ++) {

		for (int j = 0; j < i; j ++) {
			
			double dist = magnitude(sys[i].position - sys[j].position);
			
			if (dist > max_distance) {

				max_distance = dist;

			}
		}
	}

	rho = ( sys.total_mass / ( (4/3) * pi * pow((max_distance/2),3) ) );
	return ( n / sqrt(G * rho) );

}

void adaptive_step1(System &sys, double dt) {

	kick1(sys, dt);
	
	double max_dt = select1(sys);
	
	if (dt >= max_dt) {

		kick1(sys, -dt);
		
		adaptive_step1(sys, dt/2);
		//kick(sys, dt);
		adaptive_step1(sys, dt/2);

	}

	else {

		drift1(sys, dt);
		kick1(sys, dt);
		t += dt;

	}
}



// The driver for my second attempt at an adaptive time-stepping method. This method allows two (or more) 
// particles to be simultaenously evolved along different time-steps. 
void kick2(System &sys, int i, double dt) {
	
	vector<double> f = net_force_on_i(sys, i);
	sys[i].change_mom( dt * f );
	sys[i].prev_force_on = magnitude(f);

}

void drift2(System &sys, int i, double dt) {
	
	sys[i].change_pos( (dt/sys[i].mass) * sys[i].momentum );

}

bool already_initialized = false;
void initialize_prev_force_on(System &sys) {

	int size = sys.size();

	for (int i = 0; i < size; i ++) {

		sys[i].prev_force_on = magnitude(net_force_on_i(sys,i));

	}

}

double Tolerance = 0.00001;
ofstream myFile2;

//The method for selecting an appropriate step-size
//Need to assign tmin to the minimum step-size required for the system
vector<bool> select2(System &sys, double dt, double &tmin) {
	
	if (not already_initialized) {

		initialize_prev_force_on(sys);
		already_initialized = true;

	}

	int size = sys.size();
	vector<bool> correct_dt(size);
	double ideal_dt;

	tmin = dt;

	for (int i = 0; i < size; i ++) {
		
		
		ideal_dt = cbrt((6*sys[i].mass*Tolerance / pow(sys[i].prev_force_on,2)));

		if (ideal_dt < tmin) {

			tmin = ideal_dt;

		}

		if ((dt <= ideal_dt) and (not sys[i].already_has_timestep)) {

				myFile2 << "Particle, " << "," << i+1 << ", timestep, " << dt << "," << endl;
				correct_dt[i] = true;
				sys[i].already_has_timestep = true;
		
		}
		
	}
	
	return correct_dt;

}

void adaptive_step2(System &sys, double dt) {
	
	int size = sys.size();
	//This will be re-assigned in the "select" function
	double tmin = dt;
	//Drift all particles forward dt/2
	for (int i = 0; i < size; i ++) {
		drift2(sys,i,dt/2);
	}
	
	//The boolean at each index of this vector tells whether that index particle is on this time-step
	vector<bool> correct_dt = select2(sys,dt,tmin);

	//The base case, when our time-step is sufficiently small for all particles
	if (dt <= tmin) {
		for (int i = 0; i < size; i ++) {
				//Only want to kick the particles on THIS timestep 
				if (correct_dt[i]) {

					kick2(sys,i,dt);
					sys[i].already_has_timestep = false;

				}
			}

		for (int i = 0; i < size; i ++) {
				
			drift2(sys,i,dt/2);

		}
	}
	
	else {
		for (int i = 0; i < size; i ++) {
			
			drift2(sys,i,-dt/2);
		
		}

		adaptive_step2(sys, dt/2);
		
		for (int i = 0; i < size; i ++) {
			
			if (correct_dt[i]) {
				
				kick2(sys, i, dt);


			}
		}

		adaptive_step2(sys, dt/2);

		for (int i = 0; i < size; i ++) {
			if (correct_dt[i]) {
				
				//Need to reset whether a particle has already had a time-step
				sys[i].already_has_timestep = false;

			}
		}
	}
}





//Implementing the RK method described in the IAS15 paper

#include <eigen3/Eigen/Dense>
using namespace Eigen;

//Coefficients for the RK method, where the first index is the index of the particle in question
//and the second index is the dimension (in 3 dimensional space)
vector<vector<Matrix<double,8,1> > > b_values;
vector<vector<Matrix<double,8,1> > > g_values;
vector<vector<double> > init_pos;
vector<vector<double> > init_vel;
vector<vector<double> > init_acc;
Matrix<double,8,8> c_m;
double h1 = 0.0562625605369221464656521910318;
double h2 = 0.180240691736892364987579942780; 
double h3 = 0.352624717113169637373907769648;
double h4 = 0.547153626330555383001448554766;
double h5 = 0.734210177215410531523210605558; 
double h6 = 0.885320946839095768090359771030;
double h7 = 0.977520613561287501891174488626;
double h8 = 1;
//Initializes the position, velocity and acceleration for the beginning of the time-step
//Is done at the start of every time-step
void initialize_RK_step(System &sys) {

	int size = sys.size();
	for (int i = 0; i < size; i ++) {

		init_pos[i] = (sys[i].position);
		init_vel[i] = ((1/sys[i].mass) * sys[i].momentum);
		init_acc[i] = ((1/sys[i].mass) * net_force_on_i(sys, i));

	}	

}


static const double h[8]    = { 0.0, 0.0562625605369221464656521910318, 0.180240691736892364987579942780, 0.352624717113169637373907769648, 0.547153626330555383001448554766, 0.734210177215410531523210605558, 0.885320946839095768090359771030, 0.977520613561287501891174488626};
void find_matrix(int j, int k) {
	if(j > 7 or k > 7) {
		return;
	}

	else {
		if(j == k) {
			c_m(k,j) = 1;
		}
		else {
			if (j > 0 and k == 0) {
				c_m(k,j) = -h[j]*c_m(0,j-1);
			}
			else{
				if(k<j) {

					c_m(k,j) = c_m(k-1,j-1) - h[j]*c_m(k,j-1);
				}
			}
		}

		if(j < 7) {
			find_matrix(j+1,k);
		}
		else {
			find_matrix(k+1,k+1);
		}

	}

}

//Initializes the b_value and g_value vectors with i particles having initial coefficients equal to 0
//Only needs to be done once per simulation
bool very_first_timestep = true;
void initialize_first_RK_step(System &sys) {
	
	if(very_first_timestep){
		
		int size = sys.size();
		
		vector<Matrix<double,8,1> > v_coeff(3);
		Matrix<double,8,1> coeff;
		for (int i = 0; i < 8; i ++) {
			coeff[i] = 0;
		}
		vector<double> V(3);
		for (int j = 0; j < 3; j ++) {
			v_coeff[j] = coeff;
		}

		for (int i = 0; i < size; i ++) {

			b_values.push_back(v_coeff);
			g_values.push_back(v_coeff);
			init_pos.push_back(V);
			init_vel.push_back(V);
			init_acc.push_back(V);
		}

		find_matrix(0,0);

	}

	very_first_timestep = false;		


	initialize_RK_step(sys);
}


//Advances the position and velocity of all particles in the system to t = dt*h
void substepRK(System &sys, double h, double dt) {
	
	int size = sys.size();
	for (int i = 0; i < size; i++) {
		
		vector<double> pos(3);
		vector<double> vel(3);

		for (int j = 0; j < 3; j++) {
			pos[j] = init_pos[i][j]
				+ ((h*dt) * init_vel[i][j]) 
				+ ((h*h*dt*dt/2)*init_acc[i][j]) 
				+ ((h*h*h*dt*dt/6)*b_values[i][j][0]) 
				+ ((h*h*h*h*dt*dt/12)*b_values[i][j][1]) 
				+ ((h*h*h*h*h*dt*dt/20)*b_values[i][j][2])
				+ ((h*h*h*h*h*h*dt*dt/30)*b_values[i][j][3])
				+ ((h*h*h*h*h*h*h*dt*dt/42)*b_values[i][j][4])
				+ ((h*h*h*h*h*h*h*h*dt*dt/56)*b_values[i][j][5])
				+ ((h*h*h*h*h*h*h*h*h*dt*dt/72)*b_values[i][j][6])
				+ ((h*h*h*h*h*h*h*h*h*h*dt*dt/90)*b_values[i][j][7]);



			vel[j] = init_vel[i][j]
				+ ((h*dt)*init_acc[i][j]) 
				+ ((h*h*dt/2)*b_values[i][j][0]) 
				+ ((h*h*h*dt/3)*b_values[i][j][1]) 
				+ ((h*h*h*h*dt/4)*b_values[i][j][2])
				+ ((h*h*h*h*h*dt/5)*b_values[i][j][3])
				+ ((h*h*h*h*h*h*dt/6)*b_values[i][j][4])
				+ ((h*h*h*h*h*h*h*dt/7)*b_values[i][j][5])
				+ ((h*h*h*h*h*h*h*h*dt/8)*b_values[i][j][6])
				+ ((h*h*h*h*h*h*h*h*h*dt/9)*b_values[i][j][7]);
		}

		sys[i].position = pos;
		sys[i].momentum = sys[i].mass * vel;

	}
}

//Reassigns the b_values based upon the current g_values
void convert_g_to_b(System &sys) {

	int size = sys.size();			
	
	Matrix<double,8,1> new_b;
	
	for(int i = 0; i < size; i ++) {

		for(int j = 0; j < 3; j ++) {
			
			new_b = c_m * g_values[i][j];
			b_values[i][j] = new_b;
			
		}
	}


}

//Finds the g_values based on the initial conditions and the current b_values
void find_g_values(System &sys, double dt) {

	int size = sys.size();

	//Steps system forward to dt = h1;
	substepRK(sys, h1, dt);

	//Calculates all g1 values
	for (int i = 0; i < size; i ++) {

		vector<double> next_acc;
		next_acc = (1/sys[i].mass) * net_force_on_i(sys, i);

		for (int j = 0; j < 3; j ++) {

			g_values[i][j][0] = ( (next_acc[j] - init_acc[i][j]) / h1 );

		}

	}

	//Steps system forward to dt = h2;
	substepRK(sys, h2, dt);

	//Calculates all g2 values
	for (int i = 0; i < size; i ++) {

		vector<double> next_acc;
		next_acc = (1/sys[i].mass) * net_force_on_i(sys, i);

		for (int j = 0; j < 3; j ++) {

			g_values[i][j][1] = ( (next_acc[j] - init_acc[i][j] - g_values[i][j][0]*h2 ) / (h2*(h2-h1)) );

		}

	}

	//Steps system forward to dt = h3;
	substepRK(sys, h3, dt);

	//Calculates all g3 values
	for (int i = 0; i < size; i ++) {

		vector<double> next_acc;
		next_acc = (1/sys[i].mass) * net_force_on_i(sys, i);

		for (int j = 0; j < 3; j ++) {

			g_values[i][j][2] = ( (next_acc[j] - init_acc[i][j] - g_values[i][j][0]*h3 - g_values[i][j][1]*h3*(h3-h1) ) / (h3*(h3-h1)*(h3-h2)) );

		}

	}

	//Steps system forward to dt = h4
	substepRK(sys, h4, dt);

	//Calculates all g4 values
	for (int i = 0; i < size; i ++) {

		vector<double> next_acc;
		next_acc = (1/sys[i].mass) * net_force_on_i(sys, i);

		for (int j = 0; j < 3; j ++) {

			g_values[i][j][3] = ( (next_acc[j] - init_acc[i][j] - g_values[i][j][0]*h4 - g_values[i][j][1]*h4*(h4-h1) - g_values[i][j][2]*h4*(h4-h1)*(h4-h2) ) / (h4*(h4-h1)*(h4-h2)*(h4-h3)) );

		}

	}

	//Steps system forward to dt = h5
	substepRK(sys, h5, dt);

	//Calculates all g5 values
	for (int i = 0; i < size; i ++) {

		vector<double> next_acc;
		next_acc = (1/sys[i].mass) * net_force_on_i(sys, i);

		for (int j = 0; j < 3; j ++) {

			g_values[i][j][4] = ( (next_acc[j] - init_acc[i][j] - g_values[i][j][0]*h5 - g_values[i][j][1]*h5*(h5-h1) - g_values[i][j][2]*h5*(h5-h1)*(h5-h2) - g_values[i][j][3]*h5*(h5-h1)*(h5-h2)*(h5-h3) ) / (h5*(h5-h1)*(h5-h2)*(h5-h3)*(h5-h4)) );

		}

	}

	//Steps system forward to dt = h6
	substepRK(sys, h6, dt);

	//Calculates all g6 values
	for (int i = 0; i < size; i ++) {

		vector<double> next_acc;
		next_acc = (1/sys[i].mass) * net_force_on_i(sys, i);

		for (int j = 0; j < 3; j ++) {

			g_values[i][j][5] = ( (next_acc[j] - init_acc[i][j] - g_values[i][j][0]*h6 - g_values[i][j][1]*h6*(h6-h1) - g_values[i][j][2]*h6*(h6-h1)*(h6-h2) - g_values[i][j][3]*h6*(h6-h1)*(h6-h2)*(h6-h3) - g_values[i][j][4]*h6*(h6-h1)*(h6-h2)*(h6-h3)*(h6-h4) ) / (h6*(h6-h1)*(h6-h2)*(h6-h3)*(h6-h4)*(h6-h5)) );

		}

	}

	//Steps system forward to dt=h7
	substepRK(sys, h7, dt);

	//Calculates all g7 values
	for (int i = 0; i < size; i ++) {

		vector<double> next_acc;
		next_acc = (1/sys[i].mass) * net_force_on_i(sys, i);

		for (int j = 0; j < 3; j ++) {

			g_values[i][j][6] = ( (next_acc[j] - init_acc[i][j] - g_values[i][j][0]*h7 - g_values[i][j][1]*h7*(h7-h1) - g_values[i][j][2]*h7*(h7-h1)*(h7-h2) - g_values[i][j][3]*h7*(h7-h1)*(h7-h2)*(h7-h3) - g_values[i][j][4]*h7*(h7-h1)*(h7-h2)*(h7-h3)*(h7-h4) - g_values[i][j][5]*h7*(h7-h1)*(h7-h2)*(h7-h3)*(h7-h4)*(h7-h5) ) / (h7*(h7-h1)*(h7-h2)*(h7-h3)*(h7-h4)*(h7-h5)*(h7-h6)) );

		}

	}

	//Steps system forward to dt = dt
	substepRK(sys, h8, dt);

	//Calculates all g8 values
	for (int i = 0; i < size; i ++) {

		vector<double> next_acc;
		next_acc = (1/sys[i].mass) * net_force_on_i(sys, i);

		for (int j = 0; j < 3; j ++) {

			g_values[i][j][7] = ( (next_acc[j] - init_acc[i][j] - g_values[i][j][0]*h8 - g_values[i][j][1]*h8*(h8-h1) - g_values[i][j][2]*h8*(h8-h1)*(h8-h2) - g_values[i][j][3]*h8*(h8-h1)*(h8-h2)*(h8-h3) - g_values[i][j][4]*h8*(h8-h1)*(h8-h2)*(h8-h3)*(h8-h4) - g_values[i][j][5]*h8*(h8-h1)*(h8-h2)*(h8-h3)*(h8-h4)*(h8-h5) - g_values[i][j][6]*h8*(h8-h1)*(h8-h2)*(h8-h3)*(h8-h4)*(h8-h5)*(h8-h6) ) / (h8*(h8-h1)*(h8-h2)*(h8-h3)*(h8-h4)*(h8-h5)*(h8-h6)*(h8-h7)) );

		}

	}

	//Returns system to dt = 0 and converts the g values to b values
	
	substepRK(sys,0,dt);
	convert_g_to_b(sys);
	
}


//The driver for one time-step of the RK algorithm
//Recursively updates b and g values until they approach machine precision
bool first_step = true;

vector<vector<Matrix<double,8,1> > > old_b_values;

int count_ = 0;
double epsilon = 0.00000028;

void RK_Driver(System &sys, bool use_adaptive_stepsize) {

	//Needs to be done once per time-step in order get initial values
	if (first_step) {
		initialize_first_RK_step(sys);
		find_g_values(sys,dt);
		first_step = false;
	}

	
	int size = sys.size();

	old_b_values = b_values;
	
	find_g_values(sys,dt);

	//Determines whether the g_values have converged to machine precision yet
	double max_del_b_6 = 0;
	double max_y_pp = 0;

	//For use in Adaptive Stepsize Algorithm
	double max_b_6 = 0;

	for (int i = 0; i < size; i ++) {
		for (int j = 0; j < 3; j ++) {
			max_b_6 = max( abs(max_b_6) , abs(b_values[i][j][5]) );
			max_del_b_6 = max( abs(max_del_b_6) , abs(old_b_values[i][j][5] - b_values[i][j][5]) );
			max_y_pp = max( abs(max_y_pp) , abs(init_acc[i][j]) );
		
		}
	}

	double global_error = max_del_b_6 / max_y_pp ;
	//Determines if the predictor corrector loop has converged or not
	if ( ( global_error > pow(10,-16) ) and (count_ < 12)) {



		count_ ++;
		RK_Driver(sys, use_adaptive_stepsize);

	}

	else {

		if (use_adaptive_stepsize) {
			double dt_req = dt * pow ( (epsilon / (max_b_6 / max_y_pp)) , 0.14285714);
			cout << dt_req << endl;
			if (dt > dt_req) {
				dt = dt_req;
				RK_Driver(sys, use_adaptive_stepsize);
			}

			else {

				t += dt;
				substepRK(sys, 1, dt);
				dt = dt_req;
				first_step = true;

			}

		}



		t += dt;
		substepRK(sys, 1, dt);
		first_step = true;

	}


}





// Output Energy / Position Data to csv file

void output_energy(System sys, int end_time, string name, double _dt, string method) {
	
	ofstream myFile;
	myFile.open(name+".csv");

	double KE;
	double PE;

	t = 0;
	dt = _dt;
	double initial_energy = total_kinetic(sys) + total_potential(sys);

	while (t < end_time) {

		KE = total_kinetic(sys);
		PE = total_potential(sys);

		myFile << "Potential Energy , " << PE << " , "
			   << "Kinetic Energy , " << KE << " , "
			   << "Total Energy , " << KE + PE
			   << ", Delta Energy , " << initial_energy - (KE + PE) 
			   << ", Time , " << t 
			   << ", Angular Momentum , " << to_string(angular_momentum(sys)) <<  endl;

		//integration method
		if (method == "AT1") {
			adaptive_step1(sys, dt);
			t += dt;
		}

		if (method == "AT2") {
			adaptive_step2(sys, dt);
			t += dt;
		}

		if (method == "LF1") {
			leap_frog1(sys,dt);
		}

		if (method == "RK") {
			RK_Driver(sys,false);
			count_ = 0;
		}

		if (method == "RK_AT") {
			RK_Driver(sys,true);
			count_ = 0;
		}

	}
}


void output_position(System sys, int end_time, string name, double _dt, string method) {

	ofstream myFile;
	myFile.open(name+".csv");

	t = 0;
	dt = _dt;
	while (t < end_time) {
		
		//Only records every integer time-step
		myFile << sys.get_data() << " Time , " << t << "," << endl;

		//integration method
		if (method == "AT1") {
			adaptive_step1(sys, dt);
			t += dt;
		}

		if (method == "AT2") {
			adaptive_step2(sys, dt);
			t += dt;
		}

		if (method == "LF1") {
			leap_frog1(sys,dt);
		}

		if (method == "LF2") {
			leap_frog2(sys,dt);
		}

		if (method == "RK") {
			RK_Driver(sys,false);
			count_ = 0;
		}

		if (method == "RK_AT") {
			RK_Driver(sys,true);
			count_ = 0;
		}

	}
}




int main() {

	System sys;
	
	Particle p1(400,0,0,0,0.5,0,300);
	Particle p2(230,0,0,0,2,0,25);
	Particle p3(248,0,0,0,0.6,0,1);
	Particle p4(-250,0,0,0,-0.5,0,300);
	Particle p5(-400,0,0,0,1.1,0,32);
	Particle p6(-430,0,0,0,0,0,2);
	
	sys.add_particle(p1);
	sys.add_particle(p2);
	sys.add_particle(p3);
	sys.add_particle(p4);
	sys.add_particle(p5);
	sys.add_particle(p6);


	auto start1 = chrono::steady_clock::now();
	
	output_energy(sys, 50, "Data/Output/Energy", 0.5, "RK_AT");
	output_position(sys, 50, "Data/Output/Data", 0.5, "RK_AT");
	cout << "Force Evals with RK: " << force_eval_counter << endl;
	
	auto end1 = chrono::steady_clock::now();


	auto start2 = chrono::steady_clock::now();
	
	force_eval_counter = 0;
	//output_position(sys, 250, "Data/Output/Data2", 0.5, "RK_AT");
	//output_energy(sys, 250, "Data/Output/Energy2", 0.5, "RK_AT");
	cout << "Force Evals with LF: " << force_eval_counter << endl;
	
	auto end2 = chrono::steady_clock::now();

	cout << "Runga Kutta: " << chrono::duration <double, milli> (end1-start1).count() << " ms" << endl;
	cout << "LF: " << chrono::duration <double, milli> (end2-start2).count() << " ms" << endl;

	return 0;
}