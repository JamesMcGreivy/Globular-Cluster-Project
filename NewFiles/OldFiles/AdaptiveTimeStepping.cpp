#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
#include <random>
using namespace std;

double t = 0;
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
	
	cout << dt << " max: " << select1(sys) << " time: " << t << endl;
	
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

















// Output Energy / Position Data to csv file

void output_energy(System sys, int num_iterations, string name, double dt, string method) {
	
	ofstream myFile;
	myFile.open(name+".csv");

	double KE;
	double PE;

	t = 0;

	double initial_energy = total_kinetic(sys) + total_potential(sys);

	for (int i = 0; i < num_iterations; i++) {

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
			adaptive_step1(sys, dt);
			t += dt;
		}

		if (method == "LF1") {
			leap_frog1(sys,dt);
		}

		if (method == "LF2") {
			leap_frog2(sys,dt);
		}

	}
}


void output_position(System sys, int num_iterations, string name, double dt, string method) {

	ofstream myFile;
	myFile.open(name+".csv");

	t = 0;

	for (int i = 0; i < num_iterations; i++) {
		if (2*t == int(2*t)) {
			myFile << sys.get_data() << " Time , " << t << "," << endl;
		}
		
		//integration method
		if (method == "AT1") {
			adaptive_step1(sys, dt);
			t += dt;
		}

		if (method == "AT2") {
			adaptive_step1(sys, dt);
			t += dt;
		}

		if (method == "LF1") {
			leap_frog1(sys,dt);
		}

		if (method == "LF2") {
			leap_frog2(sys,dt);
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
	/*
	
	Particle p1(200,0,0,0,.2,0,100);
	Particle p2(-200,0,0,0,-.2,0,100);
	sys.add_particle(p1);
	sys.add_particle(p2);
	*/
	auto start1 = chrono::steady_clock::now();

	output_energy(sys, 2400, "Data/Output/EnergyLF", 0.5, "LF");
	output_position(sys, 2400, "Data/Output/DataLF", 0.5, "LF");

	cout << "Force Evals with LF: " << force_eval_counter << endl;
	auto end1 = chrono::steady_clock::now();

	auto start2 = chrono::steady_clock::now();
	
	myFile2.open("Data/Output/timestep_eval_data.csv");
	
	
	force_eval_counter = 0;
	
	output_energy(sys, 2400, "Data/Output/EnergyAT", 0.5, "AT");

	already_initialized = false;
	output_position(sys, 10000, "Data/Output/DataAT", 0.5, "AT");

	cout << "Force Evals with AT: " << force_eval_counter << endl;

	auto end2 = chrono::steady_clock::now();

	cout << "Leap Frog: " << chrono::duration <double, milli> (end1-start1).count() << " ms" << endl;

	cout << "Adaptive Time: " << chrono::duration <double, milli> (end2-start2).count() << " ms" << endl;

	return 0;
}

