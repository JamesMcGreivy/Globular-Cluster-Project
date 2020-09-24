#include <eigen3/Eigen/Dense>
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <vector>
#include "NewFiles/integrator.cpp"
#include <chrono>

//Electric Field function to Use
#include "test3/FieldFunctions.cpp"

std::string file_name;
void output_data(std::vector<Data> data)
{

	double num_particles = data[0][0].rows();

	std::ofstream position_data;
	position_data.open("Output/pos_"+file_name+".csv");
	
	std::ofstream velocity_data;
	velocity_data.open("Output/vel_"+file_name+".csv");

	for (int i = 0; i < data.size(); i += 1)
	{	
		
		for (int particle = 0; particle < num_particles; particle ++)
		{
			for (int component = 0; component < 3; component ++)
			{
				position_data << std::setprecision(16) << data[i][0](particle,component) << ",";
				velocity_data << std::setprecision(16) << data[i][1](particle,component) << ",";
			}
		}

		position_data << data[i].time << std::setprecision(16) << std::endl;
		velocity_data << data[i].time << std::setprecision(16) << std::endl;
	}

}

std::string config_path;
std::vector<double> configs;
void read_config()
{
    std::ifstream cFile (config_path);
    if (cFile.is_open())
    {
        std::string line;
        getline(cFile,line);
        for(int i = 0; i <= 7; i++)
        {
        	getline(cFile,line);
            line.erase(std::remove_if(line.begin(), line.end(), isspace),
                                 line.end());
            if(line[0] == '#' || line.empty())
                continue;

            auto delimiterPos = line.find("=");
            auto value = line.substr(delimiterPos + 1);
            if (i <= 1)
            {	
            	configs.push_back(atof(&value[1]));
            	
            	delimiterPos = value.find(",");
            	value = value.substr(delimiterPos + 1);
            	
            	configs.push_back(atof(&value[0]));

            	delimiterPos = value.find(",");
            	value = value.substr(delimiterPos + 1);

            	configs.push_back(atof(&value[0]));
            }
            else
            {
            	if (i == 7)
            	{
            		file_name = value;
            	}
            	else
            	{
            		configs.push_back( stof(value) );
	           	}
            }
        }
        
    }
    else {
        std::cerr << "Couldn't open config file for reading.\n";
    }
}


int main(int argc, char *argv[])
{
	Eigen::Matrix<double,1,3> pos;
	Eigen::Matrix<double,1,3> vel;

    //Reads in all of the command line arguments from the makefile
	config_path = argv[1];
    read_config();
    pos << 	configs[0],configs[1],configs[2];
	vel <<	configs[3],configs[4],configs[5];
	bool predict_b = configs[8];
	bool use_adaptive_timestep = configs[9];
	double adaptive_factor = configs[10];
 

	double t_init = 0;
	double dt = configs[6];
	double t_final = configs[7];

	std::vector<Data> data1 = LorentzForceIntegrator::Integrate(E, B, t_init, t_final, pos, vel, dt, predict_b, use_adaptive_timestep, adaptive_factor);
    //std::vector<Data> data1 = BorisIntegrator2::Integrate(E, B, t_init, t_final, pos, vel, dt, predict_b, use_adaptive_timestep, adaptive_factor);
	output_data(data1);

	return 0;

}
