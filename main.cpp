#include <eigen3/Eigen/Dense>
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <vector>

void output_data(std::vector<Data> data, std::string file_name)
{
	std::ofstream position_data;
	position_data.open(file_name);

    position_data   << "PosX" << "," << "PosY" << "," << "PosZ" << "," 
                    << "VelX" << "," << "VelY" << "," << "VelZ" << "," << "time" << std::endl;
	
	for (int i = 0; i < data.size(); i += 1)
	{	
		
		for (int component = 0; component < 3; component ++)
		{
				position_data << std::setprecision(17) << data[i][0](component) << ",";
		}

        for (int component = 0; component < 3; component ++)
        {
                position_data << std::setprecision(17) << data[i][1](component) << ",";
        }

		position_data << data[i].time << std::setprecision(17) << std::endl;
	}

}

int main(int argc, char *argv[])
{
    std::vector<Data> data = Integrator::Integrate(E, B, tInit, tFinal, posInit, velInit, 0, 
        use_adaptive_timestep, std::stod(argv[1]));
    output_data(data, argv[2]);
    return 0;
}
