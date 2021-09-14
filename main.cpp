#include <iostream>
#include <vector>
#include <tuple>
#include "magnetometer.h"
#include <cmath>
#include <armadillo>
#include <complex>
#include <iterator>
#include <ctime>
#include <filesystem>
#include <iomanip>
#include <cstdlib>

using namespace std::complex_literals;

int main(int argc, char* argv[]){
	if(argc > 4){
		std::cerr << "Too many arguments passed to main.\nFirst argument should be the table of variables.\nThe second argument, which is optional, is the name of the file with results with to generate.\nThe third argument, also optional, is the seed for armadillos RNG to use." << std::endl;
		return -1;
	}
	if(argc == 1){
		std::cerr << "As argument 1, main must have the predefined table of variables.\n" << std::endl;
		return -1;
	}

	std::string line;
	std::fstream input_variables(argv[1]);
	std::vector<double> loaded_variables;
	if (input_variables.is_open()) {
		while (std::getline(input_variables, line))
		{
			if (line[0] != '#' && line[0] != '\n'){
				std::istringstream iss(line);
				double num;

				while ((iss >> num)) loaded_variables.push_back(num);
			}
		}
	}
	input_variables.close();

	//Setting parameters
	unsigned int load_index = 0;
	const double N = loaded_variables[load_index++];
	const double T = loaded_variables[load_index++];
	const double base_rate = loaded_variables[load_index++];
	const int frames = loaded_variables[load_index++];
	const double beta = loaded_variables[load_index++];
	const double g = loaded_variables[load_index++];
	const double kappa = loaded_variables[load_index++];
	const double kappa_1 = loaded_variables[load_index++];
	const double gamma_dec = loaded_variables[load_index++];
	const double gamma_phi = loaded_variables[load_index++];
	const double LOPhi = M_PI / 2.;
	const double delta_c = 0;
	std::vector<double> delta_q(N);
	for (int i = 0; i<N; ++i) delta_q.at(i) = 4.0 * i / (N - 1) - 2.0;
	const int emission_startpoint = N / 2;
	const double eta = 1.;
	const double dt = 1./200.;
	if (frames > (T/dt)){
		std::cerr << "TOO MANY FRAMES, TOO LITTLE TIME" << std::endl;	
		return -1;
	}
	std::srand(static_cast<unsigned int>(time(NULL)));
	int seed;
       	if(argc < 4) seed = std::rand();
	if(argc == 4) seed = atoi(argv[3]);
	arma::arma_rng::set_seed(seed);
	
	
	// Main body
	const auto [r_up, r_down] = ehrenfest_chain(N, base_rate);
	const auto [c_1, c_2, c_3, c, H, J_up, J_down] = make_operators(N, g, kappa, kappa_1, delta_c, delta_q, beta, gamma_dec, gamma_phi, LOPhi, r_up, r_down);
	const auto [dY, true_state, t_emission] = homodyne_emission(N, H, J_up, J_down, c, c_1, c_2, c_3, eta, emission_startpoint, delta_q, dt, T);
	const auto [P, P_f, P_b, P_max, P_f_max, t_PQS] = homodyne_forward_backward_conditioning(T, dt, N, H, J_up, J_down, c, c_1, c_2, c_3, eta, dY, delta_q, frames);

	

	// Saving data
	std::cout << std::fixed;
	std::cout << std::setprecision(10);


	//Saving input parameters. If no resulting folder name is given, time and date is the default.
	
	std::string buffer;

	if(argc == 2){
	time_t curr_time;
	tm *curr_tm;
	time(&curr_time);
	curr_tm = localtime(&curr_time);

	char time_buffer [80];
	strftime(time_buffer, 80, "%d_%m_%Y__%H_%M_%S", curr_tm);
	buffer = time_buffer;
	}



	if(argc == 3 || argc == 4) buffer = argv[2];



	std::string run_folder = "./";
	run_folder.append(buffer);
	std::string input_variables_folder = run_folder + "/input_parameters";
	std::filesystem::create_directories(input_variables_folder);
	std::vector<std::string> input_variable_names{"N", "T", "base_rate", "frames", "beta", "g", "kappa", "kappa_1", "gamma_dec", "gamma_phi", "dt", "seed"};
	for (auto &c : input_variable_names) c.append(".txt");
	const double input_variables_array[] = {N, T, base_rate, static_cast<double>(frames), beta, g, kappa, kappa_1, gamma_dec, gamma_phi, dt, static_cast<double>(seed)};

	for (size_t i = 0; i < input_variable_names.size(); ++i){
		std::ofstream save_input_variables(input_variables_folder + "/" + input_variable_names.at(i));
		save_input_variables << std::fixed;
		save_input_variables << std::setprecision(10);
		save_input_variables << input_variables_array[i];
		save_input_variables.close();
	}
	std::ofstream delta_q_file(input_variables_folder + "/delta_q.txt");
       	for (auto &c : delta_q) delta_q_file << c << std::endl;
	delta_q_file.close();	

	// Saving generated data
	std::string data_folder = run_folder + "/data/";
	std::filesystem::create_directories(data_folder);
	std::vector<std::string> data_list_names{"dY", "t_emission", "t_PQS"};
	std::vector<std::vector<double>> data_lists{dY, t_emission, t_PQS}; 
	for (size_t i = 0; i < data_list_names.size(); ++i){
		data_list_names.at(i).append(".txt");
		std::ofstream save_data_variables(data_folder + data_list_names.at(i));
		save_data_variables << std::fixed << std::setprecision(3);
		for (auto &c : data_lists.at(i)) save_data_variables << c << std::endl;
		save_data_variables.close();
	}

	std::ofstream true_state_file(data_folder + "true_state.txt");
	for (auto &c : true_state) true_state_file << c << std::endl;
	true_state_file.close();
	
	std::ofstream P_file(data_folder + "P.txt");
	P.save(P_file, arma::raw_ascii);
	P_file.close();

	std::ofstream P_f_file(data_folder + "P_f.txt");
	P_f.save(P_f_file, arma::raw_ascii);
	P_f_file.close();

	std::ofstream P_b_file(data_folder + "P_b.txt");
	P_b.save(P_b_file, arma::raw_ascii);
	P_b_file.close();

	std::ofstream P_max_file(data_folder + "P_max.txt");
	P_max.save(P_max_file, arma::raw_ascii);
	P_max_file.close();

	std::ofstream P_f_max_file(data_folder + "P_f_max.txt");
	P_f_max.save(P_f_max_file, arma::raw_ascii);
	P_f_max_file.close();
//*/	
	return 0;
}
