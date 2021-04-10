#include <iostream>
#include <vector>
#include <tuple>

//auto ehrenfest_chain(const int size, const double base_rate){
std::tuple< std::vector<double>, std::vector<double> > ehrenfest_chain(const int size, const double base_rate){

	std::vector<double> r_up(size - 1, 0); 

	double i = 0;
	for (auto &c : r_up) c = (i++ + 1) / (size - 1) * base_rate;

	const std::vector<double> r_down(r_up.rbegin(), r_up.rend());

	return {r_up, r_down};
}
