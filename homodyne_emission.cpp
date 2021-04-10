#include <vector>
#include <tuple>
#include <cmath>
#include <armadillo>
#include <complex>
#include "choi_jamiolkowski.h"
#include <iostream>

using namespace std::complex_literals;

std::tuple< std::vector<double>, std::vector<int>, std::vector<double> > homodyne_emission(const double N, const arma::sp_cx_mat H, const arma::sp_cx_mat J_up, const arma::sp_cx_mat J_down, const arma::sp_cx_mat c, const arma::sp_cx_mat c_1, const arma::sp_cx_mat c_2, const arma::sp_cx_mat c_3, const double eta, const int emission_startpoint, const std::vector<double> delta_q, const double dt, const double T){

        // Initialize rho to 2x2 identity in initial sub-rho
        arma::sp_cx_mat rho(2 * N, 2 * N);
        rho.submat(2 * emission_startpoint, 2 * emission_startpoint, 2 * emission_startpoint + 1, 2 * emission_startpoint + 1).eye();

        // Liouvillian
        arma::cx_mat L(-1i * (sp_cx_spre(H) - sp_cx_spost(H)) \
                         - (sp_cx_spre(J_up.t() * J_up) + sp_cx_spost(J_up.t() * J_up)) / 2 \
                         - (sp_cx_spre(J_down.t() * J_down) + sp_cx_spost(J_down.t() * J_down)) / 2 \
                         + sp_cx_sprepost(c_1, c_1.t()) - (sp_cx_spre(c_1.t() * c_1) + sp_cx_spost(c_1.t() * c_1)) / 2 \
                         + sp_cx_sprepost(c_2, c_2.t()) - (sp_cx_spre(c_2.t() * c_2) + sp_cx_spost(c_2.t() * c_2)) / 2 \
                         + sp_cx_sprepost(c_3, c_3.t()) - (sp_cx_spre(c_3.t() * c_3) + sp_cx_spost(c_3.t() * c_3)) / 2);
	

	arma::sp_cx_mat eLdt(arma::expmat(L * dt));
	eLdt.clean(arma::datum::eps);
		

        // Signal generation
        int state_track = emission_startpoint;
        std::vector< int > true_state;
	true_state.reserve(T / dt);
        std::vector<double> dY;
	dY.reserve(T / dt);
	auto J_up_T = J_up.t();
	auto J_down_T = J_down.t();
	auto c_T = c.t();
	std::vector<double> t_emission;
	t_emission.reserve(T / dt);

        for(double t = 0; t<T; t += dt){
                auto first_jump_check = arma::randu();
                auto test_permutation = arma::randu();
                if (test_permutation < 0.5){
                        if (first_jump_check < (arma::trace(J_up * rho * J_up_T).real() * dt)){
                                rho = J_up * rho * J_up_T;
                                rho /= arma::trace(rho);
                                --state_track;
                        }
                        else{
                                auto second_jump_check = arma::randu();
                                if (second_jump_check < (arma::trace(J_down * rho * J_down_T).real() * dt)){
                                        rho = J_down * rho * J_down_T;
                                        rho /= arma::trace(rho);
                                        ++state_track;
                                }
                        }
                }
                else{
                        if (first_jump_check < (arma::trace(J_down * rho * J_down_T).real() * dt)){
                                rho = J_down * rho * J_down_T;
                                rho /= arma::trace(rho);
                                ++state_track;
                        }
                        else{
                                auto second_jump_check = arma::randu();
                                if (second_jump_check < (arma::trace(J_up * rho * J_up_T).real() * dt)){
                                        rho = J_up * rho * J_up_T;
                                        rho /= arma::trace(rho);
                                        --state_track;
                                }
                        }
                }
                double dW = arma::randn() * sqrt(dt);
                double dY_increment = (arma::trace(c * rho + rho * c_T) * sqrt(eta) * dt + dW).real();
                dY.push_back(dY_increment);
                rho = sp_cx_channel_to_state(eLdt * sp_cx_state_to_channel(rho)) + (c * rho + rho * c_T) * dY_increment * sqrt(eta);
                rho /= arma::trace(rho);
                true_state.push_back(state_track);
		t_emission.push_back(t);
        }
return {dY, true_state, t_emission};
}

