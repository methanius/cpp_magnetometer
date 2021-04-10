#include <iostream>
#include <armadillo>
#include <tuple>
#include "choi_jamiolkowski.h"
#include <vector>
#include <iterator>
#include <algorithm>


using namespace std::complex_literals;

std::tuple < arma::mat, arma::mat, arma::mat, arma::urowvec, arma::urowvec, std::vector<double>> homodyne_forward_backward_conditioning(const double T, const double dt, const double N, const arma::sp_cx_mat H, const arma::sp_cx_mat J_up, const arma::sp_cx_mat J_down, const arma::sp_cx_mat c, const arma::sp_cx_mat c_1, const arma::sp_cx_mat c_2, const arma::sp_cx_mat c_3, const double eta, const std::vector<double> dY, const std::vector<double> delta_q, const int frames){

        // Initialize rho and E
        arma::sp_cx_mat rho(2 * N, 2 * N), E(2 * N, 2 * N);
        rho.speye();
        E.speye();

        // Liouvillian exponent
        arma::cx_mat L( -1i * (sp_cx_spre(H) - sp_cx_spost(H)) \
                        + sp_cx_sprepost(J_down, J_down.t()) - (sp_cx_spre(J_down.t() * J_down) + sp_cx_spost(J_down.t() * J_down)) / 2.  \
                        + sp_cx_sprepost(J_up, J_up.t()) - (sp_cx_spre(J_up.t() * J_up) + sp_cx_spost(J_up.t() * J_up)) / 2. \
                        + sp_cx_sprepost(c_1, c_1.t()) - (sp_cx_spre(c_1.t() * c_1) + sp_cx_spost(c_1.t() * c_1)) / 2.  \
                        + sp_cx_sprepost(c_2, c_2.t()) - (sp_cx_spre(c_2.t() * c_2) + sp_cx_spost(c_2.t() * c_2)) / 2.  \
                        + sp_cx_sprepost(c_3, c_3.t()) - (sp_cx_spre(c_3.t() * c_3) + sp_cx_spost(c_3.t() * c_3)) / 2.);
        arma::sp_cx_mat eLdt(arma::expmat(L * dt));
	eLdt.clean(arma::datum::eps);

        // Adjoint Liouvillian exponent
        arma::cx_mat EL( 1i * (sp_cx_spre(H) - sp_cx_spost(H)) \
                        + sp_cx_sprepost(J_down.t(), J_down) - (sp_cx_spre(J_down.t() * J_down) + sp_cx_spost(J_down.t() * J_down)) / 2.  \
                        + sp_cx_sprepost(J_up.t(), J_up) - (sp_cx_spre(J_up.t() * J_up) + sp_cx_spost(J_up.t() * J_up)) / 2.  \
                        + sp_cx_sprepost(c_1.t(), c_1) - (sp_cx_spre(c_1.t() * c_1) + sp_cx_spost(c_1.t() * c_1)) / 2.  \
                        + sp_cx_sprepost(c_2.t(), c_2) - (sp_cx_spre(c_2.t() * c_2) + sp_cx_spost(c_2.t() * c_2)) / 2.  \
                        + sp_cx_sprepost(c_3.t(), c_3) - (sp_cx_spre(c_3.t() * c_3) + sp_cx_spost(c_3.t() * c_3)) / 2.);
        arma::sp_cx_mat eELdt(arma::expmat(EL * dt));
	eELdt.clean(arma::datum::eps);

// Propagation of rho and E
        arma::sp_cx_mat c_T(c.t());
        std::vector< arma::sp_cx_mat > rho_saved, E_saved;
        rho_saved.reserve(frames);
        E_saved.reserve(frames);
	double t{0};
	std::vector<double> t_PQS;
        arma::uvec record_times = arma::round(arma::linspace<arma::uvec>(0, dY.size() - 1, frames));
	unsigned int next_record_time = 0;


	// rho
	for(size_t n = 0; n < dY.size(); ++n){
                rho = sp_cx_channel_to_state( eLdt * sp_cx_state_to_channel(rho)) + (c * rho + rho * c_T) * dY.at(n) * sqrt(eta);
                rho /= arma::trace(rho);
		//std::cout << n << " == " << record_times.at(next_record_time) << " = " << (n == record_times.at(next_record_time)) << std::endl;
                if (n == record_times.at(next_record_time)){
                        rho_saved.push_back(rho);
			t_PQS.push_back(t);
                        ++next_record_time;
                }
		t += dt;
        }


        //E propagation and saving
        next_record_time = frames - 1;
        unsigned int n = dY.size();
        do{
                --n;
                E = sp_cx_channel_to_state( eELdt * sp_cx_state_to_channel(E)) + (c_T * E + E * c) * dY.at(n) * sqrt(eta);
                E /= arma::trace(E);
                if (n == record_times.at(next_record_time)){
                        E_saved.push_back(E);
                        next_record_time--;
                }
        } while(n>0);
	

	std::reverse(E_saved.begin(), E_saved.end());


	// Filtered estimate
	arma::mat P_f(N, frames);
	
	for (int j = 0; j < frames; ++j){
                for (int i = 0; i < N; ++i){
                        P_f.at(i, j) = rho_saved.at(j).at(2 * i, 2 * i).real() + rho_saved.at(j).at(2 * i + 1, 2 * i + 1).real();
                }
        }
	
	
        // PQS estimate
        arma::mat P(N, frames);
	
        for (int j = 0; j < frames; ++j){
                for (int i = 0; i < N; ++i){
                        P.at(i, j) = rho_saved.at(j).at(2 * i, 2 * i).real() * E_saved.at(j).at(2 * i, 2 * i).real() \
                        + rho_saved.at(j).at(2 * i + 1, 2 * i + 1).real() * E_saved.at(j).at(2 * i + 1, 2 * i + 1).real() \
                        + rho_saved.at(j).at(2 * i, 2 * i + 1).real() * E_saved.at(j).at(2 * i + 1, 2 * i).real()\
                        + rho_saved.at(j).at(2 * i + 1, 2 * i).real() * E_saved.at(j).at(2 * i, 2 * i + 1).real();
                }
        }
	
	arma::Row P_norms = arma::sum(P, 0);
	for (int j = 0; j < frames; ++j){
		for (int i = 0; i<N; ++i){
			P.at(i, j) /= P_norms(j);
		}
	}

	arma::mat P_b(N, frames);
	for (int j = 0; j < frames; ++j){
		for (int i = 0; i < N; ++i){
			P_b.at(i, j) = E_saved.at(j).at(2 * i, 2 * i).real() + E_saved.at(j).at(2 * i + 1, 2 * i + 1).real();
		}
	}

	arma::urowvec P_max = arma::index_max(P, 0);
	arma::urowvec P_f_max = arma::index_max(P_f, 0);






        return {P, P_f, P_b, P_max, P_f_max, t_PQS};
}

