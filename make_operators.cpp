#include <vector>
#include <tuple>
#include <cmath>
#include <armadillo>
#include <complex>

using namespace std::complex_literals;

std::tuple< arma::sp_cx_mat, arma::sp_cx_mat, arma::sp_cx_mat, arma::sp_cx_mat, arma::sp_cx_mat, arma::sp_cx_mat, arma::sp_cx_mat > make_operators(const double N, const double g, const double kappa, const double kappa_1, const double delta_c, const std::vector<double> delta_q, const double beta, const double gamma_dec, const double gamma_phi, const double LOPhi, std::vector<double> r_up, std::vector<double> r_down){
        // Initial bullshit, making necessary parameters from initial parameters
        arma::sp_cx_mat sigma_m(2, 2);
	sigma_m(1, 0) = 1.;//	= { {0. + 0i, 0. + 0i}, {1. + 0i, 0. + 0i} };
        arma::sp_cx_mat sigma_z(2, 2);
     	sigma_z(0, 0) = 1.;
	sigma_z(1, 1) = -1.;//	= { {1. + 0i, 0. + 0i}, {0. + 0i, -1. + 0i} };

        std::vector<double> gamma_p;
	for (unsigned int i = 0; i<N; ++i) gamma_p.push_back(2 * g * g * kappa / (kappa * kappa + pow(delta_c - delta_q.at(i), 2)));

        const double alpha = std::real(sqrt(2.0 * kappa_1) * beta / (kappa + 1i * delta_c));

        std::vector< arma::sp_cx_mat > a_ad(N);
        for (int i = 0; i < N; ++i)
               a_ad[i] = alpha * arma::speye(2, 2) - 1i * g * sigma_m / (kappa + 1i * (delta_c - delta_q[i]));
	
        std::vector < double > eps_q(N);
        for (int i = 0; i<N; ++i)
                eps_q[i] = (delta_c - delta_q[i]) * g * g * kappa / ( kappa * kappa + pow(delta_c - delta_q[i] , 2));
	

        std::vector < arma::sp_cx_mat > H_ap;
        for(int i = 0; i<N; ++i)
                H_ap.push_back(delta_q[i] * sigma_z / 2 + g * (alpha * sigma_m.t() + alpha * sigma_m) - eps_q[i] * sigma_m.t() * sigma_m);
	std::reverse(H_ap.begin(), H_ap.end()); //Building H from ground up, so lowest right diagonal is lowest state


        // Making the operators

        // c operators
        arma::sp_cx_mat c_1(2 * N, 2 * N);
        for (int i = 0; i<N; ++i) c_1.submat(2 * i, 2 * i, 2 * i + 1, 2 * i + 1) = sqrt(gamma_p[i]) * sigma_m;
	arma::sp_cx_mat complex_N_N_identity(N, N);
	complex_N_N_identity.speye();
        arma::sp_cx_mat c_2 = arma::kron(complex_N_N_identity, sqrt(gamma_dec) * sigma_m);

        arma::sp_cx_mat c_3 = arma::kron(complex_N_N_identity, sqrt(gamma_phi / 2) * sigma_z);

        arma::sp_cx_mat c(2 * N, 2 * N);
        for (int i = 0; i<N; ++i)
                c.submat(2 * i, 2 * i, 2 * i + 1, 2 * i + 1) = sqrt(2 * kappa_1) * (a_ad.at(i) - beta * arma::speye(2, 2)) * exp(-1i * LOPhi);
	c.clean(1e2 * arma::datum::eps);

        // H
        arma::sp_cx_mat H(2 * N, 2 * N);
        for (int i = 0; i<N; ++i)
                H.submat(2 * i, 2 * i, 2 * i + 1, 2 * i + 1) = H_ap.at(i);


        //Jump up and down operators
        arma::sp_cx_mat J_up(2 * N, 2 * N), J_down(2 * N, 2 * N);
        for (int i = 0; i < (N - 1); ++i){
                J_up.submat(2 * i, 2 * (i + 1), 2 * i + 1, 2 * (i + 1) + 1).eye();
                J_up.submat(2 * i, 2 * (i + 1), 2 * i + 1, 2 * (i + 1) + 1) *= sqrt(r_up.at(i)) ;
                J_down.submat(2 * (i + 1), 2 * i, 2 * (i + 1) + 1, 2 * i + 1).eye();
                J_down.submat(2 * (i + 1), 2 * i, 2 * (i + 1) + 1, 2 * i + 1) *= sqrt(r_down.at(i));
        }
	
        return {c_1, c_2, c_3, c, H, J_up, J_down};
}
