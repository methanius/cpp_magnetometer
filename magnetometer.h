#ifndef HAVE_MAGNETOMETER
#define HAVE_MAGNETOMER
#include <iostream>
#include <vector>
#include <tuple>
#include <cmath>
#include <armadillo>
#include <complex>
#include "choi_jamiolkowski.h"


std::tuple< std::vector<double>, std::vector<double> > ehrenfest_chain(const int size, const double base_rate);

std::tuple< arma::sp_cx_mat, arma::sp_cx_mat, arma::sp_cx_mat, arma::sp_cx_mat, arma::sp_cx_mat, arma::sp_cx_mat, arma::sp_cx_mat > make_operators(const double N, const double g, const double kappa, const double kappa_1, const double delta_c, const std::vector<double> delta_q, const double beta, const double gamma_dec, const double gamma_phi, const double LOPhi, std::vector<double> r_up, std::vector<double> r_down);

std::tuple< std::vector<double>, std::vector<int>, std::vector<double> > homodyne_emission(const double N, const arma::sp_cx_mat H, const arma::sp_cx_mat J_up, const arma::sp_cx_mat J_down, const arma::sp_cx_mat c, const arma::sp_cx_mat c_1, const arma::sp_cx_mat c_2, const arma::sp_cx_mat c_3, const double eta, const int emission_startpoint, const std::vector<double> delta_q, const double dt, const double T);

std::tuple < arma::mat, arma::mat, arma::mat, arma::urowvec, arma::urowvec, std::vector<double> > homodyne_forward_backward_conditioning(const double T, const double dt, const double N, const arma::sp_cx_mat H, const arma::sp_cx_mat J_up, const arma::sp_cx_mat J_down, const arma::sp_cx_mat c, const arma::sp_cx_mat c_1, const arma::sp_cx_mat c_2, const arma::sp_cx_mat c_3, const double eta, const std::vector<double> dY, const std::vector<double> delta_q, const int frames);

#endif
