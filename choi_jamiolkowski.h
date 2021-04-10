#ifndef HAVE_CHOI_JAMIOLKOWSKI
#define HAVE_CHOI_JAMIOLKOWSKI
#include <complex>
#include <armadillo>
#include <cmath>



arma::mat spre(const arma::mat A);
arma::mat spost(const arma::mat A);
arma::mat sprepost(const arma::mat A, const arma::mat B);
arma::vec state_to_channel(const arma::mat rho);
arma::mat channel_to_state(const arma::vec rho_vec);


arma::cx_mat cx_spre(const arma::cx_mat A);
arma::cx_mat cx_spost(const arma::cx_mat A);
arma::cx_mat cx_sprepost(const arma::cx_mat A, const arma::cx_mat B);
arma::cx_vec cx_state_to_channel(const arma::cx_mat rho);
arma::cx_mat cx_channel_to_state(const arma::cx_vec rho_vec);


arma::sp_cx_mat sp_cx_spre(const arma::sp_cx_mat A);
arma::sp_cx_mat sp_cx_spost(const arma::sp_cx_mat A);
arma::sp_cx_mat sp_cx_sprepost(const arma::sp_cx_mat A, const arma::sp_cx_mat B);
arma::sp_cx_vec sp_cx_state_to_channel(const arma::sp_cx_mat rho);
arma::sp_cx_mat sp_cx_channel_to_state(const arma::sp_cx_vec rho_vec);

#endif
