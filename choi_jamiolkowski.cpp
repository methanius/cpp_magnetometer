#include <complex>
#include <armadillo>
#include <cmath>

// doubles
arma::mat spre(const arma::mat A){
        return arma::kron(arma::eye(sqrt(A.size()), sqrt(A.size())), A.t());
}


arma::mat spost(const arma::mat A){
        return arma::kron(A, arma::eye(sqrt(A.size()), sqrt(A.size())));
}


arma::mat sprepost(const arma::mat A, const arma::mat B){
        return arma::kron(A, B.t());
}


arma::vec state_to_channel(const arma::mat rho){
	return rho.as_col();
}


arma::mat channel_to_state(arma::vec rho_vec){
	auto side_length = sqrt(rho_vec.n_elem);
	return arma::reshape(rho_vec, side_length, side_length);
}

// complex doubles
arma::cx_mat cx_spre(const arma::cx_mat A){
        return arma::kron(arma::eye(sqrt(A.size()), sqrt(A.size())), A.t());
}


arma::cx_mat cx_spost(const arma::cx_mat A){
        return arma::kron(A, arma::eye(sqrt(A.size()), sqrt(A.size())));
}


arma::cx_mat cx_sprepost(const arma::cx_mat A, const arma::cx_mat B){
        return arma::kron(A, B.t());
}


arma::cx_vec cx_state_to_channel(const arma::cx_mat rho){
	return rho.as_col();
}


arma::cx_mat cx_channel_to_state(arma::cx_vec rho_vec){
	auto side_length = sqrt(rho_vec.n_elem);
	return arma::reshape(rho_vec, side_length, side_length);
}

// sparse complex doubles
arma::sp_cx_mat sp_cx_spre(const arma::sp_cx_mat A){
	arma::sp_cx_mat I(sqrt(A.size()), sqrt(A.size()));
	I.speye();
        return arma::kron(I, A.t());
}


arma::sp_cx_mat sp_cx_spost(const arma::sp_cx_mat A){
	arma::sp_cx_mat I(sqrt(A.size()), sqrt(A.size()));
	I.speye();
        return arma::kron(A, I);
}


arma::sp_cx_mat sp_cx_sprepost(const arma::sp_cx_mat A, const arma::sp_cx_mat B){
        return arma::kron(A, B.t());
}


arma::sp_cx_vec sp_cx_state_to_channel(const arma::sp_cx_mat rho){
	return rho.as_col();
}


arma::sp_cx_mat sp_cx_channel_to_state(arma::sp_cx_vec rho_vec){
	auto side_length = sqrt(rho_vec.n_elem);
	return arma::reshape(rho_vec, side_length, side_length);
}

