#define ARMA_64BIT_WORD 1
#include <RcppDist.h>
// [[Rcpp::depends(RcppArmadillo, RcppDist, Rcpp, terra, base)]]
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;

// // [[Rcpp::export]]
// arma::rowvec mat2vec_arma(const arma::mat& m, bool byrow=false) {
//   NumericMatrix m_copy = clone(m);
//   
//   // stacks matrix columns on top of each other
//   if(byrow){
//     m_copy = transpose(m_copy);
//   }
//   NumericVector x(m_copy);
//   x.attr("dim") = R_NilValue;
//   return(x);
// }
// 
// [[Rcpp::export]]
NumericMatrix vector2matrix(NumericVector v, IntegerVector dim, bool byrow=false) {
  NumericVector v_copy = clone(v);

  // takes a vector and pulls out columns to a matrix
  NumericMatrix x(dim(0), dim(1));
  if(byrow){
    NumericMatrix x_t(dim(1), dim(0));
    std::copy(v_copy.begin(), v_copy.end(), x_t.begin());
    x = transpose(x_t);
  } else {
    std::copy(v_copy.begin(), v_copy.end(), x.begin());
  }
  return(x);
}



// // [[Rcpp::export]]
// arma::mat vector2matrix(const arma::rowvec& v, const arma::vec& dim, bool byrow=false) {
// 
//   // takes a vector and pulls out columns to a matrix
//   arma::mat x(dim(0), dim(1));
//   if(byrow){
//     NumericMatrix x_t(dim(1), dim(0));
//     std::copy(v_copy.begin(), v_copy.end(), x_t.begin());
//     x = transpose(x_t);
//   } else {
//     std::copy(v_copy.begin(), v_copy.end(), x.begin());
//   }
//   return(x);
// }




// [[Rcpp::export]]
arma::mat calc_R_theta_arma(const double& theta){
  
  // Calculate the rotation matrix R(theta)
  arma::mat R_theta = { {cos(theta), -sin(theta)}, {sin(theta), cos(theta)}};
  
  return(R_theta);
}


// [[Rcpp::export]]
arma::mat transform_s_arma(const arma::mat& s, const arma::mat& R_theta, const arma::rowvec& t, const arma::rowvec& mu_D){
  
  // Generate the rotated/translated s values
  arma::mat s_trans = s.each_row() - mu_D;
  s_trans = s_trans * R_theta.t();
  s_trans.each_row() += t + mu_D;
  
  // Return the transformed values as a matrix
  return s_trans;
}


// [[Rcpp::export]]
arma::mat transform_Y_arma(const arma::mat& Y, const arma::mat& R_theta, const arma::rowvec& t, const arma::rowvec& mu_D){
  
  // Generate the rotated/translated s values
  arma::mat Y_trans = Y.each_row() - t;
  Y_trans = Y_trans.each_row() - mu_D;
  Y_trans = Y_trans * R_theta;
  Y_trans.each_row() += mu_D;
  
  // Return the transformed values as a matrix
  return Y_trans;
}


// [[Rcpp::export]]
double update_sigma2_tune_arma(double sigma2_tune, const arma::vec& accept_reject,
                               const arma::uword& i, const double& min_val, const double& rate){
  
  double rate_mean = mean(accept_reject(arma::span((i-49), i))); // calculate the mean of the last 50 accept/reject values
  double delta_n = std::min(min_val, 1.0 / pow(i, rate));
  
  if(rate_mean < 0.44){
    sigma2_tune -= delta_n;
    if(sigma2_tune <= 0.0){
      sigma2_tune = delta_n;
    }
  } else {
    sigma2_tune += delta_n;
  }
  
  return (sigma2_tune);
}


// [[Rcpp::export]]
arma::rowvec Arma_colSums(const arma::mat& x) {
  return arma::sum(x, 0);
}

// [[Rcpp::export]]
arma::rowvec sample_t_arma(const arma::uvec& lambda, const double& sigma2, const arma::mat& Y,
                           const arma::mat& s, const arma::mat& R_theta, const double& sigma2_t,
                           const arma::rowvec& mu_D, const arma::vec m){
  
  arma::mat s_lambda = s.rows(lambda);
  
  arma::mat s_lambda_m = s_lambda.rows(arma::span(m(0), sum(m) - 1));
  s_lambda_m = s_lambda_m.each_row() - mu_D;
  s_lambda_m = s_lambda_m*R_theta.t();
  s_lambda_m = s_lambda_m.each_row() + mu_D;
  
  double sigma2_t_star = (sigma2_t/((arma::accu(m) - m(0))*sigma2_t + sigma2));
  arma::rowvec mu_t = Arma_colSums(Y.rows(arma::span(m(0), sum(m) - 1)) - s_lambda_m);
  mu_t = sigma2_t_star*mu_t;
  arma::mat Sigma_t = arma::eye(2,2);
  Sigma_t = Sigma_t*sigma2_t_star;
  
  arma::rowvec t = rmvnorm(1, as<arma::vec> (wrap(mu_t)), Sigma_t);
  return(t);
  
}



// [[Rcpp::export]]
double log_prob_arma(const double& theta, const double& theta_star, const arma::mat& S_mat, const double& kappa, const double& nu){
  
  double log_r = (kappa*cos(nu) + S_mat(0,0) + S_mat(1,1))*(cos(theta_star) - cos(theta)) +
    (kappa*sin(nu) - S_mat(0,1) + S_mat(1,0))*(sin(theta_star) - sin(theta));
  
  return(log_r);
  
}

// [[Rcpp::export]]
arma::mat sum_cube_slices(const arma::cube& x) {
  
  int n = x.n_slices;
  
  arma::mat result(2,2,arma::fill::zeros);
  
  for(int i = 0; i < n; i++) {
    result += x.slice(i);
  }
  
  return result;
}



// [[Rcpp::export]]
arma::mat S_mat_arma(const arma::mat& Y, const arma::rowvec& t, const double& sigma2, const arma::uvec& lambda,
                     const arma::mat& s, const arma::rowvec& mu_D, const arma::vec& m){
  
  arma::mat S_mat(2,2);
  
  arma::mat s_lambda = s.rows(lambda);
  arma::mat s_lambda_m = s_lambda.rows(arma::span(m(0), sum(m) - 1));
  s_lambda_m = s_lambda_m.each_row() - mu_D;
  arma::mat Y_m = Y.rows(arma::span(m(0), sum(m) - 1));
  Y_m = Y_m.each_row() - mu_D;
  Y_m = Y_m.each_row() - t;
  
  int n_m = Y_m.n_rows;
  
  arma::cube S_calc(2, 2, n_m);
  
  for(arma::uword i = 0; i < n_m; i++){
    S_calc.slice(i) = arma::trans(arma::trans(s_lambda_m.row(i))*Y_m.row(i));
  }
  
  S_mat = (0.5/sigma2)*sum_cube_slices(S_calc);
  
  return(S_mat);
  
}


// [[Rcpp::export]]
Rcpp::List sample_theta_arma(const arma::uvec& lambda, const double& sigma2, const arma::mat& Y, const double& theta,
                             const arma::mat& s, const arma::rowvec& t, const arma::rowvec& mu_D, const arma::vec m,
                             const double& kappa, const double& nu, const double& sigma2_theta){
  
  double theta_star = R::rnorm(theta, sqrt(sigma2_theta));
  double accept_reject = 0.0;
  
  if(abs(theta_star) < .01){
    
    arma::mat S_matrix = S_mat_arma(Y, t, sigma2, lambda, s, mu_D, m);
    double log_r = log_prob_arma(theta, theta_star, S_matrix, kappa, nu);
    
    if(log(R::runif(0,1)) < log_r){
      accept_reject = 1.0;
    }else {
      theta_star = theta;
    }
    
  }else {
    theta_star = theta;
  }
  
  return(Rcpp::List::create(Named("theta") = theta_star,
                            Named("accept_reject") = accept_reject));
  
}


// [[Rcpp::export]]
arma::vec rinvgammat_arma(const int& n, const arma::vec& range, const double& shape, const double& rate){
  
  double Fa = R::pgamma(1 / range(0), shape, rate, 0, 0);
  double Fb = R::pgamma(1 / range(1), shape, rate, 0, 0);
  
  arma::vec u(n);
  for(int i = 0; i < n; i++){
    u(i) = R::runif(Fa, Fb);
  }
  
  arma::vec invgamma(n);
  for(int i = 0; i < n; i++){
    invgamma(i) = 1 / R::qgamma(1 - u(i), shape, rate, 1, 0);
  }
  
  return(invgamma);
}



// // [[Rcpp::export]]
// double sample_sigma2_arma(const arma::mat& s, const arma::uvec& lambda, const double& c_sigma, const double& d_sigma,
//                           const arma::mat& Y, const arma::mat& R_theta, const arma::rowvec& t, const arma::rowvec& mu_D,
//                           const arma::vec& m, const arma::vec& range){
//   
//   // Rcout << "inside sample_sigma2 \n";
//   double n = Y.n_rows;
//   arma::mat s_lambda = s.rows(lambda);
//   // Rcout << "inside sample_sigma2 after s_lambda \n";
//   arma::mat s_lambda_m = s_lambda.rows(arma::span(m(0), sum(m) - 1));
//   // Rcout << "inside sample_sigma2 after s_lambda_m \n";
//   
//   s_lambda_m = transform_s_arma(s_lambda_m, R_theta, t, mu_D);
//   s_lambda = join_cols(s_lambda.rows(0, m(0) - 1), s_lambda_m);
//   // Rcout << "inside sample_sigma2 after s_lambda \n";
//   
//   double sigma2 = rinvgammat_arma(1, range, n + c_sigma, .5*arma::accu(pow(Y - s_lambda, 2)) + d_sigma)(0);
//   // double sigma2 = 1.0/R::rgamma(n + c_sigma, 1.0/(.5*arma::accu(pow(Y - s_lambda, 2)) + d_sigma));
//   
//   return(sigma2);
//   
// }

// [[Rcpp::export]]
double sample_sigma2_arma(const arma::mat& s, const arma::uvec& lambda, const double& c_sigma, const double& d_sigma,
                          const arma::mat& Y, const arma::mat& R_theta, const arma::rowvec& t, const arma::rowvec& mu_D,
                          const arma::vec& m, const arma::vec& range, const double& sigma2_prev){
  
  double n = Y.n_rows;
  arma::mat s_lambda = s.rows(lambda);
  arma::mat s_lambda_m = s_lambda.rows(arma::span(m(0), sum(m) - 1));
  
  s_lambda_m = transform_s_arma(s_lambda_m, R_theta, t, mu_D);
  s_lambda = join_cols(s_lambda.rows(0, m(0) - 1), s_lambda_m);
  
  // double sigma2 = rinvgammat_arma(1, range, n + c_sigma, .5*arma::accu(pow(Y - s_lambda, 2)) + d_sigma)(0);
  // double sigma2 = rinvgammat_arma2(range, n + c_sigma, .5*arma::accu(pow(Y - s_lambda, 2)) + d_sigma);
  double sigma2 = 1.0 / R::rgamma(n + c_sigma, 1.0 / (.5 * arma::accu(pow(Y - s_lambda, 2)) + d_sigma));
  if(sigma2 > range(1)){
    sigma2 = sigma2_prev;
  }
  
  return(sigma2);
  
}


// [[Rcpp::export]]
Rcpp::List create_clusters(const arma::uvec& lambda){
  
  // get unique cluster labels
  arma::uvec unique_labels = arma::unique(lambda);
  
  // create list of empty vectors for each cluster
  Rcpp::List clusters(lambda.n_elem);
  
  // loop through unique labels and add indices to corresponding cluster vector
  for (arma::uword i = 0; i < unique_labels.n_elem; ++i) {
    arma::uword label = unique_labels(i);
    arma::uvec indices = find(lambda == label);
    clusters(label) = indices;
  }
  
  return (clusters);
}

// [[Rcpp::export]]
Rcpp::List create_clusters2(const arma::uvec& lambda, const int& N) {
  
  // get unique cluster labels
  arma::uvec unique_labels = arma::unique(lambda);
  
  // check if the number of unique clusters exceeds N
  if (unique_labels.n_elem > N) {
    Rcpp::stop("Number of unique clusters exceeds the specified number of clusters N.");
  }
  
  // create list of empty vectors for each cluster
  Rcpp::List clusters(N);
  
  // loop through unique labels and add indices to corresponding cluster vector
  for (arma::uword i = 0; i < unique_labels.n_elem; ++i) {
    arma::uword label = unique_labels(i);
    arma::uvec indices = find(lambda == label);
    clusters(label) = indices;
  }
  
  return clusters;
}



// [[Rcpp::export]]
arma::uvec near_latents_arma(const arma::mat& s, const arma::rowvec& Y_k, const double& dist) {
  
  arma::vec x_diff = s.col(0) - Y_k(0);
  arma::vec y_diff = s.col(1) - Y_k(1);
  
  arma::uvec within_x = arma::find(arma::abs(x_diff) < dist);
  arma::uvec within_y = arma::find(arma::abs(y_diff) < dist);
  
  arma::uvec within_square = arma::intersect(within_x, within_y);
  
  return (within_square);
}


// // [[Rcpp::export]]
// Rcpp::List update_clusters_arma(const Rcpp::List& C_lambda, const arma::uword& lambda_k, const arma::uword& lambda_l, const arma::uword& k) {
//   
//   if(lambda_k == lambda_l){
//     return C_lambda;
//   } else{
//     Rcpp::List new_C_lambda = clone(C_lambda);
//     
//     // Remove record k from its current cluster
//     arma::uvec indices_k = Rcpp::as<arma::uvec>(new_C_lambda[lambda_k]);
//     indices_k = indices_k.elem(find(indices_k != k));
//     
//     // Check if cluster is empty and set to R_NilValue if it is
//     if(indices_k.is_empty()){
//       new_C_lambda[lambda_k] = R_NilValue;
//     } else {
//       new_C_lambda[lambda_k] = indices_k;
//     }
//     
//     // Add record k to the new cluster
//     if(new_C_lambda[lambda_l] == R_NilValue){
//       arma::uvec indices_l = {lambda_l};
//       new_C_lambda[lambda_l] = indices_l;
//     }else{
//       arma::uvec indices_l = Rcpp::as<arma::uvec>(new_C_lambda[lambda_l]);
//       indices_l = arma::join_cols(indices_l, arma::uvec({k}));
//       new_C_lambda[lambda_l] = indices_l;
//     }
//     
//     return new_C_lambda;
//   }
//   
// }


//[[Rcpp::export]]
arma::uvec sample_lambda_arma(const arma::uvec& lambda, const arma::mat& s, const double& sigma2, const arma::mat& R_theta,
                              const arma::rowvec& t, const arma::rowvec& mu_D, arma::mat Y, const arma::vec& file_id, const int& N,
                              const int& n, const arma::vec& m, const double& dist, const Rcpp::List& C_lambda){
  
  // Rcout << "inside sample_lambda \n";
  arma::mat s_trans = transform_s_arma(s, R_theta, t, mu_D);
  arma::uvec lambda_k = lambda;
  Rcpp::List new_C_lambda = Rcpp::clone(C_lambda);
  
  // Rcout << "inside sample_lambda after inits \n";
  
  // arma::wall_clock timer;
  // timer.tic();
  for(int k = 0; k < n; ++k){
    
    if(file_id(k) == 1){
      arma::uvec latent_candidates = near_latents_arma(s, Y.row(k), dist);
      // Rcout << "inside sample_lambda after latent_cands \n";
      if(latent_candidates.n_elem < 2){
        double updated_dist = dist;
        while(latent_candidates.n_elem < 2){
          updated_dist += 5;
          latent_candidates = near_latents_arma(s, Y.row(k), updated_dist);
        }
      }
      arma::mat s_k = s.rows(latent_candidates);
      arma::vec sum_sq_vec = sum(pow(s_k.each_row() - Y.row(k), 2), 1);
      sum_sq_vec = -0.5/sigma2 * sum_sq_vec;
      arma::vec gumbel = -log(-log(Rcpp::runif(latent_candidates.n_elem, 0, 1)));
      lambda_k(k) = latent_candidates(arma::index_max(sum_sq_vec + gumbel));
      // Rcout << "inside sample_lambda after lambda_k(k) \n";
      
    }else if(file_id(k) == 2){
      arma::uvec latent_candidates = near_latents_arma(s_trans, Y.row(k), dist);
      if(latent_candidates.n_elem < 2){
        double updated_dist = dist;
        while(latent_candidates.n_elem < 2){
          updated_dist += 5;
          latent_candidates = near_latents_arma(s_trans, Y.row(k), updated_dist);
        }
      }
      arma::mat s_trans_k = s_trans.rows(latent_candidates);
      arma::vec sum_sq_vec = sum(pow(s_trans_k.each_row() - Y.row(k), 2), 1);
      sum_sq_vec = -0.5/sigma2 * sum_sq_vec;
      arma::vec gumbel = -log(-log(Rcpp::runif(latent_candidates.n_elem, 0, 1)));
      lambda_k(k) = latent_candidates(arma::index_max(sum_sq_vec + gumbel));
      
    }
    
  }
  
  // double n_secs = timer.toc();
  // Rcout << "number of seconds: " << n_secs << " for full run lambda sampler \n";
  
  // Rcout << "inside sample_lambda before return(lambda_k) \n";
  return(lambda_k);
}


//[[Rcpp::export]]
List s_means_arma(const arma::uvec& lambda, const arma::mat new_Y, const int& N){
  
  arma::mat s_means(N, 2, arma::fill::zeros);
  arma::vec s_sizes(N, arma::fill::zeros);
  
  for(arma::uword i = 0; i < N; ++i){
    arma::uvec file_indices = find(lambda == i);
    if(file_indices.n_elem > 0){
      s_sizes(i) = as<double>(wrap(file_indices.n_elem));
      // Rcout << "double = " << as<double>(wrap(file_indices.n_elem)) << ", or = " << file_indices.n_elem << "\n";
      s_means.row(i) = Arma_colSums(new_Y.rows(file_indices)) / s_sizes(i);
    }
  }
  
  // Rcout << "s_sizes = " << s_sizes << "\n";
  // Rcout << "s_means = " << s_means << "\n";
  
  return Rcpp::List::create(Named("s_means") = s_means,
                            Named("s_sizes") = s_sizes);
}


//[[Rcpp::export]]
arma::rowvec propose_s_star_arma(const int& n_j, const arma::vec& s_j_mean, const double& sigma2, const arma::vec& D_bounds){
  
  arma::vec s_star(2);
  arma::mat cov(2, 2, arma::fill::eye);
  cov *= sigma2/n_j;
  s_star = arma::mvnrnd(s_j_mean, cov);
  
  while(s_star(0) < D_bounds(0) || s_star(0) > D_bounds(1) || s_star(1) < D_bounds(2) || s_star(1) > D_bounds(3)) {
    s_star = arma::mvnrnd(s_j_mean, cov);
  }
  
  return as<arma::rowvec>(wrap(s_star));
}


//[[Rcpp::export]]
arma::mat sample_s_arma(const arma::uvec& lambda, const double& sigma2, const arma::mat& Y, const arma::mat& Y_trans, 
                        const int& N, const arma::mat R_theta, const arma::rowvec& t,
                        const arma::rowvec& mu_D, const arma::vec& m, const int& n, const Rcpp::List& C_lambda,
                        const arma::vec& file_id, const arma::vec& D_bounds){
  
  // arma::wall_clock timer;
  // timer.tic();
  arma::mat Y_file1 = Y.rows(find(file_id == 1));
  // arma::mat Y_file2 = transform_Y_arma(Y.rows(find(file_id == 2)), R_theta, t, mu_D);
  arma::mat new_Y = arma::join_cols(Y_file1, Y_trans);
  // arma::mat new_Y = arma::join_cols(Y_file1, Y_file2);
  
  List s_ms = s_means_arma(lambda, new_Y, N);
  
  arma::vec s_sizes = as<arma::vec>(s_ms["s_sizes"]);
  arma::mat s_means = as<arma::mat>(s_ms["s_means"]);
  arma::mat s_star(N, 2);
  
  arma::uvec sizes_zero = find(s_sizes == 0);
  arma::uvec sizes_non_zero = find(s_sizes > 0);
  
  
  for(arma::uword i = 0; i < sizes_zero.n_elem; ++i){
    s_star.row(sizes_zero(i)) = {R::runif(D_bounds(0), D_bounds(1)), R::runif(D_bounds(2), D_bounds(3))};
  }
  
  for(arma::uword j = 0; j < sizes_non_zero.n_elem; ++j){
    s_star.row(sizes_non_zero(j)) = propose_s_star_arma(s_sizes(sizes_non_zero(j)),
                                                        as<arma::vec>(wrap(s_means.row(sizes_non_zero(j)))),
                                                        sigma2, D_bounds);
  }
  
  // double n_secs = timer.toc();
  // Rcout << "number of seconds: " << n_secs << " for full run s sampler \n";
  
  return s_star;
  
}


// [[Rcpp::export]]
void save_rds2(Rcpp::List obj, const std::string& filename) {
  
  Rcpp::Function saveRDS("saveRDS");
  saveRDS(obj, filename);
  
}


// [[Rcpp::export]]
void save_rds(Rcpp::List obj, const std::string& filename) {

  Rcpp::Environment base("package:base");
  Rcpp::Function saveRDS = base["saveRDS"];
  saveRDS(_["object"] = obj, _["file"] = filename);

}


// [[Rcpp::export]]
void write_csv_arma(const arma::rowvec& armaVector, const std::string& filePath) {

  Rcpp::Environment readr("package:readr");
  Rcpp::Function writeCsv = readr["write_csv"];

  // Convert Armadillo vector to a data frame with one column
  Rcpp::DataFrame dataFrame = Rcpp::DataFrame::create(Rcpp::Named("Column1") =  as<Rcpp::NumericVector>(wrap(armaVector)));
  writeCsv(_["x"] = dataFrame, _["file"] = filePath, _["append"] = true);

}


// [[Rcpp::export]]
Rcpp::DataFrame mat2df_arma(const arma::mat& X) {
  
  Rcpp::DataFrame X_df = as<Rcpp::DataFrame>(wrap(X));
  return(X_df);
  
} 

// [[Rcpp::export]]
arma::mat cube2mat_arma(const arma::cube& X_cube){
  
  arma::uword num_rows = X_cube.n_slices;
  arma::mat X_mat(num_rows, X_cube.n_cols * X_cube.n_rows);
  
  for(arma::uword i = 0; i < num_rows; i++){
    X_mat.row(i) = arma::conv_to<arma::rowvec>::from(arma::vectorise(X_cube.slice(i)));
    // X_mat.row(i) = arma::vectorise(X_cube.slice(i), 1);
  }
  
  return(X_mat);
}

// [[Rcpp::export]]
void write_csv_arma_df(const arma::mat& X, const std::string& filePath) {
  
  Rcpp::Environment readr("package:readr");
  Rcpp::Function writeCsv = readr["write_csv"];
  
  // Convert Armadillo vector to a data frame with one column
  Rcpp::DataFrame X_df = mat2df_arma(X);
  writeCsv(_["x"] = X_df, _["file"] = filePath, _["append"] = true);
  
}

// [[Rcpp::export]]
void write_csv_arma_vec(const arma::vec& X, const std::string& filePath) {
  
  Rcpp::Environment readr("package:readr");
  Rcpp::Function writeCsv = readr["write_csv"];
  
  // Convert Armadillo vector to a data frame with one column
  Rcpp::DataFrame X_df = as<Rcpp::DataFrame>(wrap(X));
  writeCsv(_["x"] = X_df, _["file"] = filePath, _["append"] = true);
  
}


// [[Rcpp::export]]
Rcpp::List run_mcmc_SpRL_linkage_arma(const int& n_iter, const arma::mat& Y_mat, const int& N, const arma::vec& m,
                                      const arma::vec& D_bounds, const double& dist, const Rcpp::List& init_vals,
                                      const Rcpp::List& hyperparameters, const arma::vec& sigma2_range,
                                      const std::string& file_name, bool verbose = true){
  
  // Calculate necessary quantities
  int n = Y_mat.n_rows; // Total number of records
  arma::rowvec mu_D = { (D_bounds(0) + D_bounds(1))/2.0, (D_bounds(2) + D_bounds(3))/2.0 }; // Midpoint of the spatial domain
  arma::mat Y = Y_mat.cols(0,1); // Observed spatial locations
  arma::vec file_id = Y_mat.col(3); // File IDs for each record
  
  
  // Create storage objects for the parameters
  arma::umat lambda_out(n_iter + 1, n); // Linkage structure
  arma::cube s_out(N, 2, n_iter + 1); // Latent locations
  arma::vec sigma2_out(n_iter + 1); // Location measurement error
  arma::vec theta_out(n_iter + 1); // Rotation for File 2
  arma::vec ar_theta(n_iter);
  arma::mat t_out(n_iter + 1, 2); // Translation for File 2
  
  
  // Hyperparameter values
  double c_sigma = hyperparameters["c_sigma"]; // sigma2
  double d_sigma = hyperparameters["d_sigma"]; 
  double kappa = hyperparameters["kappa"]; // theta
  double nu = hyperparameters["nu"]; 
  double sigma2_theta_min = as<double>(wrap(hyperparameters["sigma2_theta_min"])); 
  double sigma2_theta_rate = as<double>(wrap(hyperparameters["sigma2_theta_rate"])); 
  double sigma2_t = hyperparameters["sigma2_t"]; // t
  
  
  // Initialize the parameter values for the sampler
  arma::urowvec lambda_init = init_vals["lambda"];
  lambda_out.row(0) = lambda_init;
  arma::mat s_init = init_vals["s"];
  s_out.slice(0) = s_init;
  sigma2_out(0) = init_vals["sigma2"];
  theta_out(0) = init_vals["theta"];
  double sigma2_theta = as<double>(wrap(init_vals["sigma2_theta"]));
  arma::rowvec t_init = init_vals["t"];
  t_out.row(0) = t_init;
  
  
  // Create the linkage set C and the growth cluster index and initialize relevant storage
  Rcpp::List C_lambda = create_clusters(as<arma::uvec>(wrap(lambda_out.row(0))));
  Rcpp::List update_theta(2);
  arma::mat R_theta(2,2);
  
  
  if(verbose){
    Rcout << "Sampler Initialized. \n";
  }
  
  arma::vec iteration_timer(n_iter);
  arma::wall_clock timer;
  arma::wall_clock timer_full;
  timer_full.tic();
  
  for(arma::uword i = 0; i < n_iter; i++){
    
    timer.tic();
    // Calculate R(theta)
    R_theta = calc_R_theta_arma(theta_out(i));
    arma::mat Y_trans = transform_Y_arma(Y.rows(arma::span(m(0), sum(m) - 1)), R_theta, t_out.row(i), mu_D);
    

    // Sample Lambda
    lambda_out.row(i+1) = as<arma::urowvec>(wrap(sample_lambda_arma(as<arma::uvec>(wrap(lambda_out.row(i))), s_out.slice(i), sigma2_out(i),
                                                                    R_theta, t_out.row(i), mu_D, Y, file_id, N, n, m, dist, C_lambda)));
    C_lambda = create_clusters(as<arma::uvec>(wrap(lambda_out.row(i+1))));
  
    
    // Sample s
    s_out.slice(i+1) = sample_s_arma(as<arma::uvec>(wrap(lambda_out.row(i+1))), sigma2_out(i), Y, Y_trans, N, R_theta,
                                     t_out.row(i), mu_D, m, n, C_lambda, file_id, D_bounds);

    
    // Sample sigma2
    sigma2_out(i+1) = sample_sigma2_arma(s_out.slice(i+1), as<arma::uvec>(wrap(lambda_out.row(i+1))), c_sigma, d_sigma,
                                         Y, R_theta, t_out.row(i), mu_D, m, sigma2_range, sigma2_out(i));
    
    
    // Sample theta
    update_theta = sample_theta_arma(as<arma::uvec>(wrap(lambda_out.row(i+1))), sigma2_out(i+1), Y, theta_out(i), s_out.slice(i+1),
                                     t_out.row(i), mu_D, m, kappa, nu, sigma2_theta);
    theta_out(i+1) = as<double>(wrap(update_theta["theta"]));
    ar_theta(i) = update_theta["accept_reject"];
    if(ar_theta(i) == 1.0){
      R_theta = calc_R_theta_arma(theta_out(i+1));
    }
    
    
    // Sample t
    t_out.row(i+1) = sample_t_arma(as<arma::uvec>(wrap(lambda_out.row(i+1))), sigma2_out(i+1), Y,
                                   s_out.slice(i+1), R_theta, sigma2_t, mu_D, m);
  
    
    // Update MH proposal variances
    if(i > 0){
      if((i % 50) == 0){
        sigma2_theta = update_sigma2_tune_arma(sigma2_theta, ar_theta, i, sigma2_theta_min, sigma2_theta_rate);
        // save_rds(Rcpp::List::create(Named("lambda") = lambda_out,
        //                             Named("s") = s_out,
        //                             Named("theta") = theta_out,
        //                             Named("ar_theta") = ar_theta,
        //                             Named("t") = t_out,
        //                             Named("sigma2") = sigma2_out,
        //                             Named("MH_sigmas") = sigma2_theta),
        //          file_name);
      }
    }
    
  
    if(verbose & (i % 100 == 0)){
      Rcout << i+1 << " / " << n_iter << " iterations complete. \n";
    }
    iteration_timer(i) = timer.toc();
    if(sum(iteration_timer)/3600 > 143.5){
      // save_rds(Rcpp::List::create(Named("lambda") = lambda_out.row(i+1),
      //                             Named("s") = s_out.slice(i+1),
      //                             Named("theta") = theta_out(i+1),
      //                             Named("t") = t_out.row(i+1),
      //                             Named("sigma2") = sigma2_out(i+1),
      //                             Named("MH_sigmas") = sigma2_theta),
      //          file_name);
      save_rds(Rcpp::List::create(Named("lambda") = lambda_out,
                                  Named("s") = s_out,
                                  Named("theta") = theta_out,
                                  Named("ar_theta") = ar_theta,
                                  Named("t") = t_out,
                                  Named("sigma2") = sigma2_out,
                                  Named("MH_sigmas") = sigma2_theta),
               file_name);
    }
    // Rcout << "Number of seconds: " << n_secs << " for full sampler iteration \n";
  }
  
  double n_secs_full = timer_full.toc();
  Rcout << "Number of minutes: " << n_secs_full / 60.0 << " to complete full sampler \n";
  
  return(Rcpp::List::create(Named("lambda") = lambda_out,
                            Named("s") = s_out,
                            Named("theta") = theta_out,
                            Named("ar_theta") = ar_theta,
                            Named("t") = t_out,
                            Named("sigma2") = sigma2_out,
                            Named("MH_sigmas") = sigma2_theta));
  
}


//[[Rcpp::export]]
Rcpp::List run_mcmc_SpRL_linkage_arma_ft(const int& n_iter, const arma::mat& Y_mat, const int& N, const arma::vec& m,
                                         const arma::vec& D_bounds, const double& dist, const Rcpp::List& init_vals,
                                         const Rcpp::List& hyperparameters, const arma::vec& sigma2_range,
                                         const std::string& file_name_lambda, const std::string& file_name_s,
                                         const std::string& file_name_theta, const std::string& file_name_t,
                                         const std::string& file_name_sigma2, bool verbose = true){
  
  // Calculate necessary quantities
  int n = Y_mat.n_rows; // Total number of records
  arma::rowvec mu_D = { (D_bounds(0) + D_bounds(1))/2.0, (D_bounds(2) + D_bounds(3))/2.0 }; // Midpoint of the spatial domain
  arma::mat Y = Y_mat.cols(0,1); // Observed spatial locations
  arma::vec file_id = Y_mat.col(3); // File IDs for each record
  
  
  // Create storage objects for the parameters
  arma::umat lambda_out(n_iter + 1, n); // Linkage structure
  arma::cube s_out(N, 2, n_iter + 1); // Latent locations
  arma::vec sigma2_out(n_iter + 1); // Location measurement error
  arma::vec theta_out(n_iter + 1); // Rotation for File 2
  arma::mat t_out(n_iter + 1, 2); // Translation for File 2
  
  
  // Hyperparameter values
  double c_sigma = hyperparameters["c_sigma"]; // sigma2
  double d_sigma = hyperparameters["d_sigma"]; 
  double sigma2_t = hyperparameters["sigma2_t"]; // t
  
  
  // Initialize the parameter values for the sampler
  arma::urowvec lambda_init = init_vals["lambda"];
  lambda_out.row(0) = lambda_init;
  arma::mat s_init = init_vals["s"];
  s_out.slice(0) = s_init;
  sigma2_out(0) = init_vals["sigma2"];
  theta_out(0) = init_vals["theta"];
  double sigma2_theta = as<double>(wrap(init_vals["sigma2_theta"]));
  arma::rowvec t_init = init_vals["t"];
  t_out.row(0) = t_init;
  
  
  // Create the linkage set C and the growth cluster index and initialize relevant storage
  Rcpp::List C_lambda = create_clusters(as<arma::uvec>(wrap(lambda_out.row(0))));
  arma::mat R_theta = calc_R_theta_arma(theta_out(0));
  
  
  if(verbose){
    Rcout << "Sampler Initialized. \n";
  }
  
  arma::vec iteration_timer(n_iter);
  arma::wall_clock timer;
  arma::wall_clock timer_full;
  timer_full.tic();
  
  for(arma::uword i = 0; i < n_iter; i++){
    
    timer.tic();
    // Calculate transformed Y for file 2
    arma::mat Y_trans = transform_Y_arma(Y.rows(arma::span(m(0), sum(m) - 1)), R_theta, t_out.row(i), mu_D);
    
    
    // Sample Lambda
    lambda_out.row(i+1) = as<arma::urowvec>(wrap(sample_lambda_arma(as<arma::uvec>(wrap(lambda_out.row(i))), s_out.slice(i), sigma2_out(i),
                                                 R_theta, t_out.row(i), mu_D, Y, file_id, N, n, m, dist, C_lambda)));
    C_lambda = create_clusters(as<arma::uvec>(wrap(lambda_out.row(i+1))));
    
    
    // Sample s
    s_out.slice(i+1) = sample_s_arma(as<arma::uvec>(wrap(lambda_out.row(i+1))), sigma2_out(i), Y, Y_trans, N, R_theta,
                                     t_out.row(i), mu_D, m, n, C_lambda, file_id, D_bounds);
    
    
    // Sample sigma2
    sigma2_out(i+1) = sample_sigma2_arma(s_out.slice(i+1), as<arma::uvec>(wrap(lambda_out.row(i+1))), c_sigma, d_sigma,
                                         Y, R_theta, t_out.row(i), mu_D, m, sigma2_range, sigma2_out(i));
    
    
    // Sample t
    t_out.row(i+1) = sample_t_arma(as<arma::uvec>(wrap(lambda_out.row(i+1))), sigma2_out(i+1), Y,
                                   s_out.slice(i+1), R_theta, sigma2_t, mu_D, m);
    
    
    if(i == 0){
      Rcout << i+1 << " / " << n_iter << " iterations complete. \n";
      // write_csv_arma(as<arma::rowvec>(wrap(lambda_out.row(i))), file_name_lambda);
      // write_csv_arma(arma::conv_to<arma::rowvec>::from(arma::vectorise(s_out.slice(i))), file_name_s);
      // write_csv_arma(t_out.row(i), file_name_t);
      // save_rds(Rcpp::List::create(Named("lambda") = lambda_out,
      //                             Named("s") = s_out,
      //                             Named("t") = t_out,
      //                             Named("sigma2") = sigma2_out),
      //                             file_name);
      
    }
    
    if(verbose & ((i + 1) % 100 == 0)){
      Rcout << i+1 << " / " << n_iter << " iterations complete. \n";
      // write_csv_arma_df(as<arma::mat>(wrap(lambda_out.tail_rows(100))), file_name_lambda);
      // write_csv_arma_df(cube2mat_arma(s_out.tail_slices(100)), file_name_s);
      // write_csv_arma_df(t_out.tail_rows(100), file_name_t);
      if(i > 0){
        write_csv_arma_df(as<arma::mat>(wrap(lambda_out.rows(i - 99, i))), file_name_lambda);
        write_csv_arma_df(cube2mat_arma(s_out.slices(i - 99, i)), file_name_s);
        write_csv_arma_df(t_out.rows(i - 99, i), file_name_t);
        write_csv_arma_vec(sigma2_out.rows(i - 99, i), file_name_sigma2);
        // save_rds(Rcpp::List::create(Named("lambda") = lambda_out,
        //                             Named("s") = s_out,
        //                             Named("t") = t_out,
        //                             Named("sigma2") = sigma2_out),
        //                             file_name);
      }
    }
    
    
    iteration_timer(i) = timer.toc();
    if(sum(iteration_timer)/3600 > 167.5){
      
      Rcout << i+1 << " / " << n_iter << " iterations completed within time limit. \n";
      return(Rcpp::List::create(Named("lambda") = lambda_out,
                                Named("s") = s_out,
                                Named("t") = t_out,
                                Named("sigma2") = sigma2_out));
      // save_rds(Rcpp::List::create(Named("lambda") = lambda_out,
      //                             Named("s") = s_out,
      //                             Named("t") = t_out,
      //                             Named("sigma2") = sigma2_out),
      //                             file_name);
    }
    // Rcout << "Number of seconds: " << n_secs << " for full sampler iteration \n";
    
  }
  
  double n_secs_full = timer_full.toc();
  Rcout << "Number of minutes: " << n_secs_full / 60.0 << " to complete full sampler \n";
  
  return(Rcpp::List::create(Named("lambda") = lambda_out,
                            Named("s") = s_out,
                            Named("t") = t_out,
                            Named("sigma2") = sigma2_out));
  
}


//[[Rcpp::export]]
Rcpp::List run_mcmc_SpRL_linkage_arma_ft_timing(const int& n_iter, const arma::mat& Y_mat, const int& N, const arma::vec& m,
                                                const arma::vec& D_bounds, const double& dist, const Rcpp::List& init_vals,
                                                const Rcpp::List& hyperparameters, const arma::vec& sigma2_range,
                                                const std::string& file_name_lambda, const std::string& file_name_s,
                                                const std::string& file_name_theta, const std::string& file_name_t,
                                                const std::string& file_name_sigma2, bool verbose = true){
  
  // Calculate necessary quantities
  int n = Y_mat.n_rows; // Total number of records
  arma::rowvec mu_D = { (D_bounds(0) + D_bounds(1))/2.0, (D_bounds(2) + D_bounds(3))/2.0 }; // Midpoint of the spatial domain
  arma::mat Y = Y_mat.cols(0,1); // Observed spatial locations
  arma::vec file_id = Y_mat.col(3); // File IDs for each record
  
  
  // Create storage objects for the parameters
  arma::umat lambda_out(n_iter + 1, n); // Linkage structure
  arma::cube s_out(N, 2, n_iter + 1); // Latent locations
  arma::vec sigma2_out(n_iter + 1); // Location measurement error
  arma::vec theta_out(n_iter + 1); // Rotation for File 2
  arma::mat t_out(n_iter + 1, 2); // Translation for File 2
  
  
  // Hyperparameter values
  double c_sigma = hyperparameters["c_sigma"]; // sigma2
  double d_sigma = hyperparameters["d_sigma"]; 
  double sigma2_t = hyperparameters["sigma2_t"]; // t
  
  
  // Initialize the parameter values for the sampler
  arma::urowvec lambda_init = init_vals["lambda"];
  lambda_out.row(0) = lambda_init;
  arma::mat s_init = init_vals["s"];
  s_out.slice(0) = s_init;
  sigma2_out(0) = init_vals["sigma2"];
  theta_out(0) = init_vals["theta"];
  double sigma2_theta = as<double>(wrap(init_vals["sigma2_theta"]));
  arma::rowvec t_init = init_vals["t"];
  t_out.row(0) = t_init;
  
  
  // Create the linkage set C and the growth cluster index and initialize relevant storage
  Rcpp::List C_lambda = create_clusters(as<arma::uvec>(wrap(lambda_out.row(0))));
  arma::mat R_theta = calc_R_theta_arma(theta_out(0));
  
  
  if(verbose){
    Rcout << "Sampler Initialized. \n";
  }
  
  arma::vec iteration_timer(n_iter);
  arma::wall_clock timer;
  arma::wall_clock timer_full;
  timer_full.tic();
  
  for(arma::uword i = 0; i < n_iter; i++){
    
    timer.tic();
    // Calculate transformed Y for file 2
    arma::mat Y_trans = transform_Y_arma(Y.rows(arma::span(m(0), sum(m) - 1)), R_theta, t_out.row(i), mu_D);
    
    
    // Sample Lambda
    lambda_out.row(i+1) = as<arma::urowvec>(wrap(sample_lambda_arma(as<arma::uvec>(wrap(lambda_out.row(i))), s_out.slice(i), sigma2_out(i),
                                                 R_theta, t_out.row(i), mu_D, Y, file_id, N, n, m, dist, C_lambda)));
    C_lambda = create_clusters(as<arma::uvec>(wrap(lambda_out.row(i+1))));
    
    
    // Sample s
    s_out.slice(i+1) = sample_s_arma(as<arma::uvec>(wrap(lambda_out.row(i+1))), sigma2_out(i), Y, Y_trans, N, R_theta,
                t_out.row(i), mu_D, m, n, C_lambda, file_id, D_bounds);
    
    
    // Sample sigma2
    sigma2_out(i+1) = sample_sigma2_arma(s_out.slice(i+1), as<arma::uvec>(wrap(lambda_out.row(i+1))), c_sigma, d_sigma,
               Y, R_theta, t_out.row(i), mu_D, m, sigma2_range, sigma2_out(i));
    
    
    // Sample t
    t_out.row(i+1) = sample_t_arma(as<arma::uvec>(wrap(lambda_out.row(i+1))), sigma2_out(i+1), Y,
              s_out.slice(i+1), R_theta, sigma2_t, mu_D, m);
    
    
    if(i == 0){
      Rcout << i+1 << " / " << n_iter << " iterations complete. \n";
      // write_csv_arma(as<arma::rowvec>(wrap(lambda_out.row(i))), file_name_lambda);
      // write_csv_arma(arma::conv_to<arma::rowvec>::from(arma::vectorise(s_out.slice(i))), file_name_s);
      // write_csv_arma(t_out.row(i), file_name_t);
      // save_rds(Rcpp::List::create(Named("lambda") = lambda_out,
      //                             Named("s") = s_out,
      //                             Named("t") = t_out,
      //                             Named("sigma2") = sigma2_out),
      //                             file_name);
      
    }
    
    if(verbose & ((i + 1) % 100 == 0)){
      Rcout << i+1 << " / " << n_iter << " iterations complete. \n";
      // write_csv_arma_df(as<arma::mat>(wrap(lambda_out.tail_rows(100))), file_name_lambda);
      // write_csv_arma_df(cube2mat_arma(s_out.tail_slices(100)), file_name_s);
      // write_csv_arma_df(t_out.tail_rows(100), file_name_t);
      if(i > 0){
        // write_csv_arma_df(as<arma::mat>(wrap(lambda_out.rows(i - 99, i))), file_name_lambda);
        // write_csv_arma_df(cube2mat_arma(s_out.slices(i - 99, i)), file_name_s);
        // write_csv_arma_df(t_out.rows(i - 99, i), file_name_t);
        // write_csv_arma_vec(sigma2_out.rows(i - 99, i), file_name_sigma2);
        // save_rds(Rcpp::List::create(Named("lambda") = lambda_out,
        //                             Named("s") = s_out,
        //                             Named("t") = t_out,
        //                             Named("sigma2") = sigma2_out),
        //                             file_name);
      }
    }
    
    
    iteration_timer(i) = timer.toc();
    if(sum(iteration_timer)/3600 > 167.5){
      
      Rcout << i+1 << " / " << n_iter << " iterations completed within time limit. \n";
      return(Rcpp::List::create(Named("lambda") = lambda_out,
                                Named("s") = s_out,
                                Named("t") = t_out,
                                Named("sigma2") = sigma2_out));
      // save_rds(Rcpp::List::create(Named("lambda") = lambda_out,
      //                             Named("s") = s_out,
      //                             Named("t") = t_out,
      //                             Named("sigma2") = sigma2_out),
      //                             file_name);
    }
    // Rcout << "Number of seconds: " << n_secs << " for full sampler iteration \n";
    
  }
  
  double n_secs_full = timer_full.toc();
  Rcout << "Number of minutes: " << n_secs_full / 60.0 << " to complete full sampler \n";
  
  write_csv_arma_df(as<arma::mat>(wrap(lambda_out)), file_name_lambda);
  write_csv_arma_df(cube2mat_arma(s_out), file_name_s);
  write_csv_arma_df(t_out, file_name_t);
  write_csv_arma_vec(sigma2_out, file_name_sigma2);
  
  return(Rcpp::List::create(Named("lambda") = lambda_out,
                            Named("s") = s_out,
                            Named("t") = t_out,
                            Named("sigma2") = sigma2_out,
                            Named("iter_timing") = iteration_timer));
  
}



//[[Rcpp::export]]
Rcpp::List run_mcmc_SpRL_linkage_arma_ft_ns(const int& n_iter, const arma::mat& Y_mat, const int& N, const arma::vec& m,
                                            const arma::vec& D_bounds, const double& dist, const Rcpp::List& init_vals,
                                            const Rcpp::List& hyperparameters, const arma::vec& sigma2_range,
                                            const std::string& file_name_lambda, const std::string& file_name_s,
                                            const std::string& file_name_theta, const std::string& file_name_t,
                                            const std::string& file_name_sigma2, bool verbose = true){
  
  // Calculate necessary quantities
  int n = Y_mat.n_rows; // Total number of records
  arma::rowvec mu_D = { (D_bounds(0) + D_bounds(1))/2.0, (D_bounds(2) + D_bounds(3))/2.0 }; // Midpoint of the spatial domain
  arma::mat Y = Y_mat.cols(0,1); // Observed spatial locations
  arma::vec file_id = Y_mat.col(3); // File IDs for each record
  
  
  // Create storage objects for the parameters
  arma::umat lambda_out(n_iter + 1, n); // Linkage structure
  arma::cube s_out(N, 2, n_iter + 1); // Latent locations
  arma::vec sigma2_out(n_iter + 1); // Location measurement error
  arma::vec theta_out(n_iter + 1); // Rotation for File 2
  arma::mat t_out(n_iter + 1, 2); // Translation for File 2
  
  
  // Hyperparameter values
  double c_sigma = hyperparameters["c_sigma"]; // sigma2
  double d_sigma = hyperparameters["d_sigma"]; 
  double sigma2_t = hyperparameters["sigma2_t"]; // t
  
  
  // Initialize the parameter values for the sampler
  arma::urowvec lambda_init = init_vals["lambda"];
  lambda_out.row(0) = lambda_init;
  arma::mat s_init = init_vals["s"];
  s_out.slice(0) = s_init;
  sigma2_out(0) = init_vals["sigma2"];
  theta_out(0) = init_vals["theta"];
  double sigma2_theta = as<double>(wrap(init_vals["sigma2_theta"]));
  arma::rowvec t_init = init_vals["t"];
  t_out.row(0) = t_init;
  
  
  // Create the linkage set C and the growth cluster index and initialize relevant storage
  Rcpp::List C_lambda = create_clusters(as<arma::uvec>(wrap(lambda_out.row(0))));
  arma::mat R_theta = calc_R_theta_arma(theta_out(0));
  
  
  if(verbose){
    Rcout << "Sampler Initialized. \n";
  }
  
  arma::vec iteration_timer(n_iter);
  arma::wall_clock timer;
  arma::wall_clock timer_full;
  timer_full.tic();
  
  for(arma::uword i = 0; i < n_iter; i++){
    
    timer.tic();
    // Calculate transformed Y for file 2
    arma::mat Y_trans = transform_Y_arma(Y.rows(arma::span(m(0), sum(m) - 1)), R_theta, t_out.row(i), mu_D);
    
    
    // Sample Lambda
    lambda_out.row(i+1) = as<arma::urowvec>(wrap(sample_lambda_arma(as<arma::uvec>(wrap(lambda_out.row(i))), s_out.slice(i), sigma2_out(i),
                                                 R_theta, t_out.row(i), mu_D, Y, file_id, N, n, m, dist, C_lambda)));
    C_lambda = create_clusters(as<arma::uvec>(wrap(lambda_out.row(i+1))));
    
    
    // Sample s
    s_out.slice(i+1) = sample_s_arma(as<arma::uvec>(wrap(lambda_out.row(i+1))), sigma2_out(i), Y, Y_trans, N, R_theta,
                t_out.row(i), mu_D, m, n, C_lambda, file_id, D_bounds);
    
    
    // Sample sigma2
    sigma2_out(i+1) = sample_sigma2_arma(s_out.slice(i+1), as<arma::uvec>(wrap(lambda_out.row(i+1))), c_sigma, d_sigma,
               Y, R_theta, t_out.row(i), mu_D, m, sigma2_range, sigma2_out(i));
    
    
    // Sample t
    t_out.row(i+1) = sample_t_arma(as<arma::uvec>(wrap(lambda_out.row(i+1))), sigma2_out(i+1), Y,
              s_out.slice(i+1), R_theta, sigma2_t, mu_D, m);
    
    
    if(i == 0){
      Rcout << i+1 << " / " << n_iter << " iterations complete. \n";
      // write_csv_arma(as<arma::rowvec>(wrap(lambda_out.row(i))), file_name_lambda);
      // write_csv_arma(arma::conv_to<arma::rowvec>::from(arma::vectorise(s_out.slice(i))), file_name_s);
      // write_csv_arma(t_out.row(i), file_name_t);
      // save_rds(Rcpp::List::create(Named("lambda") = lambda_out,
      //                             Named("s") = s_out,
      //                             Named("t") = t_out,
      //                             Named("sigma2") = sigma2_out),
      //                             file_name);
      
    }
    
    if(verbose & ((i + 1) % 100 == 0)){
      Rcout << i+1 << " / " << n_iter << " iterations complete. \n";
      // write_csv_arma_df(as<arma::mat>(wrap(lambda_out.tail_rows(100))), file_name_lambda);
      // write_csv_arma_df(cube2mat_arma(s_out.tail_slices(100)), file_name_s);
      // write_csv_arma_df(t_out.tail_rows(100), file_name_t);
      if(i > 0){
        // write_csv_arma_df(as<arma::mat>(wrap(lambda_out.rows(i - 99, i))), file_name_lambda);
        // write_csv_arma_df(cube2mat_arma(s_out.slices(i - 99, i)), file_name_s);
        // write_csv_arma_df(t_out.rows(i - 99, i), file_name_t);
        // write_csv_arma_vec(sigma2_out.rows(i - 99, i), file_name_sigma2);
        // save_rds(Rcpp::List::create(Named("lambda") = lambda_out,
        //                             Named("s") = s_out,
        //                             Named("t") = t_out,
        //                             Named("sigma2") = sigma2_out),
        //                             file_name);
      }
    }
    
    
    iteration_timer(i) = timer.toc();
    if(sum(iteration_timer)/3600 > 167.5){
      
      Rcout << i+1 << " / " << n_iter << " iterations completed within time limit. \n";
      return(Rcpp::List::create(Named("lambda") = lambda_out,
                                Named("s") = s_out,
                                Named("t") = t_out,
                                Named("sigma2") = sigma2_out));
      // save_rds(Rcpp::List::create(Named("lambda") = lambda_out,
      //                             Named("s") = s_out,
      //                             Named("t") = t_out,
      //                             Named("sigma2") = sigma2_out),
      //                             file_name);
    }
    // Rcout << "Number of seconds: " << n_secs << " for full sampler iteration \n";
    
  }
  
  double n_secs_full = timer_full.toc();
  Rcout << "Number of minutes: " << n_secs_full / 60.0 << " to complete full sampler \n";
  
  return(Rcpp::List::create(Named("lambda") = lambda_out,
                            Named("s") = s_out,
                            Named("t") = t_out,
                            Named("sigma2") = sigma2_out));
  
}

// Function to calculate Euclidean distance between two points
// [[Rcpp::export]]
double euclidean_dist_arma(double x1, double y1, double x2, double y2) {
  return sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2));
}

// [[Rcpp::export]]
arma::uword closest_point_arma(const arma::mat& file_1, const arma::rowvec& file_2_k, const double& dist) {
  
  double check_dist = dist;
  
  arma::vec x_diff = file_1.col(0) - file_2_k(0);
  arma::vec y_diff = file_1.col(1) - file_2_k(1);
  
  arma::uvec within_x = arma::find(arma::abs(x_diff) <= dist);
  arma::uvec within_y = arma::find(arma::abs(y_diff) <= dist);
  
  arma::uvec within_square = arma::intersect(within_x, within_y);
  
  while(within_square.is_empty()){

    check_dist *= 2.0;
    within_x = arma::find(arma::abs(x_diff) <= check_dist);
    within_y = arma::find(arma::abs(y_diff) <= check_dist);

    within_square = arma::intersect(within_x, within_y);

  }
  
  arma::vec near_dists(within_square.n_elem);
  
  for(arma::uword i = 0; i < within_square.n_elem; i++){
    
    near_dists(i) = euclidean_dist_arma(file_2_k(0), file_2_k(1), file_1(within_square(i), 0), file_1(within_square(i), 1));
    
  }
  
  arma::uword nearest_point = within_square(arma::index_min(near_dists));
  
  return (nearest_point);
}


// [[Rcpp::export]]
arma::uvec nearest_distance_matching(const arma::mat& file_1, const arma::mat& file_2,
                                    const double& dist){
  
  arma::uvec lambda(file_1.n_rows + file_2.n_rows);
  
  for(arma::uword i = 0; i < file_1.n_rows; i ++){
    lambda(i) = i;
  }
  
  for(arma::uword j = 0; j < file_2.n_rows; j++){
    lambda(file_1.n_rows + j) = closest_point_arma(file_1, file_2.row(j), dist);
  }
  
  return(lambda);
  
}  

// // [[Rcpp::export]]
// arma::uvec growth_cluster_check(const Rcpp::List& clusters, const arma::mat& V, double q) {
//   
//   arma::uword n_clusters = clusters.size();
//   arma::uvec growth_clusters(n_clusters);
//   arma::vec ids = V.col(1);
//   arma::vec vols = V.col(0);
//   
//   for (arma::uword i = 0; i < n_clusters; ++i) {
//     
//     // Check if cluster is empty
//     if (clusters[i] == R_NilValue) {
//       growth_clusters(i) = 0;
//     } else {
//       arma::uvec indices = clusters[i];
//       arma::vec file_ids = ids(indices);
//       
//       arma::vec unique_files = arma::unique(file_ids);
//       if (unique_files.n_elem > 1) {
//         
//         arma::vec cluster_volumes = vols(indices);
//         arma::vec file_volumes(unique_files.n_elem, arma::fill::zeros);
//         
//         for (arma::uword j = 0; j < unique_files.n_elem; ++j) {
//           arma::uvec file_indices = find(file_ids == unique_files(j));
//           file_volumes(j) = arma::sum(cluster_volumes(file_indices));
//         }
//         
//         arma::uword min_file_index = arma::index_min(unique_files);
//         arma::uword max_file_index = arma::index_max(unique_files);
//         
//         double min_file_volume = file_volumes(min_file_index);
//         double max_file_volume = file_volumes(max_file_index);
//         
//         if (max_file_volume >= q * min_file_volume) {
//           growth_clusters(i) = 1;
//         } else {
//           growth_clusters(i) = 0;
//         }
//         
//       } else {
//         growth_clusters(i) = 0;
//       }
//     }
//   }
//   
//   return growth_clusters;
// }
// 
// 
// // Define the Q(phi) function
// // [[Rcpp::export]]
// double Q_func(const double& d, const double& phi) {
//   if (d <= phi) {
//     return 1.0 - (1.5) * (d / phi) + (.5) * pow((d / phi), 3.0);
//   } else {
//     return 0.0;
//   }
// }
// 
// 
// // Compute spatial covariance matrix for growth clusters
// // [[Rcpp::export]]
// arma::mat spatial_covariance(const arma::mat& s, const arma::uvec& gc_tracker, const double& phi) {
//   
//   // Get number of growth clusters
//   arma::uword n_gc = arma::sum(gc_tracker);
//   
//   // Initialize spatial covariance matrix
//   arma::mat Q(n_gc, n_gc, arma::fill::zeros);
//   
//   // Loop over all pairs of growth clusters
//   arma::uword i = 0;
//   arma::uword j = 0;
//   for (arma::uword c = 0; c < gc_tracker.n_elem; ++c) {
//     if (gc_tracker(c) == 1) {
//       j = 0;
//       for (arma::uword d = 0; d < gc_tracker.n_elem; ++d) {
//         if (gc_tracker(d) == 1) {
//           double distance = arma::norm(s.row(c) - s.row(d), 2);
//           Q(i, j) = Q_func(distance, phi);
//           ++j;
//         } 
//       }
//       ++i;
//     }  
//   }
//   
//   return (Q);
// } 
// 
// 
// // [[Rcpp::export]]
// arma::mat calculate_A(const Rcpp::List& C_lambda, const arma::mat& V, const arma::uvec& gc_tracker, const double& b) {
//   
//   arma::uword n_gc = arma::sum(gc_tracker);
//   arma::mat A(n_gc, n_gc, arma::fill::zeros);
//   arma::vec vols = V.col(0);
//   arma::vec ids = V.col(1);
//   
//   arma::uword c_index = 0;
//   
//   for (arma::uword c1 = 0; c1 < C_lambda.size(); ++c1) {
//     
//     if (gc_tracker(c1) == 1) {
//       
//       arma::uvec indices1 = C_lambda[c1];
//       arma::vec file_ids1 = ids(indices1);
//       arma::vec vols_ids1 = vols(indices1);
//       arma::vec unique_files1 = arma::unique(file_ids1);
//       
//       arma::uword min_file_index1 = arma::index_min(unique_files1);
//       double v_min1 = arma::sum(vols_ids1(find(file_ids1 == unique_files1(min_file_index1))));
//       
//       double a_c1 = v_min1 / (b + v_min1);
//       
//       A(c_index, c_index) = a_c1;
//       
//       c_index++;
//     }
//   }
//   
//   return (A);
// } 
// 
// 
// // [[Rcpp::export]]
// arma::mat calculate_W(const Rcpp::List& C_lambda, const arma::mat& V, const arma::uvec& gc_tracker, const double& v_max) {
//   
//   arma::uword n_gc = arma::sum(gc_tracker);
//   arma::mat W(n_gc, n_gc, arma::fill::zeros);
//   arma::vec vols = V.col(0);
//   arma::vec ids = V.col(1);
//   
//   arma::uword c_index = 0;
//   
//   for (arma::uword c1 = 0; c1 < C_lambda.size(); ++c1) {
//     
//     if (gc_tracker(c1) == 1) {
//       
//       arma::uvec indices1 = C_lambda[c1];
//       arma::vec file_ids1 = ids(indices1);
//       arma::vec vols_ids1 = vols(indices1);
//       arma::vec unique_files1 = arma::unique(file_ids1);
//       
//       arma::uword min_file_index1 = arma::index_min(unique_files1);
//       double v_min1 = arma::sum(vols_ids1(find(file_ids1 == unique_files1(min_file_index1))));
//       
//       // double w_c1 = v_min1 / v_max;
//       double w_c1 = v_min1;
//       
//       // W(c_index, c_index) = w_c1;
//       W(c_index, c_index) = 1.0;
//       
//       c_index++;
//     }
//   }
//   
//   return (W);
// } 
// 
// 
// 
// // [[Rcpp::export]]
// arma::vec calculate_growth(const Rcpp::List& C_lambda, const arma::uvec& gc_tracker, const arma::mat& V, const double& int_length) {
//   
//   arma::uword n_clusters = C_lambda.size();
//   arma::vec growth(n_clusters, arma::fill::zeros);
//   arma::vec ids = V.col(1);
//   arma::vec vols = V.col(0);
//   
//   for (arma::uword i = 0; i < n_clusters; ++i) {
//     
//     if (gc_tracker(i) == 1) {
//       arma::uvec indices = C_lambda[i];
//       arma::vec file_ids = ids(indices);
//       
//       arma::vec unique_files = arma::unique(file_ids);
//       arma::vec cluster_volumes = vols(indices);
//       arma::vec file_volumes(unique_files.n_elem, arma::fill::zeros);
//       
//       for (arma::uword j = 0; j < unique_files.n_elem; ++j) {
//         arma::uvec file_indices = find(file_ids == unique_files(j));
//         file_volumes(j) = arma::sum(cluster_volumes(file_indices));
//       }
//       
//       arma::uword min_file_index = arma::index_min(unique_files);
//       arma::uword max_file_index = arma::index_max(unique_files);
//       
//       double min_file_volume = file_volumes(min_file_index);
//       double max_file_volume = file_volumes(max_file_index);
//       
//       growth(i) = (max_file_volume - min_file_volume) / (int_length * (max_file_index - min_file_index));
//       
//     }
//   } 
//   
//   return (growth(find(gc_tracker == 1)));
// }
// 
// 
//[[Rcpp::export]]
arma::mat update_covars_arma(const arma::mat& s, const Rcpp::List& rasters){

  arma::mat X_s(s.n_rows, rasters.size() + 1);
  X_s.col(0).fill(1);

  Rcpp::Environment terra("package:terra");
  Rcpp::Function extract = terra["extract"];

  for(int i = 1; i < rasters.size() + 1; ++i){
    Rcpp::List vals = extract(_["x"] = rasters[i - 1], _["y"] = s, _["method"] = "bilinear");
    X_s.col(i) = as<arma::colvec>(wrap(vals[0]));
  }

  return(X_s);
}


//[[Rcpp::export]]
arma::mat pred_means(const arma::mat& X, const arma::vec& beta, const double& gamma,
                     const double& alpha, const arma::vec& V, const arma::uword N){
  
  arma::mat A(N, N, arma::fill::zeros);
  
  for (arma::uword i = 0; i < N; ++i) {

    double a_ii = pow(V(i), alpha) / (pow(gamma, alpha) + pow(V(i), alpha));
    A(i, i) = a_ii;

  }
  
  arma::mat pred_means = A * X * beta;
  
  return(pred_means);
  
}


//[[Rcpp::export]]
arma::vec pred_means_add(const arma::mat& X, const arma::vec& beta, const double& gamma,
                         const double& alpha, const arma::vec& V, const arma::uword N){
  
  arma::vec A(N, arma::fill::zeros);
  
  for (arma::uword i = 0; i < N; ++i) {
    
    double a_i = (beta(0) * pow(V(i), alpha)) / (pow(gamma, alpha) + pow(V(i), alpha));
    A(i) = a_i;
    
  }
  
  arma::vec pred_means = A + X * beta;
  
  return(pred_means);
  
}


//[[Rcpp::export]]
arma::mat generate_replicates(const arma::mat& X, const arma::mat& beta, const arma::vec& tau, const arma::vec& gamma,
                              const arma::vec& alpha, const arma::vec& V, const arma::uword N){
  
  arma::mat replicates(beta.n_rows, N, arma::fill::zeros);
  arma::mat covariance_mat(N, N, arma::fill::eye);
  arma::vec rep_vec(N, arma::fill::zeros);
  arma::vec norm_vec(N, arma::fill::zeros);
  
  for(arma::uword i = 0; i < beta.n_rows; ++i){
    rep_vec = pred_means(X, as<arma::vec>(wrap(beta.row(i))), gamma(i),
                            alpha(i), V, N);
    norm_vec = arma::mvnrnd(rep_vec, tau(i)*covariance_mat);
    // Rcout << "after norm_vec \n";
    replicates.row(i) = as<arma::rowvec>(wrap(norm_vec));
    // Rcout << "after norm_vec \n";
  }
  
  return replicates;
  
}


//[[Rcpp::export]]
Rcpp::List calculate_crps(const arma::mat& y_rep, const arma::mat& y_rep_star, const arma::vec& y){
  
  arma::vec crps_hat(y.n_elem, arma::fill::zeros);
  
  for(arma::uword i = 0; i < y.n_elem; ++i){
    crps_hat(i) = arma::mean(arma::abs(y(i) - y_rep.col(i))) - 0.5 * arma::mean(arma::abs(y_rep.col(i) - y_rep_star.col(i)));
  }
  
  double crps_est = arma::mean(crps_hat);
  double crps_se = arma::stddev(crps_hat) / sqrt(as<double>(wrap(y.n_elem)));
  
  return(Rcpp::List::create(Named("CRPS_est") = crps_est,
                            Named("CRPS_est_se") = crps_se,
                            Named("pointwise_estimates") = crps_hat));
           
} 







// 
// 
// //[[Rcpp::export]]
// arma::rowvec sample_beta_arma(const arma::mat& X_s, const arma::mat& A, const arma::mat& W, 
//                               const arma::mat Q_phi, const arma::colvec& G, const arma::mat& Sigma_beta_inv, const double& gamma2,
//                               const double& tau2){
//   
//   arma::vec mu = inv_sympd(X_s.t()*A*inv_sympd(gamma2*Q_phi + tau2*W)*A*X_s + Sigma_beta_inv)*X_s.t()*A*inv_sympd(gamma2*Q_phi + tau2*W)*G;
//   arma::mat cov = inv_sympd(X_s.t()*A*inv_sympd(gamma2*Q_phi + tau2*W)*A*X_s + Sigma_beta_inv);
//   
//   arma::vec beta = arma::mvnrnd(mu, cov);
//   
//   return(as<arma::rowvec>(wrap(beta)));
// }
// 
// 
// //[[Rcpp::export]]
// Rcpp::List sample_b_arma(const Rcpp::List& C_lambda, const arma::mat& V, const arma::uvec& gc_tracker,
//                          const arma::mat& X_s, const arma::mat& A, const arma::mat& Q_phi, const arma::mat& W,
//                          const arma::colvec& G, const arma::vec& beta, const double& gamma2, const double& tau2,
//                          const double& u, const double& w, const double& b, const double& sigma2_b){
//   
//   double b_star = R::rnorm(b, sqrt(sigma2_b));
//   double accept_reject = 0.0;
//   if(b_star < 0){
//     return(Rcpp::List::create(Named("b") = b,
//                               Named("accept_reject") = accept_reject));
//   }else {
//     arma::mat new_A = calculate_A(C_lambda, V, gc_tracker, b_star); 
//     double log_r = (u - 1)*(log(b_star) - log(b)) - w*(b_star - b) -
//       0.5*as_scalar(trans(G - new_A*X_s*beta)*inv_sympd(gamma2*Q_phi + tau2*W)*(G - new_A*X_s*beta) - 
//       trans(G - A*X_s*beta)*inv_sympd(gamma2*Q_phi + tau2*W)*(G - A*X_s*beta));
//     if(log(R::runif(0,1)) < log_r){
//       accept_reject = 1.0;
//     }else {
//       b_star = b;
//     }
//     
//   }
//   
//   return(Rcpp::List::create(Named("b") = b_star,
//                             Named("accept_reject") = accept_reject));
// } 
// 
// 
// //[[Rcpp::export]]
// Rcpp::List sample_phi_arma(const Rcpp::List& C_lambda, const arma::mat& V, const arma::uvec& gc_tracker,
//                            const arma::mat& X_s, const arma::mat& A, const arma::mat& Q_phi, const arma::mat& W,
//                            const arma::colvec& G, const arma::vec& beta, const double& gamma2, const double& tau2,
//                            const arma::mat& s, const double& r, const double& phi, const double& sigma2_phi){
//   
//   double phi_star = R::rnorm(phi, sqrt(sigma2_phi));
//   double accept_reject = 0.0;
//   if(phi_star < 0 || phi_star > r){
//     return(Rcpp::List::create(Named("phi") = phi,
//                               Named("accept_reject") = accept_reject));
//   }else {
//     arma::mat new_Q_phi = spatial_covariance(s, gc_tracker, phi_star); 
//     double log_r = -0.5 * as_scalar(trans(G - A*X_s*beta)*inv_sympd(gamma2*new_Q_phi + tau2*W)*(G - A*X_s*beta) - 
//                                     trans(G - A*X_s*beta)*inv_sympd(gamma2*Q_phi + tau2*W)*(G - A*X_s*beta)) -
//                     0.5 * (log(det(gamma2 * new_Q_phi + tau2 * W)) - log(det(gamma2 * Q_phi + tau2 * W)));
//     // Rcout << "log_r val = " << log_r << "\n";
//     if(log(R::runif(0,1)) < log_r){
//       accept_reject = 1.0;
//     }else {
//       phi_star = phi;
//     }
//     
//   }
//   
//   return(Rcpp::List::create(Named("phi") = phi_star,
//                             Named("accept_reject") = accept_reject));
// }
// 
// 
// //[[Rcpp::export]]
// Rcpp::List sample_tau2_arma(const arma::uvec& gc_tracker, const arma::mat& X_s, const arma::mat& A,
//                             const arma::mat& Q_phi, const arma::mat& W, const arma::colvec& G, const arma::vec& beta,
//                             const double& gamma2, const double& tau2, const double& c_tau, const double& d_tau,
//                             const double& sigma2_tau, const double& max_tau2){
//   
//   double tau2_star = R::rnorm(tau2, sqrt(sigma2_tau));
//   double accept_reject = 0.0;
//   if(tau2_star < 0 || tau2_star > max_tau2){
//     return(Rcpp::List::create(Named("tau2") = tau2,
//                               Named("accept_reject") = accept_reject));
//   }else {
//     double log_r = (-c_tau - 1)*(log(tau2_star) - log(tau2)) - d_tau*(1.0/tau2_star - 1.0/tau2) - 
//       0.5 * as_scalar(trans(G - A*X_s*beta)*inv_sympd(gamma2*Q_phi + tau2_star*W)*(G - A*X_s*beta) - 
//       trans(G - A*X_s*beta)*inv_sympd(gamma2*Q_phi + tau2*W)*(G - A*X_s*beta)) -
//       0.5 * (log(det(gamma2 * Q_phi + tau2_star * W)) - log(det(gamma2 * Q_phi + tau2 * W)));
//     // Rcout << "log_r val = " << log_r << "\n";
//     if(log(R::runif(0,1)) < log_r){
//       accept_reject = 1.0;
//     }else {
//       tau2_star = tau2;
//     }
//     
//   }
//   
//   return(Rcpp::List::create(Named("tau2") = tau2_star,
//                             Named("accept_reject") = accept_reject));
// }
// 
// 
// //[[Rcpp::export]]
// Rcpp::List sample_gamma2_arma(const arma::uvec& gc_tracker, const arma::mat& X_s, const arma::mat& A,
//                               const arma::mat& Q_phi, const arma::mat& W, const arma::colvec& G, const arma::vec& beta,
//                               const double& gamma2, const double& tau2, const double& c_gamma, const double& d_gamma,
//                               const double& sigma2_gamma, const double& max_gamma2){
//   
//   double gamma2_star = R::rnorm(gamma2, sqrt(sigma2_gamma));
//   double accept_reject = 0.0;
//   if(gamma2_star < 0 || gamma2_star > max_gamma2){
//     return(Rcpp::List::create(Named("gamma2") = gamma2,
//                               Named("accept_reject") = accept_reject));
//   }else {
//     double log_r = (-c_gamma - 1)*(log(gamma2_star) - log(gamma2)) - d_gamma*(1.0/gamma2_star - 1.0/gamma2) - 
//       0.5 * as_scalar(trans(G - A*X_s*beta)*inv_sympd(gamma2_star*Q_phi + tau2*W)*(G - A*X_s*beta) - 
//       trans(G - A*X_s*beta)*inv_sympd(gamma2*Q_phi + tau2*W)*(G - A*X_s*beta)) - 
//       0.5 * (log(det(gamma2_star * Q_phi + tau2 * W)) - log(det(gamma2 * Q_phi + tau2 * W)));
//     // Rcout << "log_r val = " << log_r << "\n";
//     if(log(R::runif(0,1)) < log_r){
//       accept_reject = 1.0;
//     }else {
//       gamma2_star = gamma2;
//     }
//     
//   }
//   
//   return(Rcpp::List::create(Named("gamma2") = gamma2_star,
//                             Named("accept_reject") = accept_reject));
// }
// 
// 
// 
// 
// 
// 
// //[[Rcpp::export]]
// Rcpp::List run_mcmc_SpRL_growth_arma(const int& n_iter, const arma::mat& Y_mat, const int& N, const double& q,
//                                      const arma::vec& m, const arma::vec& D_bounds, const double& dist, const Rcpp::List& init_vals,
//                                      const Rcpp::List& hyperparameters, const Rcpp::List& rasters, const double& int_length, 
//                                      const std::string& file_name, bool verbose = true){
// 
//   // Calculate necessary quantities
//   int n = Y_mat.n_rows; // Total number of records
//   arma::rowvec mu_D = { (D_bounds(0) + D_bounds(1))/2.0, (D_bounds(2) + D_bounds(3))/2.0 }; // Midpoint of the spatial domain
//   arma::mat Y = Y_mat.cols(0,1); // Observed spatial locations
//   arma::mat V = Y_mat.cols(2,3); // Observed canopy volumes with file designation
//   double v_max = max(V.col(0));
//   arma::vec file_id = Y_mat.col(3); // File IDs for each record
// 
// 
//   // Create storage objects for the parameters
//   arma::vec b_out(n_iter + 1); // Controls how quickly growth mean approaches asymptote
//   arma::vec ar_b(n_iter);
//   arma::mat beta_out(n_iter + 1, rasters.size() + 1); // Linear coefficients for covariates w/ intercept
//   arma::vec phi_out(n_iter + 1); // Interaction distance for spatial covariance of growths
//   arma::vec ar_phi(n_iter);
//   arma::vec gamma2_out(n_iter + 1); // Spatial variance for growths
//   arma::vec ar_gamma2(n_iter);
//   arma::vec tau2_out(n_iter + 1); // Measurement error variance for growths
//   arma::vec ar_tau2(n_iter);
// 
// 
//   // Hyperparameter values
//   double u = hyperparameters["u"]; // b
//   double w = hyperparameters["w"];
//   double sigma2_b_min = as<double>(wrap(hyperparameters["sigma2_b_min"]));
//   double sigma2_b_rate = as<double>(wrap(hyperparameters["sigma2_b_rate"]));
//   arma::mat Sigma_beta_inv = arma::inv_sympd(as<arma::mat>(wrap(hyperparameters["Sigma_beta"]))); // Beta
//   double r = hyperparameters["r"]; // phi
//   double sigma2_phi_min = as<double>(wrap(hyperparameters["sigma2_phi_min"]));
//   double sigma2_phi_rate = as<double>(wrap(hyperparameters["sigma2_phi_rate"]));
//   double c_gamma = hyperparameters["c_gamma"]; // gamma2
//   double d_gamma = hyperparameters["d_gamma"];
//   double max_gamma2 = hyperparameters["max_gamma2"];
//   double sigma2_gamma_min = as<double>(wrap(hyperparameters["sigma2_gamma_min"]));
//   double sigma2_gamma_rate = as<double>(wrap(hyperparameters["sigma2_gamma_rate"]));
//   double c_tau = hyperparameters["c_tau"]; // tau2
//   double d_tau = hyperparameters["d_tau"];
//   double sigma2_tau_min = as<double>(wrap(hyperparameters["sigma2_tau_min"]));
//   double sigma2_tau_rate = as<double>(wrap(hyperparameters["sigma2_tau_rate"]));
//   double max_tau2 = hyperparameters["max_tau2"];
// 
// 
//   // Initialize the parameter values for the sampler
//   arma::uvec lambda = init_vals["lambda"];
//   arma::mat s = init_vals["s"];
//   b_out(0) = init_vals["b"];
//   double sigma2_b = as<double>(wrap(init_vals["sigma2_b"]));
//   arma::rowvec beta_init = init_vals["beta"];
//   beta_out.row(0) = beta_init;
//   phi_out(0) = init_vals["phi"];
//   double sigma2_phi = init_vals["sigma2_phi"];
//   gamma2_out(0) = init_vals["gamma2"];
//   double sigma2_gamma = as<double>(wrap(init_vals["sigma2_gamma"]));
//   tau2_out(0) = init_vals["tau2"];
//   double sigma2_tau = as<double>(wrap(init_vals["sigma2_tau"]));
// 
// 
//   // Create the linkage set C and the growth cluster index and initialize relevant storage
//   List C_lambda = create_clusters(lambda);
//   arma::uvec gc_tracker = growth_cluster_check(C_lambda, V, q);
//   Rcpp::List update_b(2);
//   Rcpp::List update_phi(2);
//   Rcpp::List update_tau2(2);
//   Rcpp::List update_gamma2(2);
//   arma::mat X_s = update_covars_arma(s, rasters);
//   X_s = X_s.rows(find(gc_tracker == 1));
//   arma::colvec G = calculate_growth(C_lambda, gc_tracker, V, int_length);
//   arma::mat W = calculate_W(C_lambda, V, gc_tracker, v_max);
// 
// 
//   if(verbose){
//     Rcout << "Sampler initialized. \n";
//   }
// 
//   // arma::wall_clock timer;
//   arma::wall_clock timer_full;
//   timer_full.tic();
// 
//   for(arma::uword i = 0; i < n_iter; i++){
// 
//     // timer.tic();
//     arma::mat A = calculate_A(C_lambda, V, gc_tracker, b_out(i));
//     arma::mat Q_phi = spatial_covariance(s, gc_tracker, phi_out(i));
// 
// 
//     // Sample Beta
//     beta_out.row(i+1) = sample_beta_arma(X_s, A, W, Q_phi, G, Sigma_beta_inv, gamma2_out(i), tau2_out(i));
// 
// 
//     // Sample b
//     update_b = sample_b_arma(C_lambda, V, gc_tracker, X_s, A, Q_phi, W, G, as<arma::vec>(wrap(beta_out.row(i+1))),
//                              gamma2_out(i), tau2_out(i), u, w, b_out(i), sigma2_b);
//     b_out(i+1) = update_b["b"];
//     ar_b(i) = update_b["accept_reject"];
//     if(ar_b(i) == 1){
//       A = calculate_A(C_lambda, V, gc_tracker, b_out(i+1));
//     }
// 
// 
//     // Sample phi
//     update_phi = sample_phi_arma(C_lambda, V, gc_tracker, X_s, A, Q_phi, W, G, as<arma::vec>(wrap(beta_out.row(i+1))),
//                                  gamma2_out(i), tau2_out(i), s, r, phi_out(i), sigma2_phi);
//     phi_out(i+1) = update_phi["phi"];
//     ar_phi(i) = update_phi["accept_reject"];
//     if(ar_phi(i) == 1){
//       Q_phi = spatial_covariance(s, gc_tracker, phi_out(i+1));
//     }
// 
// 
//     // Sample tau2
//     update_tau2 = sample_tau2_arma(gc_tracker, X_s, A, Q_phi, W, G, as<arma::vec>(wrap(beta_out.row(i+1))),
//                                    gamma2_out(i), tau2_out(i), c_tau, d_tau, sigma2_tau, max_tau2);
//     tau2_out(i+1) = update_tau2["tau2"];
//     ar_tau2(i) = update_tau2["accept_reject"];
// 
// 
//     // Sample gamma2
//     update_gamma2 = sample_gamma2_arma(gc_tracker, X_s, A, Q_phi, W, G, as<arma::vec>(wrap(beta_out.row(i+1))),
//                                        gamma2_out(i), tau2_out(i), c_gamma, d_gamma, sigma2_gamma, max_gamma2);
//     gamma2_out(i+1) = update_gamma2["gamma2"];
//     ar_gamma2(i) = update_gamma2["accept_reject"];
// 
// 
//     // Update MH proposal variances
//     if(i > 0){
//       if(i % 50 == 0){
//         sigma2_b = update_sigma2_tune_arma(sigma2_b, ar_b, i, sigma2_b_min, sigma2_b_rate);
//         sigma2_phi = update_sigma2_tune_arma(sigma2_phi, ar_phi, i, sigma2_phi_min, sigma2_phi_rate);
//         sigma2_tau = update_sigma2_tune_arma(sigma2_tau, ar_tau2, i, sigma2_tau_min, sigma2_tau_rate);
//         sigma2_gamma = update_sigma2_tune_arma(sigma2_gamma, ar_gamma2, i, sigma2_gamma_min, sigma2_gamma_rate);
//         // save_rds(Rcpp::List::create(Named("beta") = beta_out,
//         //                             Named("b") = b_out,
//         //                             Named("ar_b") = ar_b,
//         //                             Named("phi") = phi_out,
//         //                             Named("ar_phi") = ar_phi,
//         //                             Named("tau2") = tau2_out,
//         //                             Named("ar_tau2") = ar_tau2,
//         //                             Named("gamma2") = gamma2_out,
//         //                             Named("ar_gamma2") = ar_gamma2,
//         //                             Named("MH_sigmas") = arma::vec({sigma2_b, sigma2_phi, sigma2_tau, sigma2_gamma})),
//         //                             file_name);
//       }
//     }
// 
//     if(verbose & (i % 100 == 0)){
//       Rcout << i+1 << " / " << n_iter << " iterations complete. \n";
//     }
// 
//     // double n_secs = timer.toc();
//     // Rcout << "Number of seconds: " << n_secs << " for full sampler iteration \n";
//   }
// 
//   double n_secs_full = timer_full.toc();
//   Rcout << "Number of minutes: " << n_secs_full / 60.0 << " to complete full sampler \n";
// 
//   return(Rcpp::List::create(Named("beta") = beta_out,
//                             Named("b") = b_out,
//                             Named("ar_b") = ar_b,
//                             Named("phi") = phi_out,
//                             Named("ar_phi") = ar_phi,
//                             Named("tau2") = tau2_out,
//                             Named("ar_tau2") = ar_tau2,
//                             Named("gamma2") = gamma2_out,
//                             Named("ar_gamma2") = ar_gamma2,
//                             Named("MH_sigmas") = arma::vec({sigma2_b, sigma2_phi, sigma2_tau, sigma2_gamma})));
// 
// }
// 
// 
// 
// 
// 
// // Non-spatial growth model
// 
// //[[Rcpp::export]]
// arma::rowvec sample_beta_arma_ns(const arma::mat& X_s, const arma::mat& A, const arma::mat& W, 
//                                  const arma::colvec& G, const arma::mat& Sigma_beta_inv, const double& tau2){
//   
//   arma::vec mu = inv_sympd(X_s.t()*A*inv_sympd(tau2*W)*A*X_s + Sigma_beta_inv)*X_s.t()*A*inv_sympd(tau2*W)*G;
//   arma::mat cov = inv_sympd(X_s.t()*A*inv_sympd(tau2*W)*A*X_s + Sigma_beta_inv);
//   
//   arma::vec beta = arma::mvnrnd(mu, cov);
//   
//   return(as<arma::rowvec>(wrap(beta)));
// }
// 
// 
// //[[Rcpp::export]]
// Rcpp::List sample_b_arma_ns(const Rcpp::List& C_lambda, const arma::mat& V, const arma::uvec& gc_tracker,
//                             const arma::mat& X_s, const arma::mat& A, const arma::mat& W,
//                             const arma::colvec& G, const arma::vec& beta, const double& tau2,
//                             const double& u, const double& w, const double& b, const double& sigma2_b){
//   
//   double b_star = R::rnorm(b, sqrt(sigma2_b));
//   double accept_reject = 0.0;
//   if(b_star < 0){
//     return(Rcpp::List::create(Named("b") = b,
//                               Named("accept_reject") = accept_reject));
//   }else {
//     arma::mat new_A = calculate_A(C_lambda, V, gc_tracker, b_star); 
//     double log_r = (u - 1)*(log(b_star) - log(b)) - w*(b_star - b) -
//       0.5*as_scalar(trans(G - new_A*X_s*beta)*inv_sympd(tau2*W)*(G - new_A*X_s*beta) - 
//       trans(G - A*X_s*beta)*inv_sympd(tau2*W)*(G - A*X_s*beta));
//     if(log(R::runif(0,1)) < log_r){
//       accept_reject = 1.0;
//     }else {
//       b_star = b;
//     }
//     
//   }
//   
//   return(Rcpp::List::create(Named("b") = b_star,
//                             Named("accept_reject") = accept_reject));
// } 
// 
// 
// 
// // //[[Rcpp::export]]
// // double sample_tau2_arma_ns(const arma::uvec& gc_tracker, const arma::mat& X_s, const arma::mat& A,
// //                            const arma::mat& W, const arma::colvec& G, const arma::vec& beta,
// //                            const double& tau2, const double& c_tau, const double& d_tau,
// //                            const double& max_tau2){
// //   
// //   double tau2_star = 1.0 / R::rgamma(G.n_elem / 2.0 + c_tau, 1.0 / (.5 * as_scalar(trans(G - A*X_s*beta)*inv_sympd(W)*(G - A*X_s*beta)) + d_tau));
// //   // Rcout << "tau2_star = " << tau2_star << "\n";
// //   if(tau2_star < 0 || tau2_star > max_tau2){
// //     return(tau2);
// //   }else {
// //     return(tau2_star);
// //   }
// //   
// // }
// 
// //[[Rcpp::export]]
// double sample_tau2_arma_ns(const arma::uvec& gc_tracker, const arma::mat& X_s, const arma::mat& A,
//                            const arma::mat& W, const arma::colvec& G, const arma::vec& beta,
//                            const double& tau2, const double& c_tau, const double& d_tau,
//                            const double& max_tau2){
//   
//   arma::uvec gc_index = find(gc_tracker == 1);
//   double c_tau_star = c_tau + gc_index.n_elem/2.0;
//   // Rcout << "Value of gc_index.n_elem = " << gc_index.n_elem << " and /2 = " << gc_index.n_elem/2.0 << "\n";
//   double d_tau_star = d_tau + 0.5 * as_scalar((trans(G - A*X_s*beta)*inv_sympd(W)*(G - A*X_s*beta)));
//   double tau2_star = 1.0 / R::rgamma(c_tau_star, 1.0 / d_tau_star);
//   if(tau2_star < max_tau2){
//     return(tau2_star);
//   }else {
//     return(tau2);
//   }
//   
// }
// 
// //[[Rcpp::export]]
// Rcpp::List run_mcmc_SpRL_growth_arma_ns(const int& n_iter, const arma::mat& Y_mat, const int& N, const double& q,
//                                         const arma::vec& m, const arma::vec& D_bounds, const double& dist, const Rcpp::List& init_vals,
//                                         const Rcpp::List& hyperparameters, const Rcpp::List& rasters, const double& int_length, 
//                                         const std::string& file_name, bool verbose = true){
//   
//   // Calculate necessary quantities
//   int n = Y_mat.n_rows; // Total number of records
//   arma::rowvec mu_D = { (D_bounds(0) + D_bounds(1))/2.0, (D_bounds(2) + D_bounds(3))/2.0 }; // Midpoint of the spatial domain
//   arma::mat Y = Y_mat.cols(0,1); // Observed spatial locations
//   arma::mat V = Y_mat.cols(2,3); // Observed canopy volumes with file designation
//   double v_max = max(V.col(0));
//   arma::vec file_id = Y_mat.col(3); // File IDs for each record
//   
//   
//   // Create storage objects for the parameters
//   arma::vec b_out(n_iter + 1); // Controls how quickly growth mean approaches asymptote
//   arma::vec ar_b(n_iter);
//   arma::mat beta_out(n_iter + 1, rasters.size() + 1); // Linear coefficients for covariates w/ intercept
//   arma::vec tau2_out(n_iter + 1); // Measurement error variance for growths
//   
//   
//   // Hyperparameter values
//   double u = hyperparameters["u"]; // b
//   double w = hyperparameters["w"];
//   double sigma2_b_min = as<double>(wrap(hyperparameters["sigma2_b_min"]));
//   double sigma2_b_rate = as<double>(wrap(hyperparameters["sigma2_b_rate"]));
//   arma::mat Sigma_beta_inv = arma::inv_sympd(as<arma::mat>(wrap(hyperparameters["Sigma_beta"]))); // Beta
//   double c_tau = hyperparameters["c_tau"]; // tau2
//   double d_tau = hyperparameters["d_tau"];
//   double max_tau2 = hyperparameters["max_tau2"];
//   
//   
//   // Initialize the parameter values for the sampler
//   arma::uvec lambda = init_vals["lambda"];
//   arma::mat s = init_vals["s"];
//   b_out(0) = init_vals["b"];
//   double sigma2_b = as<double>(wrap(init_vals["sigma2_b"]));
//   arma::rowvec beta_init = init_vals["beta"];
//   beta_out.row(0) = beta_init;
//   tau2_out(0) = init_vals["tau2"];
//   
//   
//   // Create the linkage set C and the growth cluster index and initialize relevant storage
//   List C_lambda = create_clusters(lambda);
//   arma::uvec gc_tracker = growth_cluster_check(C_lambda, V, q);
//   Rcpp::List update_b(2);
//   arma::mat X_s = update_covars_arma(s, rasters);
//   // Rcout << "X_s size " << X_s.n_rows << "\n";
//   X_s = X_s.rows(find(gc_tracker == 1));
//   // Rcout << "X_s size " << X_s.n_rows << "\n";
//   arma::colvec G = calculate_growth(C_lambda, gc_tracker, V, int_length);
//   // Rcout << "G size " << G.n_elem << "\n";
//   // arma::mat W = calculate_W(C_lambda, V, gc_tracker, v_max);
//   arma::mat W(G.n_elem, G.n_elem, arma::fill::eye);
//   
//   if(verbose){
//     Rcout << "Sampler initialized. \n";
//   }
//   
//   // arma::wall_clock timer;
//   arma::wall_clock timer_full;
//   timer_full.tic();
//   
//   for(arma::uword i = 0; i < n_iter; i++){
//     
//     // timer.tic();
//     arma::mat A = calculate_A(C_lambda, V, gc_tracker, b_out(i));
//     
//     
//     // Sample Beta
//     beta_out.row(i+1) = sample_beta_arma_ns(X_s, A, W, G, Sigma_beta_inv, tau2_out(i));
//     
//     
//     // Sample b
//     update_b = sample_b_arma_ns(C_lambda, V, gc_tracker, X_s, A, W, G, as<arma::vec>(wrap(beta_out.row(i+1))),
//                                 tau2_out(i), u, w, b_out(i), sigma2_b);
//     b_out(i+1) = update_b["b"];
//     ar_b(i) = update_b["accept_reject"];
//     if(ar_b(i) == 1){
//       A = calculate_A(C_lambda, V, gc_tracker, b_out(i+1));
//     }
//     
//     
//     // Sample tau2
//     tau2_out(i+1) = sample_tau2_arma_ns(gc_tracker, X_s, A, W, G, as<arma::vec>(wrap(beta_out.row(i+1))),
//                                         tau2_out(i), c_tau, d_tau, max_tau2);
//     
//     
//     // Update MH proposal variances
//     if(i > 0){
//       if(i % 50 == 0){
//         sigma2_b = update_sigma2_tune_arma(sigma2_b, ar_b, i, sigma2_b_min, sigma2_b_rate);
//         // Rcout << "Acceptance rate = " << mean(ar_b(arma::span((i-49), i))) << " and new sigma2_b = " << sigma2_b << "\n";
//         // save_rds(Rcpp::List::create(Named("beta") = beta_out,
//         //                             Named("b") = b_out,
//         //                             Named("ar_b") = ar_b,
//         //                             Named("phi") = phi_out,
//         //                             Named("ar_phi") = ar_phi,
//         //                             Named("tau2") = tau2_out,
//         //                             Named("ar_tau2") = ar_tau2,
//         //                             Named("gamma2") = gamma2_out,
//         //                             Named("ar_gamma2") = ar_gamma2,
//         //                             Named("MH_sigmas") = arma::vec({sigma2_b, sigma2_phi, sigma2_tau, sigma2_gamma})),
//         //                             file_name);
//       }
//     }
//     
//     if(verbose & (i % 100 == 0)){
//       Rcout << i+1 << " / " << n_iter << " iterations complete. \n";
//     }
//     
//     // double n_secs = timer.toc();
//     // Rcout << "Number of seconds: " << n_secs << " for full sampler iteration \n";
//   }
//   
//   double n_secs_full = timer_full.toc();
//   Rcout << "Number of minutes: " << n_secs_full / 60.0 << " to complete full sampler \n";
//   
//   return(Rcpp::List::create(Named("beta") = beta_out,
//                             Named("b") = b_out,
//                             Named("ar_b") = ar_b,
//                             Named("tau2") = tau2_out,
//                             Named("MH_sigmas") = arma::vec({sigma2_b})));
//   
// }

