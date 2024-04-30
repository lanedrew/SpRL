#define ARMA_64BIT_WORD 1
#include <RcppDist.h>
// [[Rcpp::depends(RcppArmadillo, RcppDist, Rcpp, terra, base)]]
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;


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


// [[Rcpp::export]]
double sample_sigma2_arma(const arma::mat& s, const arma::uvec& lambda, const double& c_sigma, const double& d_sigma,
                          const arma::mat& Y, const arma::mat& R_theta, const arma::rowvec& t, const arma::rowvec& mu_D,
                          const arma::vec& m, const arma::vec& range, const double& sigma2_prev){
  
  double n = Y.n_rows;
  arma::mat s_lambda = s.rows(lambda);
  arma::mat s_lambda_m = s_lambda.rows(arma::span(m(0), sum(m) - 1));
  
  s_lambda_m = transform_s_arma(s_lambda_m, R_theta, t, mu_D);
  s_lambda = join_cols(s_lambda.rows(0, m(0) - 1), s_lambda_m);
  
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


//[[Rcpp::export]]
arma::uvec sample_lambda_arma(const arma::uvec& lambda, const arma::mat& s, const double& sigma2, const arma::mat& R_theta,
                              const arma::rowvec& t, const arma::rowvec& mu_D, arma::mat Y, const arma::vec& file_id, const int& N,
                              const int& n, const arma::vec& m, const double& dist, const Rcpp::List& C_lambda){
  
  // Rcout << "inside sample_lambda \n";
  arma::mat s_trans = transform_s_arma(s, R_theta, t, mu_D);
  arma::uvec lambda_k = lambda;
  Rcpp::List new_C_lambda = Rcpp::clone(C_lambda);
  
  // arma::wall_clock timer;
  // timer.tic();
  for(int k = 0; k < n; ++k){
    
    if(file_id(k) == 1){
      arma::uvec latent_candidates = near_latents_arma(s, Y.row(k), dist);
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
      s_means.row(i) = Arma_colSums(new_Y.rows(file_indices)) / s_sizes(i);
    }
  }
  
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
  arma::mat new_Y = arma::join_cols(Y_file1, Y_trans);
  
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


// Function to calculate Larger Neighbor Volume (LNV), Neighborhood Density (ND), and Relative Spacing Index (RSI)
// within a specified radius
// [[Rcpp::export]]
arma::mat calc_comp_indices_arma2(const arma::mat& treeData, const double& radius = 15.0) {
  
  int n = treeData.n_rows;
  arma::mat indices(n, 3); // Initialize matrix for comp metrics
  
  for (int i = 0; i < n; i++) {
    double x1 = treeData(i, 0);
    double y1 = treeData(i, 1);
    double size1 = treeData(i, 2);
    
    double LNV = 0.0; // Initialize LNV
    double AND = 0.0; // Initialize AND
    double ND = 0.0; // Initialize ND
    double dnn = radius; // Initialize DNN to a large value
    double RSI = 0.0; // Initialize RSI
    
    double neighbor_counter = 0.0;
    
    for (int j = 0; j < n; j++) {
      if (i != j) { // Exclude the current tree from calculations
        double x2 = treeData(j, 0);
        double y2 = treeData(j, 1);
        double size2 = treeData(j, 2);
        
        double distance = euclidean_dist_arma(x1, y1, x2, y2);
        
        if (distance <= radius) { // Consider trees within the specified radius
          // Update BAN by adding the basal area of the neighbor
          if(size2 > size1){
            LNV += size2;
          }
          // Update AND if the current neighbor is within the specified radius
          AND += distance;
          neighbor_counter += 1.0;
          
          if(dnn < distance){
            dnn = distance;
          }
        }
      }
    }
    
    if(neighbor_counter > 0.0){
      AND = AND / neighbor_counter;
      ND = neighbor_counter / (arma::datum::pi * pow(radius, 2));
    } else{
      AND = radius;
    }
    
    RSI = dnn / AND;
    
    indices(i, 0) = LNV;
    indices(i, 1) = RSI;
    indices(i, 2) = ND;
  }
  
  return indices;
}