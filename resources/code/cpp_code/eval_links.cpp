#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector get_tp_fp_fn(IntegerVector z, IntegerVector true_id) {
  int n = z.length();
  IntegerVector res(3);
  
  int tp, fp, fn;
  tp = 0;
  fp = 0;
  fn = 0;
  
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++) {
      if(i < j) {
        bool est_link = (z(i) == z(j));
        bool true_link = (true_id(i) == true_id(j));
        
        if(est_link) {
          if(true_link) {
            tp++;
          } else{
            fp++;
          }
        } else {
          if(true_link) {
            fn++;
          }
        }
      }
    }
  }
  res(0) = tp;
  res(1) = fp;
  res(2) = fn;
  return(res);
}

// [[Rcpp::export]]
List eval_links(IntegerMatrix z, IntegerVector true_id) {
  List res;
  int M = z.nrow(); // first row is init
  int n = z.ncol();
  
  IntegerVector tp_fp_fn(3);
  int tp;
  int fp;
  int fn;

  NumericVector precision(M);
  NumericVector recall(M);
  NumericVector f1(M);

  // loop through iterations
  for(int m = 0; m < M; m++) {
    IntegerVector z_m = z(m, _);

    tp_fp_fn = get_tp_fp_fn(z_m, true_id);
    tp = tp_fp_fn(0);
    fp = tp_fp_fn(1);
    fn = tp_fp_fn(2);

    precision(m) = 1.0 * tp / (tp + fp);
    recall(m) = 1.0 * tp / (tp + fn);
    f1(m) = 2.0 * tp / (2.0 * tp + fn + fp);
  }

  res["precision"] = precision;
  res["recall"] = recall;
  res["F1_Score"] = f1;
  return(res);
}

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector get_prec_rec_rcpp(IntegerVector lambda, IntegerVector true_id) {
  int n = lambda.size();
  double tp = 0, tp_fn = 0, tp_fp = 0;
  
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++) {
      bool true_link = true_id[i] == true_id[j];
      bool est_link = lambda[i] == lambda[j];
      tp += true_link && est_link;
      tp_fn += true_link;
      tp_fp += est_link;
    }
  }
  
  double precision = tp/tp_fp;
  double recall = tp/tp_fn;
  double F1_score = (2 * precision * recall) / (precision + recall);
  
  return NumericVector::create(
    _["precision"] = precision,
    _["recall"] = recall,
    _["F1_score"] = F1_score
  );
}


// [[Rcpp::export]]
int ss_growth2(NumericVector lambda, DataFrame V) {
  IntegerVector files = V["file"];
  NumericVector sizes = V["size"];
  int n = lambda.size();
  int N = max(lambda) + 1;
  int growth_clusters = 0;
  
  for (int i = 0; i < N; i++) {
    double v_min = 0;
    double v_max = 0;
    int min_file = INT_MAX;
    int max_file = INT_MIN;
    int cluster_size = 0;
    
    for (int j = 0; j < n; j++) {
      if (lambda[j] == i) {
        cluster_size++;
        if (files[j] < min_file) {
          min_file = files[j];
          v_min = sizes[j];
        } else if (files[j] == min_file) {
          v_min += sizes[j];
        }
        if (files[j] > max_file) {
          max_file = files[j];
          v_max = sizes[j];
        } else if (files[j] == max_file) {
          v_max += sizes[j];
        }
      }
    }
    
    if (cluster_size >= 2 && min_file != max_file && v_max > 0.8 * v_min) {
      growth_clusters++;
    }
  }
  
  return growth_clusters;
}


