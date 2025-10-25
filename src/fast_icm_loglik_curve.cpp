#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::depends(RcppEigen)]]

//' Fast ICM Log-likelihood Computation using Eigen (Fixed Version)
//'
//' Computes the log-likelihood for a single curve using precomputed cache from
//' .build_icm_cache(). Uses simplified Eigen operations to avoid compilation warnings.
//'
//' @param yresid Matrix of residuals (P x M) where P is number of time points and M is number of channels
//' @param Uxx Eigenvectors from kernel eigendecomposition (P x P)
//' @param chollist List of Cholesky factors for each eigenvalue block
//' @param logdetsum Sum of log determinants for likelihood computation
//'
//' @return Numeric scalar representing the log-likelihood value
// [[Rcpp::export]]
double fast_icm_loglik_curve_eigen(
  const Eigen::MatrixXd& yresid,
  const Eigen::MatrixXd& Uxx,
  const Rcpp::List& chollist,
  double logdetsum
) {
  int P = yresid.rows();
  int M = yresid.cols();
  int Pch = chollist.size();

  // Compute Ytil = Uxx^T * yresid using simple matrix multiplication
  Eigen::MatrixXd Ytil = Uxx.transpose() * yresid;
  
  double quad = 0.0;
  
  // Process each block
  for (int j = 0; j < Pch; ++j) {
    // Extract Cholesky factor for this block
    Eigen::MatrixXd Lj = Rcpp::as<Eigen::MatrixXd>(chollist[j]);
    
    // Extract row j from Ytil
    Eigen::VectorXd v = Ytil.row(j);
    
    // Solve Lj * w = v using forward substitution (simplified)
    Eigen::VectorXd w = Lj.triangularView<Eigen::Lower>().solve(v);
    
    // Add squared norm to quadratic form
    quad += w.dot(w);  // Use dot product instead of squaredNorm()
  }
  
  // Compute final log-likelihood
  double loglik = -0.5 * (Pch * M * std::log(2.0 * M_PI) + logdetsum + quad);
  return loglik;
}

//' Batch Fast ICM Log-likelihood Computation (Fixed Version)
//'
//' Computes log-likelihoods for multiple curves in a single C++ call for better performance.
//'
//' @param curves List of residual curve matrices
//' @param Uxx Eigenvectors from kernel eigendecomposition (P x P)
//' @param chollist List of Cholesky factors for each eigenvalue block
//' @param logdetsum Sum of log determinants for likelihood computation
//'
//' @return Numeric vector of log-likelihood values
// [[Rcpp::export]]
Rcpp::NumericVector fast_icm_loglik_curves_batch(
  const Rcpp::List& curves,
  const Eigen::MatrixXd& Uxx,
  const Rcpp::List& chollist,
  double logdetsum
) {
  int n_curves = curves.size();
  Rcpp::NumericVector results(n_curves);
  
  for (int i = 0; i < n_curves; ++i) {
    Eigen::MatrixXd yresid = Rcpp::as<Eigen::MatrixXd>(curves[i]);
    results[i] = fast_icm_loglik_curve_eigen(yresid, Uxx, chollist, logdetsum);
  }
  
  return results;
}

// Wavelet filter structure
struct WaveletFilter {
  std::vector<double> low_pass;
  std::vector<double> high_pass;
  int length;
};

// Get wavelet filter coefficients
WaveletFilter get_wavelet_filter(const std::string& wf) {
  WaveletFilter filter;
  
  if (wf == "la8") {
    // Daubechies 8-tap filter
    filter.low_pass = {
      -0.010597401784997278, 0.032883011666982945, 0.030841381835986965,
      -0.18703481171888114, -0.027983769416983849, 0.63088076792959036,
      0.71484657055254153, 0.23037781330885523
    };
    filter.high_pass = {
      -0.23037781330885523, 0.71484657055254153, -0.63088076792959036,
      -0.027983769416983849, 0.18703481171888114, 0.030841381835986965,
      -0.032883011666982945, -0.010597401784997278
    };
  } else if (wf == "haar") {
    // Haar wavelet
    filter.low_pass = {0.7071067811865476, 0.7071067811865476};
    filter.high_pass = {-0.7071067811865476, 0.7071067811865476};
  } else {
    // Default to Haar if unknown
    filter.low_pass = {0.7071067811865476, 0.7071067811865476};
    filter.high_pass = {-0.7071067811865476, 0.7071067811865476};
  }
  
  filter.length = filter.low_pass.size();
  return filter;
}

// Fast 1D wavelet forward transform
// [[Rcpp::export]]
Rcpp::List fast_wavelet_forward_1d(
  const Rcpp::NumericVector& x,
  const std::string& wf = "la8",
  int J = 4,
  const std::string& boundary = "periodic"
) {
  int n = x.size();
  if (n <= 0) {
    return Rcpp::List::create(
      Rcpp::Named("coeff") = Rcpp::NumericVector(),
      Rcpp::Named("map") = Rcpp::List::create()
    );
  }
  
  // Get wavelet filter
  WaveletFilter filter = get_wavelet_filter(wf);
  
  // Initialize coefficients
  std::vector<double> coeff;
  std::vector<int> map;
  
  // Copy input data
  std::vector<double> data(x.begin(), x.end());
  
  // Perform wavelet decomposition
  for (int j = 0; j < J && data.size() > filter.length; ++j) {
    int n_curr = data.size();
    int n_new = n_curr / 2;
    
    if (n_new == 0) break;
    
    std::vector<double> low(n_new);
    std::vector<double> high(n_new);
    
    // Convolution with low-pass and high-pass filters
    for (int i = 0; i < n_new; ++i) {
      double low_sum = 0.0;
      double high_sum = 0.0;
      
      for (int k = 0; k < filter.length; ++k) {
        int idx = (2 * i + k) % n_curr;  // Periodic boundary
        low_sum += filter.low_pass[k] * data[idx];
        high_sum += filter.high_pass[k] * data[idx];
      }
      
      low[i] = low_sum;
      high[i] = high_sum;
    }
    
    // Store high-frequency coefficients
    for (int i = 0; i < n_new; ++i) {
      coeff.push_back(high[i]);
    }
    map.push_back(n_new);
    
    // Continue with low-frequency part
    data = low;
  }
  
  // Store final low-frequency coefficients
  for (size_t i = 0; i < data.size(); ++i) {
    coeff.push_back(data[i]);
  }
  map.push_back(data.size());
  
  // Create R list
  Rcpp::NumericVector r_coeff(coeff.begin(), coeff.end());
  Rcpp::List r_map;
  for (size_t i = 0; i < map.size(); ++i) {
    r_map.push_back(map[i]);
  }
  
  return Rcpp::List::create(
    Rcpp::Named("coeff") = r_coeff,
    Rcpp::Named("map") = r_map
  );
}

// Fast 1D wavelet inverse transform
// [[Rcpp::export]]
Rcpp::NumericVector fast_wavelet_inverse_1d(
  const Rcpp::NumericVector& coeff,
  const Rcpp::List& map
) {
  if (coeff.size() == 0) {
    return Rcpp::NumericVector();
  }
  
  // Get wavelet filter (use la8 as default)
  WaveletFilter filter = get_wavelet_filter("la8");
  
  // Extract coefficients and map
  std::vector<double> coeff_vec(coeff.begin(), coeff.end());
  std::vector<int> map_vec;
  for (int i = 0; i < map.size(); ++i) {
    map_vec.push_back(Rcpp::as<int>(map[i]));
  }
  
  // Reconstruct signal
  std::vector<double> data;
  int coeff_idx = 0;
  
  // Start with the coarsest level
  for (int i = map_vec.size() - 1; i >= 0; --i) {
    int n_level = map_vec[i];
    
    if (i == map_vec.size() - 1) {
      // Coarsest level - just copy coefficients
      for (int j = 0; j < n_level; ++j) {
        data.push_back(coeff_vec[coeff_idx + j]);
      }
    } else {
      // Reconstruct from low and high frequency components
      std::vector<double> low(data.begin(), data.end());
      std::vector<double> high;
      
      // Extract high-frequency coefficients
      for (int j = 0; j < n_level; ++j) {
        high.push_back(coeff_vec[coeff_idx + j]);
      }
      
      // Reconstruct signal
      std::vector<double> reconstructed(2 * n_level);
      for (int j = 0; j < n_level; ++j) {
        for (int k = 0; k < filter.length; ++k) {
          int idx = (2 * j + k) % (2 * n_level);
          if (k < filter.length / 2) {
            reconstructed[idx] += filter.low_pass[k] * low[j];
          } else {
            reconstructed[idx] += filter.high_pass[k] * high[j];
          }
        }
      }
      
      data = reconstructed;
    }
    
    coeff_idx += n_level;
  }
  
  return Rcpp::NumericVector(data.begin(), data.end());
}