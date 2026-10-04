#include <Rcpp.h>
#include <cmath>
#include <algorithm>

using namespace Rcpp;

// Helper to compute log integrand
inline double log_integrand(double aa, int x, int T, double N1, double N0,
                            double shape, double scale) {
  if (aa <= 0.0) return R_NegInf;
  double p = (aa * N1) / (aa * N1 + N0);
  if (p < 0.0) p = 0.0;
  if (p > 1.0) p = 1.0;
  double log_binom = R::dbinom(x, T, p, 1);
  double log_gamma = R::dgamma(aa, shape, scale, 1);
  return log_binom + log_gamma;
}

double compute_log_integral(int x, int T, double N1, double N0,
                            double shape, double scale) {
  double a = (shape < 1.0) ? 1e-8 : 0.0;
  double b = 100.0;  // Original mirage limit

  // Composite Simpson's rule with fine grid
  int n = 1000;  // Must be even
  double h = (b - a) / n;

  // Find max for log-sum-exp stability
  double M = R_NegInf;
  for (int i = 0; i <= n; ++i) {
    double aa = a + i * h;
    double val = log_integrand(aa, x, T, N1, N0, shape, scale);
    if (val > M) M = val;
  }
  if (M == R_NegInf) return R_NegInf;

  // Composite Simpson's on shifted integrand
  double sum = 0.0;
  for (int i = 0; i <= n; ++i) {
    double aa = a + i * h;
    double val = std::exp(log_integrand(aa, x, T, N1, N0, shape, scale) - M);
    if (i == 0 || i == n) {
      sum += val;
    } else if (i % 2 == 1) {
      sum += 4.0 * val;
    } else {
      sum += 2.0 * val;
    }
  }
  double integral = (h / 3.0) * sum;

  if (integral <= 0.0) return R_NegInf;
  return M + std::log(integral);
}

// [[Rcpp::export]]
NumericVector calc_log_bf_var_cpp(IntegerVector var_case,
                                  IntegerVector var_contr,
                                  NumericVector gamma_vec,
                                  NumericVector sigma_vec,
                                  double N1,
                                  double N0) {
  int n = var_case.size();
  NumericVector log_bf(n);

  for (int i = 0; i < n; ++i) {
    int x = var_case[i];
    int T = x + var_contr[i];
    if (T == 0) {
      log_bf[i] = 0.0;
      continue;
    }

    double shape = gamma_vec[i] * sigma_vec[i];
    double scale = 1.0 / sigma_vec[i];

    double log_marg1 = compute_log_integral(x, T, N1, N0, shape, scale);
    double p0 = N1 / (N1 + N0);
    double log_marg0 = R::dbinom(x, T, p0, 1);

    log_bf[i] = log_marg1 - log_marg0;
  }
  return log_bf;
}
