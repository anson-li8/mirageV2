#include <Rcpp.h>
#include <cmath>
#include <algorithm>

using namespace Rcpp;

// Helper to compute log integrand
inline double log_integrand(double aa, int x, int T, double N1, double N0, double shape, double scale) {
  if (aa <= 0.0) return R_NegInf;
  double p = (aa * N1) / (aa * N1 + N0);
  if (p < 0.0) p = 0.0;
  if (p > 1.0) p = 1.0;
  double log_binom = R::dbinom(x, T, p, 1);
  double log_gamma = R::dgamma(aa, shape, scale, 1);
  return log_binom + log_gamma;
}

// Adaptive Simpson's integration
double adaptive_simpson(double a, double b, double fa, double fb, double fm, double S, double eps, int depth,
                        int x, int T, double N1, double N0, double shape, double scale, double M) {
    double m = (a + b) / 2.0;
    double h = (b - a) / 2.0;
    double fml = std::exp(log_integrand((a + m) / 2.0, x, T, N1, N0, shape, scale) - M);
    double fmr = std::exp(log_integrand((m + b) / 2.0, x, T, N1, N0, shape, scale) - M);

    double S_left = (h / 6.0) * (fa + 4.0 * fml + fm);
    double S_right = (h / 6.0) * (fm + 4.0 * fmr + fb);
    double S_new = S_left + S_right;

    if (depth <= 0 || std::abs(S_new - S) <= 15.0 * eps) {
        return S_new + (S_new - S) / 15.0;
    }

    return adaptive_simpson(a, m, fa, fm, fml, S_left, eps/2.0, depth-1, x, T, N1, N0, shape, scale, M) +
           adaptive_simpson(m, b, fm, fb, fmr, S_right, eps/2.0, depth-1, x, T, N1, N0, shape, scale, M);
}

double compute_log_integral(int x, int T, double N1, double N0, double shape, double scale) {
    // In case shape < 1 where Gamma PDF approaches infinity at 0
    double a = (shape < 1.0) ? 1e-8 : 0.0;
    double b = 100.0; // Original mirage truncation limit

    // Find approximate max to use as shift M (prevents overflow)
    double M = R_NegInf;
    int steps = 100;
    for (int i = 0; i <= steps; ++i) {
        double aa = a + i * (b - a) / steps;
        double val = log_integrand(aa, x, T, N1, N0, shape, scale);
        if (val > M) M = val;
    }

    if (M == R_NegInf) return R_NegInf;

    double fa = std::exp(log_integrand(a, x, T, N1, N0, shape, scale) - M);
    double fb = std::exp(log_integrand(b, x, T, N1, N0, shape, scale) - M);
    double fm = std::exp(log_integrand((a+b)/2.0, x, T, N1, N0, shape, scale) - M);

    double S = ((b - a) / 6.0) * (fa + 4.0 * fm + fb);

    double eps = 1e-5;
    int depth = 50;

    double integral = adaptive_simpson(a, b, fa, fb, fm, S, eps, depth, x, T, N1, N0, shape, scale, M);

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
