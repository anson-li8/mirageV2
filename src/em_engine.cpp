#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace Rcpp;

// Stable log(exp(a) + exp(b))
inline double log_sum_exp(double a, double b) {
  if (a == R_NegInf) return b;
  if (b == R_NegInf) return a;
  double m = std::max(a, b);
  return m + std::log(std::exp(a - m) + std::exp(b - m));
}

// log((1-eta) + eta * exp(log_bf))
inline double log_mixture(double log_bf, double eta) {
  if (eta <= 0.0) return 0.0;
  if (eta >= 1.0) return log_bf;
  return log_sum_exp(std::log1p(-eta), std::log(eta) + log_bf);
}

// [[Rcpp::export]]
List em_mirage_cpp(IntegerVector category,      // 0-indexed category per variant
                   IntegerVector gene_index,   // 0-indexed gene per variant
                   NumericVector log_var_bf,
                   int n_genes,
                   int n_categories,
                   double delta_init,
                   NumericVector eta_init,
                   bool estimate_delta,
                   bool estimate_eta,
                   NumericVector fixed_eta,    // length n_categories (used if !estimate_eta)
                   int max_iter,
                   double tol) {

  int n_variants = log_var_bf.size();

  // Init parameters
  double delta = delta_init;
  std::vector<double> eta(n_categories);
  for (int g = 0; g < n_categories; ++g) {
    eta[g] = estimate_eta ? eta_init[g] : fixed_eta[g];
  }

  // Precompute: variants per gene
  std::vector<std::vector<int>> gene_variants(n_genes);
  for (int j = 0; j < n_variants; ++j) {
    gene_variants[gene_index[j]].push_back(j);
  }

  // Precompute: count of category g variants in gene i
  // counts[i * n_categories + g]
  std::vector<int> counts(n_genes * n_categories, 0);
  for (int j = 0; j < n_variants; ++j) {
    counts[gene_index[j] * n_categories + category[j]]++;
  }

  // Working arrays
  std::vector<double> gene_log_bf(n_genes);
  std::vector<double> EUi(n_genes);
  std::vector<double> UiZij(n_variants);

  // History
  NumericMatrix eta_history(max_iter, n_categories);
  NumericVector delta_history(max_iter);

  // Store initial
  for (int g = 0; g < n_categories; ++g) eta_history(0, g) = eta[g];
  delta_history[0] = delta;

  bool converged = false;
  int final_iter = 1;

  for (int iter = 1; iter < max_iter; ++iter) {
    std::vector<double> eta_prev = eta;
    double delta_prev = delta;

    // Gene-level log BFs
    for (int i = 0; i < n_genes; ++i) {
      double lb = 0.0;
      for (int j : gene_variants[i]) {
        lb += log_mixture(log_var_bf[j], eta_prev[category[j]]);
      }
      gene_log_bf[i] = lb;
    }

    // E-step: EUi
    double log_delta = std::log(delta_prev);
    double log_1m_delta = std::log1p(-delta_prev);
    for (int i = 0; i < n_genes; ++i) {
      double log_numer = log_delta + gene_log_bf[i];
      double log_denom = log_sum_exp(log_numer, log_1m_delta);
      EUi[i] = std::exp(log_numer - log_denom);
    }

    // E-step: UiZij
    for (int j = 0; j < n_variants; ++j) {
      int g = category[j];
      int i = gene_index[j];
      double log_eta = std::log(eta_prev[g]);
      double log_1m_eta = std::log1p(-eta_prev[g]);
      double a = log_eta + log_var_bf[j];
      double b = log_1m_eta;
      double log_P_j = a - log_sum_exp(a, b);
      UiZij[j] = EUi[i] * std::exp(log_P_j);
    }

    // M-step: delta
    if (estimate_delta) {
      double sum_eui = 0.0;
      for (int i = 0; i < n_genes; ++i) sum_eui += EUi[i];
      delta = sum_eui / n_genes;
      delta = std::max(delta, 1e-300);
      delta = std::min(delta, 1.0 - 1e-15);
    }

    // M-step: eta
    if (estimate_eta) {
      for (int g = 0; g < n_categories; ++g) {
        double numerator = 0.0;
        for (int j = 0; j < n_variants; ++j) {
          if (category[j] == g) numerator += UiZij[j];
        }
        if (numerator <= 0.0) {
          eta[g] = 0.0;
        } else {
          double denominator = 0.0;
          for (int i = 0; i < n_genes; ++i) {
            denominator += counts[i * n_categories + g] * EUi[i];
          }
          eta[g] = numerator / denominator;
          eta[g] = std::min(eta[g], 1.0 - 1e-15);
        }
      }
    }

    // Store history
    for (int g = 0; g < n_categories; ++g) eta_history(iter, g) = eta[g];
    delta_history[iter] = delta;

    // Convergence check
    double diff = 0.0;
    for (int g = 0; g < n_categories; ++g) diff += std::abs(eta[g] - eta_prev[g]);

    final_iter = iter + 1;
    if (diff < tol) {
      converged = true;
      break;
    }
  }

  // Final posteriors
  // Recompute with final eta/delta
  for (int i = 0; i < n_genes; ++i) {
    double lb = 0.0;
    for (int j : gene_variants[i]) {
      lb += log_mixture(log_var_bf[j], eta[category[j]]);
    }
    gene_log_bf[i] = lb;
  }

  NumericVector post_prob(n_variants);
  NumericVector gene_log_bf_out(n_genes);
  NumericVector EUi_out(n_genes);

  double log_delta_f = std::log(delta);
  double log_1m_delta_f = std::log1p(-delta);
  for (int i = 0; i < n_genes; ++i) {
    double log_numer = log_delta_f + gene_log_bf[i];
    double log_denom = log_sum_exp(log_numer, log_1m_delta_f);
    EUi_out[i] = std::exp(log_numer - log_denom);
    gene_log_bf_out[i] = gene_log_bf[i];
  }

  for (int j = 0; j < n_variants; ++j) {
    int g = category[j];
    int i = gene_index[j];
    double log_eta = std::log(eta[g]);
    double log_1m_eta = std::log1p(-eta[g]);
    double a = log_eta + log_var_bf[j];
    double b = log_1m_eta;
    double log_P_j = a - log_sum_exp(a, b);
    post_prob[j] = EUi_out[i] * std::exp(log_P_j);
  }

  // Trim history
  NumericMatrix eta_hist_out = eta_history(Range(0, final_iter - 1), Range(0, n_categories - 1));
  NumericVector delta_hist_out = delta_history[Range(0, final_iter - 1)];

  NumericVector eta_out(n_categories);
  for (int g = 0; g < n_categories; ++g) eta_out[g] = eta[g];

  return List::create(
    _["delta.est"] = delta,
    _["eta.est"] = eta_out,
    _["gene.log.bf"] = gene_log_bf_out,
    _["post.prob"] = post_prob,
    _["EUi"] = EUi_out,
    _["eta.history"] = eta_hist_out,
    _["delta.history"] = delta_hist_out,
    _["converged"] = converged,
    _["n.iter"] = final_iter
  );
}
