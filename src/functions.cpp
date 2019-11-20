#include <vector>
#include <iostream>

#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix estimateNewNegativeScores(NumericMatrix expression, const std::vector<double> &max_scores,
                                        const std::vector<double> &neg_scores) {
  if (neg_scores.size() != expression.nrow() || max_scores.size() != expression.nrow())
    stop("max_scores and neg_scores must have the same length as expression");

  NumericMatrix res(expression.nrow(), expression.ncol());

  for (int j = 0; j < expression.ncol(); ++j) {
    for (int i = 0; i < expression.nrow(); ++i) {
      double t_exp = expression(i, j);
      if (t_exp < -1e-30)
        stop("Expression must be positive");

      if (max_scores[i] <= t_exp) {
        res(i, j) = 1;
        continue;
      }

      res(i, j) = std::max(t_exp / max_scores[i], neg_scores[i]);
    }
  }

  return(res);
}

// [[Rcpp::export]]
NumericMatrix estimatePariwiseNegativeScoreChange(int base_id, NumericMatrix neg_scores, const std::vector<double> &pos_scores) {
  base_id -= 1;
  if (neg_scores.nrow() != pos_scores.size())
    stop("neg_scores and pos_scores must have the same length");

  if (base_id < 0 || base_id >= neg_scores.ncol())
    stop("base_id is out of bounds: " + std::to_string(base_id));

  NumericMatrix score_change(Rcpp::clone(neg_scores));
  NumericVector v = neg_scores(_, base_id);
  std::vector<double> base_neg_scores = as<std::vector<double>>(v);
  for (int j = 0; j < score_change.ncol(); ++j) {
    for (int i = 0; i < score_change.nrow(); ++i) {
      score_change(i, j) = (1 - std::max(score_change(i, j), base_neg_scores[i])) * pos_scores[i];
    }
  }

  return(score_change);
}

// [[Rcpp::export]]
std::vector<double> estimateDNegativeScores(NumericMatrix d_scores, const std::vector<double> &pos_scores,
                                            const std::vector<double> &sum_scores, const std::vector<bool> is_positive) {
  if (sum_scores.size() != d_scores.nrow() || pos_scores.size() != d_scores.nrow() || is_positive.size() != d_scores.nrow())
    stop("sum_scores, pos_scores and is_positive must have the same length as expression");

  // NumericMatrix res(d_scores.nrow(), d_scores.ncol());
  std::vector<double> res(d_scores.ncol(), 0);

  for (int j = 0; j < d_scores.ncol(); ++j) {
    for (int i = 0; i < d_scores.nrow(); ++i) {
      double td = d_scores(i, j);
      double ds = (pos_scores[i] + td) / std::max(sum_scores[i] + td, 1e-30) - (pos_scores[i] / sum_scores[i]);
      if (is_positive[i]) {
        res[j] += ds;
      } else {
        res[j] -= ds;
      }
    }
  }

  return(res);
}
