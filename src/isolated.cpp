// Copyright (c) 2018-2023, University of Bremen, M. Farzalipour Tabriz
// Copyrights licensed under the 2-Clause BSD License.
// See the accompanying LICENSE.txt file for terms.

#include "isolated.hpp"

std::vector<double> nonlinear_fit(const double &opt_tol,
                                  nonlinear_fit_data &fit_data) {
  auto log = spdlog::get("loggers");
  double fit_MSE = 0;
  std::vector<double> fit_parameters = {1, 1, 1, 1};
  const auto opt_algorithm = nlopt::LN_COBYLA;

  nlopt::opt opt(opt_algorithm, 4);

  opt.set_min_objective(fit_eval, &fit_data);
  opt.set_xtol_rel(opt_tol); // tolerance for error value

  try {
    opt.optimize(fit_parameters, fit_MSE);
  } catch (const std::exception &e) {
    log->error("Nonlinear fitting failed: {}", std::string(e.what()));
  }

  return fit_parameters;
}

double fit_eval(const std::vector<double> &c, std::vector<double> &grad,
                void *data) {
  const auto d = static_cast<const nonlinear_fit_data *>(data);
  const arma::rowvec scales = d->sizes;
  const arma::rowvec energies = d->energies;
  const auto madelung_term = d->madelung_term;
  const arma::rowvec model_energies =
      c.at(0) + c.at(1) * scales + c.at(2) * arma::square(scales) +
      (c.at(1) - madelung_term) / c.at(3) * arma::exp(-c.at(3) * scales);
  const auto fit_MSE = arma::accu(arma::square(energies - model_energies));
  return fit_MSE;
}
