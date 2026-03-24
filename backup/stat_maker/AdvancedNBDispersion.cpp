#include "AdvancedNBDispersion.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <iostream>

/**
 * @brief computeMu
 * 
 *  mu_{i,r} = exp( offset[r] + logAbundance[i] + group[r]*logFC[i] )
 * 
 * We'll store them in a 2D array for convenience: mu[i][r].
 */
static std::vector<std::vector<double> >
computeMu(const CountData &data) {
    std::vector<std::vector<double> > mu(data.nGenes, std::vector<double>(data.nReplicates, 0.0));
    for (int i = 0; i < data.nGenes; ++i) {
        for (int r = 0; r < data.nReplicates; ++r) {
            double eta = data.offset[r] + data.logAbundance[i];
            if (data.group[r] == 1) {
                eta += data.logFC[i];
            }
            mu[i][r] = std::exp(eta);
        }
    }
    return mu;
}

/**
 * @brief logNB Likelihood
 *
 * The negative binomial pmf:
 *   P(Y = y) = Gamma(y + 1/alpha) / [Gamma(1/alpha)*Gamma(y+1)] * (1/(1+alpha*mu))^(1/alpha) * (alpha*mu/(1+alpha*mu))^y
 *
 * We define the log-likelihood for a single observation y_{i,r} with mean mu_{i,r}:
 *   ln L(y_{i,r}; alpha)
 *     = lgamma(y_{i,r} + 1/alpha) - lgamma(1/alpha) - lgamma(y_{i,r}+1)
 *       + (1/alpha)*ln[1/(1+alpha*mu_{i,r})]
 *       + y_{i,r}*ln[ alpha*mu_{i,r}/(1+alpha*mu_{i,r}) ]
 *
 * In practice, we skip constant terms w.r.t. alpha (like -lgamma(y_{i,r}+1)).
 */
static double logNB(double y, double mu, double alpha) {
    // We skip -lgamma(y+1) since it doesn't depend on alpha
    // Also watch for alpha->0 edge cases
    if (alpha <= 0.0 || mu <= 0.0) {
        return 0.0; // not strictly correct, but for safety
    }

    double term1 = std::lgamma(y + 1.0 / alpha) - std::lgamma(1.0 / alpha);
    double term2 = (1.0 / alpha) * std::log(1.0 / (1.0 + alpha * mu));
    double term3 = y * std::log((alpha * mu) / (1.0 + alpha * mu));

    return term1 + term2 + term3;
}

/**
 * @brief coxReidAdjustment
 * 
 * The Cox-Reid adjustment involves the log-determinant of the information matrix
 * w.r.t. nuisance parameters. In a simplified form, we might approximate this by:
 *     penalty(alpha) = 0.5 * sum_{i,r} [ log(1 + alpha * mu_{i,r}) ]
 * or other expansions. 
 * 
 * In edgeR, the approach is more nuanced, but here's a commonly referenced simplified version.
 */
static double coxReidAdjustment(const std::vector<std::vector<double> > &mu,
                                double alpha) {
    double sumLog = 0.0;
    // Summation over all i,r
    for (size_t i = 0; i < mu.size(); ++i) {
        for (size_t r = 0; r < mu[i].size(); ++r) {
            double val = 1.0 + alpha * mu[i][r];
            if (val > 0.0) {
                sumLog += std::log(val);
            }
        }
    }
    // Factor 0.5 in many references 
    return 0.5 * sumLog;
}

/**
 * @brief fitCommonDispersion
 *
 * Implementation details:
 * 1) Compute mu[i][r] from the data's offsets, logAbundance, logFC.
 * 2) Define:
 *      L(alpha) = sum_{i,r} logNB(y_{i,r}, mu_{i,r}, alpha)  -  coxReidAdjustment(...) 
 * 3) Numerically differentiate L(alpha) w.r.t. alpha to get L'(alpha) and L''(alpha).
 * 4) Perform Newton steps or a fallback approach if alpha leaves domain.
 * 5) Stop if convergence or maxIter.
 */
double fitCommonDispersion(CountData &data,
                           double alphaStart,
                           int maxIter,
                           double tol) {
    // 1) Precompute mu
    auto mu = computeMu(data);

    // Start alpha
    double alpha = (alphaStart > 0.0) ? alphaStart : 0.1;

    for (int iter = 0; iter < maxIter; ++iter) {
        // Evaluate log-likelihood with Cox-Reid
        double logLik = 0.0;
        // We'll accumulate derivatives wrt alpha: first (d1) and second (d2)
        double d1 = 0.0;
        double d2 = 0.0;

        for (int i = 0; i < data.nGenes; ++i) {
            for (int r = 0; r < data.nReplicates; ++r) {
                double y = data.y[i][r];
                double m = mu[i][r];
                if (m <= 0.0) continue;
                // let's define helper terms
                // We want partial derivatives of logNB wrt alpha
                //   logNB = lgamma(y + 1/alpha) - lgamma(1/alpha)
                //          + (1/alpha)*ln(1/(1+alpha*m)) + y ln(alpha*m/(1+alpha*m))

                // We'll define some partial derivatives:
                //    d/dalpha of logNB( y, m, alpha )
                // We can do direct symbolic derivation or numeric approach.
                // For advanced usage, see: Robinson & Smyth 2008, eqn(5),(6).
                //
                // We'll do direct partial derivatives:
                // derivative wrt alpha:
                //   d/dalpha [lgamma(y + 1/alpha)] 
                //     = digamma(y + 1/alpha)*(-1/alpha^2)
                //   d/dalpha [-lgamma(1/alpha)] 
                //     = -digamma(1/alpha)*(-1/alpha^2)
                //   d/dalpha [(1/alpha)*ln(1/(1+alpha*m))]
                //     = derivative of [ (1/alpha) * -ln(1+alpha*m ) ]
                //     = ...
                // For simplicity, we might do a numeric derivative. 
                // But let's do a partial symbolic approach.

                // We'll do partial symbolic => define a small helper function:
                auto logNBandDerivs = [&](double y_, double mu_, double alpha_) {
                    // This returns (value, d/dalpha, d2/dalpha^2)
                    // We'll do small steps or partial symbolic
                    // For advanced, you'd do full symbolic expansions.

                    // Numerically:
                    const double eps = 1e-8;
                    double valP = logNB(y_, mu_, alpha_ + eps);
                    double valM = logNB(y_, mu_, alpha_ - eps > 1e-15 ? alpha_ - eps : alpha_);
                    double val0 = logNB(y_, mu_, alpha_);

                    double firstDeriv = (valP - valM) / (2.0 * eps);
                    // second derivative => central difference of first derivative
                    double valPP = logNB(y_, mu_, alpha_ + 2.0 * eps);
                    double valMM = logNB(y_, mu_, alpha_ - 2.0 * eps > 1e-15 ? alpha_ - 2.0 * eps : alpha_);
                    double secondDeriv = (valPP - 2.0 * val0 + valMM) / (4.0 * eps * eps);

                    return std::make_tuple(val0, firstDeriv, secondDeriv);
                };

                auto tup = logNBandDerivs(y, m, alpha);
                double val0 = std::get<0>(tup);
                double der1 = std::get<1>(tup);
                double der2 = std::get<2>(tup);

                logLik += val0;
                d1 += der1;
                d2 += der2;
            }
        }

        // Subtract the Cox-Reid adjustment from logLik
        double crPen = coxReidAdjustment(mu, alpha);
        logLik -= crPen;

        // We also need derivative of the CR penalty:
        //   derivative of 0.5 * sum log(1 + alpha*mu_{i,r})
        // = 0.5 * sum [ mu_{i,r}/(1 + alpha*mu_{i,r}) ]
        double dCR = 0.0;
        double d2CR = 0.0;
        for (size_t i = 0; i < mu.size(); ++i) {
            for (size_t r = 0; r < mu[i].size(); ++r) {
                double val = 1.0 + alpha * mu[i][r];
                if (val > 0.0) {
                    double tmp = mu[i][r] / val; // derivative wrt alpha
                    dCR += tmp;
                    // second derivative wrt alpha => - mu[i][r]^2 / val^2
                    d2CR += -(mu[i][r] * mu[i][r]) / (val * val);
                }
            }
        }
        dCR *= 0.5;
        d2CR *= 0.5;

        // Incorporate CR derivative terms
        double fullD1 = d1 - dCR;
        double fullD2 = d2 - d2CR;

        // Newton step
        // alpha_new = alpha - (fullD1 / fullD2)
        // watch out for sign of fullD2
        if (std::fabs(fullD2) < 1e-15) {
            // fallback to small step or break
            // we can do a simple line search or secant approach
            double step = (fullD1 > 0) ? -1e-3 : 1e-3;
            double alphaOld = alpha;
            alpha += step;
            if (alpha <= 0) alpha = alphaOld * 0.5;
        } else {
            double delta = fullD1 / fullD2;
            double alphaOld = alpha;
            alpha = alpha - delta;
            if (alpha <= 1e-15) {
                alpha = alphaOld * 0.5; // fallback
            }
        }

        if (std::fabs(fullD1) < tol) {
            // converged
            break;
        }
    }

    return alpha;
}