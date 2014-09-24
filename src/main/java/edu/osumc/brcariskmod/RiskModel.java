/*
 * This source file is part of the BRCA Risk Modifiers software package.
 * 
 * Copyright 2014 The Ohio State University Wexner Medical Center
 * 
 * Licensed under the Apache License, Version 2.0 (the "License"); you may not
 * use this file except in compliance with the License. You may obtain a copy of
 * the License at
 * 
 * http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
 * WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
 * License for the specific language governing permissions and limitations under
 * the License.
 */

package edu.osumc.brcariskmod;

import java.io.Serializable;
import java.util.Arrays;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.solvers.RiddersSolver;
import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.util.ArithmeticUtils;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.MathArrays;
import org.apache.commons.math3.util.MathArrays.OrderDirection;
import org.apache.commons.math3.util.MultidimensionalCounter;
import org.apache.commons.math3.util.Precision;

/**
 * Genetic risk prediction model for BRCA1/2 carriers.
 * <p>
 * This class implements a model that predicts breast and ovarian cancer risk
 * for a BRCA1/2 carrier as a function of her genotypes at several modifier
 * SNPs. The model is a variation on that proposed by: Antoniou AC, Beesley J,
 * McGuffog L, Sinilnikova OM, Healey S, Neuhausen SL, et al. Common breast
 * cancer susceptibility alleles and the risk of breast cancer for BRCA1 and
 * BRCA2 mutation carriers: implications for risk prediction. Cancer Res. 2010
 * Dec 1;70(23):9742–54. The primary difference is that the baseline survivor
 * function in this model is obtained directly at each age as a solution to a
 * well-behaved equation amenable to standard numerical root-finding techniques.
 * This class is serializable so that instances can be initialized with data and
 * then saved for later use.
 * 
 * @author Daniel Kinnamon
 * @version 2014-09-23
 * @since 2014-09-23
 */
public final class RiskModel implements Serializable {

  /**
   * The {@code serialVersionUID} should be incremented by 1 every time an
   * incompatible change (for serialization) is made to the class. If all
   * changes are compatible, then the {@code serialVersionUID} should not be
   * changed.
   */
  private static final long serialVersionUID = 1L;
  /** @serial {@code SNP} array containing information on all model SNPs. */
  private final SNP[] modelSNPs;
  /**
   * @serial {@code double} array of times at which the the survivor function is
   *         evaluated (abscissae).
   */
  private final double[] times;
  /**
   * @serial {@code double} array of marginal survivor function values, where
   *         {@code margSurv[i]} is the survivor function evaluated at
   *         {@code times[i]}.
   */
  private final double[] margSurv;
  /**
   * @serial {@code double} array of baseline survivor function values, where
   *         {@code baseSurv[i]} is the survivor function evaluated at
   *         {@code times[i]}.
   */
  private final double[] baseSurv;
  /**
   * @serial {@code double} array in which {@code genoEtaVals[i]} contains the
   *         value of the Cox model linear predictor (based on the ln hazard
   *         ratios in {@code modelSNPs}) for multivariant genotype i in the
   *         population or Monte Carlo sample.
   */
  private final double[] genoEtaVals;
  /**
   * @serial (@code double} array in which {@code genoLnProbs[i]} contains the
   *         ln probability of multivariant genotype i based on the allele
   *         frequencies in {@code modelSNPs} and the assumptions of HWE and
   *         linkage equilibrium. This array will be {@code null} for Monte
   *         Carlo samples.
   */
  private final double[] genoLnProbs;
  /**
   * @serial {@code boolean} indicator of whether calculations are based on the
   *         exact multivariant genotype distribution based on the allele
   *         frequencies in {@code modelSNPs} and the assumptions of HWE and
   *         linkage equilbrium (TRUE) or a Monte Carlo sample from it (FALSE).
   */
  private final boolean exactGenoDist;
  /**
   * @serial {@code int} constant number of SNPs above which a Monte Carlo
   *         sample from the multivariant genotype distribution must be used.
   */
  private static final int MAX_SNPS_EXACT = 15;
  /**
   * @serial {@code long} constant seed for random number generator (first 9
   *         digits of PI) to render Monte Carlo samples effectively
   *         deterministic.
   */
  private static final long RNG_SEED = 314159265L;
  /**
   * @serial {@code int} constant Monte Carlo sample size.
   *         <p>
   *         A sample size of 1e7 was chosen so that (1) the Monte Carlo
   *         approximation to the marginal survivor function for a given t
   *         should be within 6.16e-4 of its true value in 99.9% of samples
   *         (Hoeffding Inequality) and (2) the Monte Carlo estimate of the CDF
   *         of the linear predictor should be within 6.16e-4 of the true CDF
   *         everywhere in 99.9% of samples (Dvoretzky-Kiefer-Wolfowitz
   *         Inequality).
   */
  private static final int MONTE_CARLO_SAMP_SIZE = 10000000;
  /**
   * @serial {@code double} constant absolute accuracy within which two
   *         probabilities are considered to be equal--that is, p1 = p2 if
   *         abs(p1 - p2) <= this value.
   */
  private static final double PROB_CMP_EPSILON = 1e-10;
  /**
   * @serial {@code int} constant maximum number of function evaluations after
   *         which the {@link RiddersSolver} throws an exception if it has not
   *         converged.
   */
  private static final int SOLVER_MAX_EVAL = 100;

  /**
   * @serial {@code ObjFunException} is a private constant subclass of
   *         {@code RuntimeException} that is used to throw unchecked local
   *         exceptions from inside a {@link UnivariateFunction} that will not
   *         be caught by Apache Commons Math.
   */
  private static final class ObjFunException extends RuntimeException {

    private static final long serialVersionUID = 1L;
  }

  /**
   * Constructs a {@code RiskModel} object.
   * <p>
   * Automatically chooses whether to use the exact genotype distribution or a
   * Monte Carlo approximation to it based on whether the number of SNPs in
   * {@code modelSNPs} is <= MAX_SNPS_EXACT or not. The elements of the
   * {@code times} array must be strictly increasing with the index number and
   * non-negative. The elements of the {@code margSurv} array must be
   * non-increasing with the index number and in [0, 1].
   * 
   * @param modelSNPs {@link SNP} array containing all SNPs to be included in
   *          the model
   * @param times {@code double} array containing the times at which marginal
   *          survivor function values are available
   * @param margSurv {@code double} array in which {@code margSurv[i]} is the
   *          value of the marginal survivor function at {@code time[i]}
   */
  public RiskModel(final SNP[] modelSNPs, final double[] times,
    final double[] margSurv) {
    this(modelSNPs, times, margSurv, (modelSNPs.length <= MAX_SNPS_EXACT)
      ? true : false);
  }

  /**
   * Constructs a {@code RiskModel} object.
   * <p>
   * The elements of the {@code times} array must be strictly increasing with
   * the index number and non-negative. The elements of the {@code margSurv}
   * array must be non-increasing with the index number and in [0, 1].
   * 
   * @param modelSNPs {@link SNP} array containing all SNPs to be included in
   *          the model
   * @param times {@code double} array containing the times at which marginal
   *          survivor function values are available
   * @param margSurv {@code double} array in which {@code margSurv[i]} is the
   *          value of the marginal survivor function at {@code time[i]}
   * @param exactGenoDist {@code boolean} indicating whether the exact genotype
   *          distribution (TRUE) or a Monte Carlo sample from it (FALSE) should
   *          be used
   */
  public RiskModel(final SNP[] modelSNPs, final double[] times,
    final double[] margSurv, final boolean exactGenoDist) {
    // Copy constructor argument references to instance members.
    this.modelSNPs = modelSNPs;
    this.times = times;
    this.margSurv = margSurv;
    this.exactGenoDist = exactGenoDist;
    // Check that arguments resulted in a properly constructed instance.
    checkConstructorArgs();
    // Based on the value of exactGenoDist...
    if (exactGenoDist) {
      /*
       * Initialize genoEtaVals and genoLnProbs arrays of appropriate sizes with
       * NaN values.
       */
      final int nGenos = ArithmeticUtils.pow(3, modelSNPs.length);
      genoEtaVals = new double[nGenos];
      Arrays.fill(genoEtaVals, Double.NaN);
      genoLnProbs = new double[nGenos];
      Arrays.fill(genoLnProbs, Double.NaN);
      /*
       * Populate these arrays with linear predictor values and exact
       * multivariant genotype probabilities for all possible genotypes.
       */
      exactGenoDist();
    } else {
      /*
       * Initialize genoEtaVals array of appropriate size with NaN values and
       * unused genoLnProbs array with a null reference.
       */
      genoEtaVals = new double[MONTE_CARLO_SAMP_SIZE];
      Arrays.fill(genoEtaVals, Double.NaN);
      genoLnProbs = null;
      /*
       * Populate this array with linear predictor values for a random sample of
       * size MONTE_CARLO_SAMP_SIZE drawn from the multivariant genotype
       * distribution.
       */
      monteCarloGenoDist();
    }
    // Initialize baseSurv array of appropriate size with NaN values.
    baseSurv = new double[times.length];
    Arrays.fill(baseSurv, Double.NaN);
    /*
     * Populate this array with baseline survivor function values obtained by
     * numerically finding roots of an objective function based on margSurv,
     * genoEtaVals, and genoLnProbs (if applicable).
     */
    solveBaseSurv();
  }

  /** Checks arguments to RiskModel constructor. */
  private final void checkConstructorArgs() {
    // Check that times array elements are strictly increasing.
    if (!MathArrays.isMonotonic(times, OrderDirection.INCREASING, true)) {
      throw new IllegalArgumentException("RiskModel: times array "
        + "elements are not in strictly increasing order.");
    }
    // Check that times array elements are all non-negative.
    for (double time : times) {
      if (time < 0) {
        throw new IllegalArgumentException("RiskModel: times array "
          + "elements are not all non-negative.");
      }
    }
    // Check that margSurv array elements are non-increasing.
    if (!MathArrays.isMonotonic(margSurv, OrderDirection.DECREASING, false)) {
      throw new IllegalArgumentException("RiskModel: margSurv array "
        + "elements are not in non-increasing order.");
    }
    // Check that margSurv array elements are in [0, 1].
    for (double margSurvT : margSurv) {
      if (margSurvT < 0.0D || margSurvT > 1.0D) {
        throw new IllegalArgumentException("RiskModel: margSurv array "
          + "elements are not all in [0,1].");
      }
    }
    /*
     * Check that exactGenoDist argument is compatible with length of modelSNPs
     * array.
     */
    if ((modelSNPs.length > MAX_SNPS_EXACT) && exactGenoDist) {
      throw new IllegalArgumentException("RiskModel: Cannot calculate " +
        "exact genotype distribution with more than" + MAX_SNPS_EXACT +
        "SNPs. Please change exactGenoDist argument to false.");
    }
  }

  /**
   * Calculates linear predictor values and exact probabilities for all possible
   * multivariant genotypes.
   */
  private final void exactGenoDist() {
    /*
     * Initialize array of 3s, one for each SNP, to initialize
     * MultidimensionalCounter object. This object translates each univariate
     * index in [1, 3^modelSNPs.length] to a unique multivariant genotype
     * represented as an int array with modelSNPs.length slots each taking
     * values in [0,2].
     */
    int[] dimSizes = new int[modelSNPs.length];
    Arrays.fill(dimSizes, 3);
    MultidimensionalCounter genoCounter = new MultidimensionalCounter(dimSizes);
    /*
     * Check that genoCounter represents correct number of distinct multivariant
     * genotypes.
     */
    if (genoCounter.getSize() != genoEtaVals.length) {
      throw new RuntimeException("MultidimensionalCounter size does not "
        + "match calculated number of genotypes.");
    }
    // Initialize variable to check cumulative probability over all genotypes.
    double cumProb = 0.0D;
    // For each multivariant genotype...
    for (int genoIdx = 0; genoIdx < genoEtaVals.length; genoIdx++) {
      /*
       * Initialize linear predictor and probability for multivariant genotype
       * to zero.
       */
      genoEtaVals[genoIdx] = 0.0D;
      genoLnProbs[genoIdx] = 0.0D;
      // For each locus...
      for (int locusIdx = 0; locusIdx < modelSNPs.length; locusIdx++) {
        /*
         * Get genotype for the current locusIdx (as number of 2 alleles) from
         * the appropriate position of the multivariant genotype int array for
         * the current genoIdx and add X_j*Beta to the linear predictor.
         */
        genoEtaVals[genoIdx] +=
          genoCounter.getCounts(genoIdx)[locusIdx] *
            modelSNPs[locusIdx].getAllele2LnHR();
        /*
         * Get genotype for the current locusIdx (as number of 2 alleles) from
         * the appropriate position of the multivariant genotype int array for
         * the current genoIdx and add the ln probability of that genotype under
         * HWE returned from the SNP object for locusIdx to genoLnProbs. This
         * method of calculating the multivariant genotype probability assumes
         * that SNPs are in linkage equilibrium.
         */
        genoLnProbs[genoIdx] +=
          modelSNPs[locusIdx]
            .getLnProbGeno(genoCounter.getCounts(genoIdx)[locusIdx]);
      }
      /*
       * Exponentiate sum of ln genotype probabilities over loci to get the
       * multivariant genotype probability (i.e., the product of the
       * single-variant genotype probabilities across loci assuming linkage
       * equilibrium) and add to cumulative total.
       */
      cumProb += FastMath.exp(genoLnProbs[genoIdx]);
    }
    // Check that cumulative probability over all genotypes is numerically 1.
    if (!Precision.equals(cumProb, 1.0D, PROB_CMP_EPSILON)) {
      throw new RuntimeException("Cumulative probability over all possible "
        + "multivariant genotypes is not numerically 1.");
    }
  }

  /**
   * Calculates linear predictor values for a random Monte Carlo sample of
   * multivariant genotypes drawn from the appropriate distribution.
   */
  private final void monteCarloGenoDist() {
    /*
     * Initialize single Mersenne-Twister RNG reference to be used for all
     * random number generation. This ensures that sequential (i.e.,
     * independent) values from the underlying random number stream are used.
     */
    final MersenneTwister rngMT = new MersenneTwister(RNG_SEED);
    // For each multivariant genotype...
    for (int genoIdx = 0; genoIdx < MONTE_CARLO_SAMP_SIZE; genoIdx++) {
      // Initialize linear predictor for multivariant genotype to zero.
      genoEtaVals[genoIdx] = 0.0D;
      // For each locus...
      for (int locusIdx = 0; locusIdx < modelSNPs.length; locusIdx++) {
        /*
         * Randomly sample single-variant genotype (as number of 2 alleles) from
         * SNP object for current locusIdx using rngMT as the underlying random
         * number source and add X_j*Beta to the linear predictor. This method
         * of simulating multivariant genotypes by independently drawing
         * single-variant genotypes for each locusIdx assumes linkage
         * equilibrium.
         */
        genoEtaVals[genoIdx] +=
          modelSNPs[locusIdx].getRandGeno(rngMT) *
            modelSNPs[locusIdx].getAllele2LnHR();
      }
    }
  }

  /**
   * Solves for the baseline survivor function at each
   * {@code (times[i], margSurv[i])} point.
   * <p>
   * Obtains the baseline survivor function {@code (times[i], baseSurv[i])}
   * points that would lead to the observed marginal survivor function points
   * {@code (times[i], margSurv[i])}, assuming the given allele 2 frequencies,
   * HWE, linkage equilibrium, and a Cox model with the given allele 2 ln hazard
   * ratios.
   */
  private final void solveBaseSurv() {
    /*
     * Initialize new RiddersSolver object with absolute accuracy tolerance
     * given by PROB_CMP_EPSILON (default relative and function value tolerances
     * are adequately small). Ridders' method is a regula falsi (bracketing)
     * root finding algorithm that is guaranteed to converge and has good speed
     * properties. See: (1) Ridders C. A new algorithm for computing a single
     * root of a real continuous function. IEEE Transactions on Circuits and
     * Systems. 1979 Nov;26(11):979–80. (2) Press et al. Numerical Recipes, 3rd
     * Ed. New York: Cambridge University Press, 2007. Section 9.2.1, pp. 452-4.
     */
    final RiddersSolver solver = new RiddersSolver(PROB_CMP_EPSILON);
    // For each time...
    for (int timeIdx = 0; timeIdx < times.length; timeIdx++) {
      /*
       * Get marginal survivor function for the current timeIdx. Must use final
       * double declared within loop in order to capture this local variable in
       * the anonymous class definition for objFun.
       */
      final double margSurvT = margSurv[timeIdx];
      /*
       * Define anonymous class implementing the objective function for the
       * current time, which is the difference between the estimated marginal
       * survivor function calculated from the current estimate of the baseline
       * survivor function and the observed marginal survivor function.
       */
      UnivariateFunction objFun = new UnivariateFunction() {

        @Override
        public double value(final double s0T) {
          // Throw ObjFunException if s0T is not in the appropriate range.
          if (s0T < 0.0D || s0T > 1.0D) {
            throw new ObjFunException();
          }
          // Try to calculate objective function value at s0T.
          try {
            if (Precision.equals(s0T, 0.0D)) {
              /*
               * If s0T is numerically zero (within 1 ulp), then objective
               * function value is -margSurvT.
               */
              return -margSurvT;
            } else if (Precision.equals(s0T, 1.0D)) {
              /*
               * Otherwise, if s0T is numerically 1 (within 1 ulp), then
               * objective function value is 1 - margSurvT.
               */
              return 1.0D - margSurvT;
            } else {
              /*
               * Otherwise, must take (expectation/simple average) of
               * s0T^exp(genoEtaVal) over all (possible/sampled) genotypes for
               * (exact/Monte Carlo) calculation of objective function value.
               */
              double objFunVal = 0.0D;
              for (int genoIdx = 0; genoIdx < genoEtaVals.length; genoIdx++) {
                // @formatter:off
                /* EXACT SUM ELEMENTS:
                 * exp(genoLnProb)*s0T^exp(genoEtaVal)
                 * = exp(genoLnProb)*exp(ln(s0T))^exp(genoEtaVal)
                 * = exp(genoLnProb+ln(s0T)*exp(genoEtaVal))
                 * 
                 * MONTE CARLO SUM ELEMENTS:
                 * s0T^exp(genoEtaVal)
                 * = exp(ln(s0T))^exp(genoEtaVal)
                 * = exp(0.0D + ln(s0T)*exp(genoEtaVal))
                 */
                // @formatter:on
                objFunVal +=
                  FastMath.exp((exactGenoDist ? genoLnProbs[genoIdx] : 0.0D) +
                    FastMath.log(s0T) * FastMath.exp(genoEtaVals[genoIdx]));
              }
              /*
               * If (exact/Monte Carlo) calculation, obtain final objective
               * function value by subtracting margSurvT from the (sum/sum
               * divided by the Monte Carlo sample size), which is the
               * (expectation/simple average) of s0T^exp(genoEtaVal) over all
               * (possible/sampled) genotypes.
               */
              objFunVal =
                exactGenoDist ? (objFunVal - margSurvT) : (objFunVal /
                  MONTE_CARLO_SAMP_SIZE - margSurvT);
              return objFunVal;
            }
          } catch (Exception e) {
            /*
             * A locally defined exception is thrown if any of the underlying
             * calculations throw an exception. No value is returned.
             */
            throw new ObjFunException();
          }
        }
      };
      if (Precision.equals(margSurvT, 1.0D)) {
        /*
         * The observed marginal survivor function is exactly 1 if and only if
         * the baseline survivor function is exactly 1.
         */
        baseSurv[timeIdx] = 1.0D;
      } else if (Precision.equals(margSurvT, 0.0D)) {
        /*
         * The observed marginal survivor function is exactly 0 if and only if
         * the baseline survivor function is exactly 0.
         */
        baseSurv[timeIdx] = 0.0D;
      } else {
        /*
         * Otherwise, numerically solve for baseline survivor function at time t
         * in the search interval [0.0D, 1.0D]. Note that bracketing guarantees
         * that the root will be within the search interval, so range checks of
         * the solution are unnecessary.
         */
        baseSurv[timeIdx] = solver.solve(SOLVER_MAX_EVAL, objFun, 0.0D, 1.0D);
        /*
         * If this is not the first time, is the solved baseline survivor
         * function within PROB_CMP_EPSILON of the value at the previous time?
         * If so, it should exactly equal this value.
         */
        if ((timeIdx > 0) &&
          Precision.equals(baseSurv[timeIdx], baseSurv[timeIdx - 1],
            PROB_CMP_EPSILON)) {
          baseSurv[timeIdx] = baseSurv[timeIdx - 1];
        }
      }
    }
    // Check that solved baseline survivor function is non-increasing over time.
    if (!MathArrays.isMonotonic(baseSurv, OrderDirection.DECREASING, false)) {
      throw new RuntimeException("RiskModel: Solved baseline survivor "
        + "function is increasing as a function of time.");
    }
  }

  /** Returns {@code double} array of baseline survivor function values. */
  final double[] getBaseSurv() {
    return baseSurv;
  }

}