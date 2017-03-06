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

package edu.osumc.PRISM;

import java.io.PrintWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import org.apache.commons.lang3.StringUtils;
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
 * @version 2014-10-31
 * @since 2014-09-23
 */
public final class RiskModel implements Serializable {

  /**
   * The {@code serialVersionUID} should be automatically regenerated within the
   * IDE whenever the fields of the class change. It should not be changed if
   * only the methods change.
   */
  private static final long serialVersionUID = 837026756986969984L;
  /** @serial {@code String} risk model name. */
  private final String modelName;
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
   *         linkage equilbrium {@code (true)} or a Monte Carlo sample from it
   *         {@code (false)}.
   */
  private final boolean useExactGenoDist;
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
   * @serial {@code ObjFunException} is a private static nested class extending
   *         {@code RuntimeException} that is used to throw local exceptions
   *         from inside a {@link UnivariateFunction} so that they will not be
   *         caught by Apache Commons Math.
   */
  private static final class ObjFunException extends RuntimeException {

    private static final long serialVersionUID = 4577459979033121052L;

    /**
     * Constructs an ObjFunException object that packages the underlying
     * exception.
     */
    private ObjFunException(final Throwable cause) {
      super(cause);
    }
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
   * @param modelName {@code String} risk model name
   * @param modelSNPs {@link SNP} array containing all SNPs to be included in
   *          the model
   * @param times {@code double} array containing the times at which marginal
   *          survivor function values are available
   * @param margSurv {@code double} array in which {@code margSurv[i]} is the
   *          value of the marginal survivor function at {@code time[i]}
   */
  public RiskModel(final String modelName, final SNP[] modelSNPs,
    final double[] times, final double[] margSurv) {
    this(modelName, modelSNPs, times, margSurv,
      (modelSNPs.length <= MAX_SNPS_EXACT) ? true : false);
  }

  /**
   * Constructs a {@code RiskModel} object.
   * <p>
   * The elements of the {@code times} array must be strictly increasing with
   * the index number and non-negative. The elements of the {@code margSurv}
   * array must be non-increasing with the index number and in [0, 1].
   * 
   * @param modelName {@code String} risk model name
   * @param modelSNPs {@link SNP} array containing all SNPs to be included in
   *          the model
   * @param times {@code double} array containing the times at which marginal
   *          survivor function values are available
   * @param margSurv {@code double} array in which {@code margSurv[i]} is the
   *          value of the marginal survivor function at {@code time[i]}
   * @param useExactGenoDist {@code boolean} indicating whether the exact
   *          genotype distribution {@code (true)} or a Monte Carlo sample from
   *          it {@code (false)} should be used.
   */
  public RiskModel(final String modelName, final SNP[] modelSNPs,
    final double[] times, final double[] margSurv,
    final boolean useExactGenoDist) {
    // Copy constructor argument references to instance members.
    this.modelName = modelName;
    this.modelSNPs = modelSNPs;
    this.times = times;
    this.margSurv = margSurv;
    this.useExactGenoDist = useExactGenoDist;
    // Check that arguments resulted in a properly constructed instance.
    checkConstructorArgs();
    // Based on the value of useExactGenoDist...
    if (useExactGenoDist) {
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
     * Check that useExactGenoDist argument is compatible with length of
     * modelSNPs array.
     */
    if ((modelSNPs.length > MAX_SNPS_EXACT) && useExactGenoDist) {
      throw new IllegalArgumentException("RiskModel: Cannot calculate " +
        "exact genotype distribution with more than" + MAX_SNPS_EXACT +
        "SNPs. Please change useExactGenoDist argument to false.");
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
          // Try to calculate objective function value at s0T.
          try {
            if (s0T < 0.0D || s0T > 1.0D) {
              /*
               * If s0T is not in the appropriate range, throw
               * IllegalArgumentException.
               */
              throw new IllegalArgumentException(
                "Baseline survivor function guess not in [0,1].");
            } else if (Precision.equals(s0T, 0.0D)) {
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
                  FastMath
                    .exp((useExactGenoDist ? genoLnProbs[genoIdx] : 0.0D) +
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
                useExactGenoDist ? (objFunVal - margSurvT) : (objFunVal /
                  MONTE_CARLO_SAMP_SIZE - margSurvT);
              return objFunVal;
            }
          } catch (Exception e) {
            /*
             * If any of the underlying calculations throw an exception, this is
             * caught, packaged in a local unchecked ObjFunException, and
             * rethrown. No value is returned.
             */
            throw new ObjFunException(e);
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
        try {
          baseSurv[timeIdx] = solver.solve(SOLVER_MAX_EVAL, objFun, 0.0D, 1.0D);
        } catch (ObjFunException e) {
          /*
           * If an ObjFunException is thrown, repackage the underlying cause in
           * a RuntimeException recognizable outside of this class and throw it
           * up the call stack.
           */
          throw new RuntimeException(e.getCause());
        }
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

  /**
   * Returns an {@link Individual.RiskPrediction} object based on this
   * {@code RiskModel} object and the information provided in the
   * {@code inputGenos} argument.
   * <p>
   * {@code inputGenos} can contain genotypes for all typed SNPs in an
   * individual, including those that are not used in the risk model. These are
   * simply ignored. SNPs used in the risk model but not included in
   * {@code inputGenos} are treated as missing, and the individual's prognostic
   * index in the {@code Individual.RiskPrediction} object is given by the
   * expected Cox model linear predictor calculated by integrating over all
   * missing genotypes under the assumptions of HWE and linkage equilibrium.
   * Absolute cumulative risk predictions for the individual are based on this
   * expected prognostic index. A {@code LinkedHashMap} of SNP genotypes
   * actually used in the risk model can be extracted from the resulting
   * {@code Individual.RiskPrediction} object. SNPs are added to this
   * {@code LinkedHashMap} in the same order that they are stored in the
   * {@code RiskModel} object.
   * 
   * @param inputGenos {@link Individual.Genotypes} object containing an
   *          individual's genotype information
   */
  public final Individual.RiskPrediction getRiskPrediction(
    final Individual.Genotypes inputGenos) {
    // Initialize Individual.RiskPrediction object.
    Individual.RiskPrediction riskPrediction =
      new Individual.RiskPrediction(inputGenos.getIndivID(), modelName);
    // Calculate (expected) Cox model linear predictor.
    double etaVal = 0.0D;
    for (SNP curSNP : modelSNPs) {
      final String rsID = curSNP.getRsID();
      final String allele1 = inputGenos.getAllele1(rsID);
      final String allele2 = inputGenos.getAllele2(rsID);
      final SNP.AlleleOrientation orientRs =
        (inputGenos.getOrientRs(rsID) != null) ? inputGenos.getOrientRs(rsID)
          : curSNP.getOrientRs();
      riskPrediction.addUsedGenotype(rsID, allele1, allele2);
      etaVal += curSNP.genoScore(allele1, allele2, orientRs);
    }
    /*
     * Calculate linear predictor percentile in population, assuming HWE and
     * linkage equilibrium between SNPs.
     */
    double etaPctl = 0.0D;
    for (int genoIdx = 0; genoIdx < genoEtaVals.length; genoIdx++) {
      if (genoEtaVals[genoIdx] <= etaVal) {
        etaPctl += useExactGenoDist ? FastMath.exp(genoLnProbs[genoIdx]) : 1.0D;
      }
    }
    if (!useExactGenoDist) {
      etaPctl /= genoEtaVals.length;
    }
    riskPrediction.setPrognosticIndex(etaVal, etaPctl);
    /*
     * Calculate predicted cumulative risk function (1 minus predicted survivor
     * function under assumed continuous time model).
     */
    for (int timeIdx = 0; timeIdx < times.length; timeIdx++) {
      final double predCumRiskT;
      if (Precision.equals(baseSurv[timeIdx], 1.0D)) {
        /*
         * If the baseline survivor function is exactly 1, then the predicted
         * cumulative risk function is exactly 1 - 1 = 0.
         */
        predCumRiskT = 0.0D;
      } else if (Precision.equals(baseSurv[timeIdx], 0.0D)) {
        /*
         * If the baseline survivor function is exactly 0, then the predicted
         * cumulative risk function is exactly 1 - 0 = 1.
         */
        predCumRiskT = 1.0D;
      } else {
        predCumRiskT =
          1.0D - FastMath.exp(FastMath.log(baseSurv[timeIdx]) *
            FastMath.exp(etaVal));
      }
      riskPrediction.addPredCumRiskT(times[timeIdx], predCumRiskT);
    }
    return riskPrediction;
  }

  /**
   * Prints {@code RiskModel} object information using a {@code PrintWriter}.
   * 
   * @param pw {@code PrintWriter} for output
   */
  public final void print(final PrintWriter pw) {
    final ArrayList<String> sourcePubList = new ArrayList<String>();
    pw.println(StringUtils.repeat("=", modelName.length()));
    pw.println(modelName);
    pw.println(StringUtils.repeat("=", modelName.length()));
    pw.println();
    pw.println("SUMMARY");
    pw.println("-------");
    pw.println("Number of SNPs Included: " + modelSNPs.length);
    pw.println("Genotype Distribution: " +
      (useExactGenoDist ? "Direct Enumeration" : "Monte Carlo Approximation"));
    if (!useExactGenoDist) {
      pw.println("Monte Carlo Sample Size: " + MONTE_CARLO_SAMP_SIZE);
    }
    pw.println("Number of Survivor Function Evaluation Times: " + times.length);
    pw.println();
    pw.println("MODEL SNP DETAILS");
    final String snpTableHeader =
      String.format(Locale.US,
        "%1$-14s | %2$-8s | %3$-8s | %4$-7s | %5$-7s | %6$-7s | %7$-3s",
        "RS #", "A1", "A2", "ORIENT", "A2 FREQ", "A2 HR", "REF");
    final String snpTableHrule =
      snpTableHeader.replaceAll("[^|]", "-").replace("|", "+");
    pw.println(snpTableHrule);
    pw.println(snpTableHeader);
    pw.println(snpTableHrule);
    for (SNP curSNP : modelSNPs) {
      /*
       * If source publication is not in the list, add it to the end. The
       * reference number will be the index of the source publication in the
       * list + 1.
       */
      if (!sourcePubList.contains(curSNP.getSourcePub())) {
        sourcePubList.add(curSNP.getSourcePub());
      }
      pw.format(Locale.US,
        "%1$-14s | %2$-8s | %3$-8s | %4$-7s | %5$-7.3f | %6$-7.3f | %7$-3d%n",
        curSNP.getRsID(), curSNP.getAllele1(), curSNP.getAllele2(), ((curSNP
          .getOrientRs() == SNP.AlleleOrientation.FORWARD) ? "Forward"
          : "Reverse"), curSNP.getAllele2Freq(), FastMath.exp(curSNP
          .getAllele2LnHR()), sourcePubList.indexOf(curSNP.getSourcePub()) + 1);
    }
    pw.println(snpTableHrule);
    /*
     * Loop over reference list, outputting reference number and publication
     * information.
     */
    for (int refNum = 1; refNum <= sourcePubList.size(); refNum++) {
      pw.println(refNum + ") " + sourcePubList.get(refNum - 1));
    }
    pw.println();
    pw.println("SURVIVOR FUNCTIONS");
    final String lifeTableHeader =
      String
        .format(Locale.US, "%1$-7s | %2$-7s | %3$-7s", "t", "S(t)", "So(t)");
    final String lifeTableHrule =
      lifeTableHeader.replaceAll("[^|]", "-").replace("|", "+");
    pw.println(lifeTableHrule);
    pw.println(lifeTableHeader);
    pw.println(lifeTableHrule);
    for (int timeIdx = 0; timeIdx < times.length; timeIdx++) {
      pw.format(Locale.US, "%1$-7.3f | %2$-7.3f | %3$-7.3f%n", times[timeIdx],
        margSurv[timeIdx], baseSurv[timeIdx]);
    }
    pw.println(lifeTableHrule);
    // Flush PrintWriter to ensure that output is produced.
    pw.flush();
  }

  /** Returns {@code String} model name. */
  final String getModelName() {
    return modelName;
  }

  /** Returns {@code double} array of baseline survivor function values. */
  public final double[] getBaseSurv() {
    return baseSurv;
  }

}