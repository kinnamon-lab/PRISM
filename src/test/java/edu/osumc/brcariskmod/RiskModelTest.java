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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.util.ArithmeticUtils;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.MathArrays;
import org.junit.Assert;
import org.junit.Test;
import org.junit.experimental.runners.Enclosed;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameter;
import org.junit.runners.Parameterized.Parameters;

/**
 * JUnit test class for {@code RiskModel} class.
 * 
 * @author Daniel Kinnamon
 * @version 2014-10-07
 * @since 2014-10-05
 */
@RunWith(Enclosed.class)
public final class RiskModelTest {

  /**
   * Parameterized inner test class for {@code RiskModel} constructor argument
   * checks.
   */
  @RunWith(Parameterized.class)
  public static final class ConstructorArgTest {

    /**
     * Produces a collection of {@code Object} arrays containing parameters
     * defining each test case.
     */
    @Parameters(name = "ConstructorArgTest: Case {index}")
    public static final Collection<Object[]> cases() {
      // Declare ArrayList to hold test cases.
      final ArrayList<Object[]> casesList = new ArrayList<Object[]>();
      // Create test SNP object and arrays.
      final SNP testSNP =
        new SNP("rs1234", "Foo et al. Title. AJHG 2012; 21(3):1-5", "A", "T",
          SNP.AlleleOrientation.FORWARD, 0.14, -0.123);
      final SNP[] shortSNPArr = new SNP[10];
      Arrays.fill(shortSNPArr, testSNP);
      final SNP[] longSNPArr = new SNP[20];
      Arrays.fill(longSNPArr, testSNP);
      // Strictly increasing survivor function evaluation times.
      final double[] goodTimes = { 0.0D, 1.0D, 2.0D, 3.0D };
      // Non-decreasing survivor function evaluation times.
      final double[] nonDecTimes = { 0.0D, 1.0D, 1.0D, 2.0D };
      /*
       * Strictly increasing survivor function evaluation times with a negative
       * value.
       */
      final double[] negTimes = { -1.0D, 0.0D, 1.0D, 2.0D };
      // Non-increasing marginal survivor function values.
      final double[] goodMargSurv = { 1.0D, 0.75D, 0.75D, 0.5D };
      // Increasing marginal survivor function values.
      final double[] incMargSurv = { 1.0D, 0.75D, 0.8D, 0.5D };
      // Non-increasing marginal survivor function values, value above upper
      // bound.
      final double[] aboveUpperMargSurv = { 1.1D, 0.75D, 0.75D, 0.0D };
      // Non-increasing marginal survivor function values, value below lower
      // bound.
      final double[] belowLowerMargSurv = { 1.0D, 0.75D, 0.75D, -0.01D };
      // IllegalArgumentException: nonMonoTimes.
      casesList.add(new Object[] { shortSNPArr, nonDecTimes, goodMargSurv,
        true, IllegalArgumentException.class, "strictly increasing" });
      // IllegalArgumentException: negTimes.
      casesList.add(new Object[] { shortSNPArr, negTimes, goodMargSurv, true,
        IllegalArgumentException.class, "non-negative" });
      // IllegalArgumentException: incMargSurv.
      casesList.add(new Object[] { shortSNPArr, goodTimes, incMargSurv, true,
        IllegalArgumentException.class, "non-increasing" });
      // IllegalArgumentException: aboveUpperMargSurv.
      casesList.add(new Object[] { shortSNPArr, goodTimes, aboveUpperMargSurv,
        true, IllegalArgumentException.class, "[0,1]" });
      // IllegalArgumentException: belowLowerMargSurv.
      casesList.add(new Object[] { shortSNPArr, goodTimes, belowLowerMargSurv,
        true, IllegalArgumentException.class, "[0,1]" });
      /*
       * IllegalArgumentException: Too many SNPs for exact genotype
       * distribution.
       */
      casesList.add(new Object[] { longSNPArr, goodTimes, goodMargSurv, true,
        IllegalArgumentException.class, "exactGenoDist" });
      return casesList;
    }

    /*
     * Inject Object parameter array values into test class members for a
     * particular case.
     */
    /** {@code SNP} array of test SNPs. */
    @Parameter(0)
    public SNP[] testSNPs;
    /** {@code double} array of test survivor function evaluation times. */
    @Parameter(1)
    public double[] testTimes;
    /**
     * {@code double} array of test marginal survivor function values (i.e.,
     * {@code testMargSurv[i]} = S({@code testTimes[i]})).
     */
    @Parameter(2)
    public double[] testMargSurv;
    /**
     * {@code boolean} indicator of whether the exact genotype distribution
     * (TRUE) or a Monte Carlo approximation (FALSE) should be used for this
     * test.
     */
    @Parameter(3)
    public boolean testExactGenoDist;
    /** {@code Class<Exception>} expected exception type. */
    @Parameter(4)
    public Class<Exception> exceptType;
    /** {@code String} expected exception message substring. */
    @Parameter(5)
    public String exceptMessageSubstr;

    /**
     * Runs test on {@code ConstructorArgTest} object instantiated for a
     * particular test case.
     */
    @Test
    public final void test() {
      // Use try-catch block for exception checking.
      try {
        new RiskModel(testSNPs, testTimes, testMargSurv, testExactGenoDist);
        /*
         * Code in the try block below this comment should not be executed
         * because all test cases result in exceptions.
         */
        Assert.fail("Expected exception to be thrown.");
      } catch (Exception e) {
        // Check for expected exception.
        Assert.assertTrue("Expected " + exceptType.getSimpleName() + " to be " +
          "thrown.", e.getClass() == exceptType);
        Assert.assertTrue("Expected exception message to contain substring '" +
          exceptMessageSubstr + "'.",
          e.getMessage().contains(exceptMessageSubstr));
      }
    }

  }

  /**
   * Parameterized inner test class for {@code RiskModel} constructor baseline
   * survivor function solving.
   */
  @RunWith(Parameterized.class)
  public static final class SolveBaseSurvTest {

    /**
     * Produces a collection of {@code Object} arrays containing parameters
     * defining each test case.
     */
    @Parameters(name = "SolveBaseSurvTest: Case {index}")
    public static final Collection<Object[]> cases() {
      // Declare ArrayList to hold test cases.
      final ArrayList<Object[]> casesList = new ArrayList<Object[]>();
      /*
       * Declare new MersenneTwister RandomGenerator with known seed given by
       * the first 9 digits in the Euler-Mascheroni constant.
       */
      final MersenneTwister rngMT = new MersenneTwister(577215664);
      // Array of numbers of SNPs for various test cases.
      final int[] nSNPsArr =
        { 1, 1, 1, 1, 1, 5, 5, 5, 5, 5, 10, 10, 10, 10, 10, 15, 15, 15, 15, 15 };
      // Number of survivor function evaluation times to test.
      final int numTimes = 10;
      // Create test cases for each number of SNPs...
      for (int nSNPs : nSNPsArr) {
        // Declare Objects to hold simulated data for current test case.
        final double[] times = new double[numTimes];
        final double[] baseSurv = new double[numTimes];
        final SNP[] modelSNPs = new SNP[nSNPs];
        final BinomialDistribution[] genoDists =
          new BinomialDistribution[nSNPs];
        // For each locus...
        for (int locusIdx = 0; locusIdx < nSNPs; locusIdx++) {
          /*
           * Add new SNP with random allele 2 frequency uniformly distributed on
           * (0.05,0.95) and random per-allele 2 ln hazard ratio uniformly
           * distributed on (-0.5,0.5).
           */
          modelSNPs[locusIdx] =
            new SNP("rs" + locusIdx, "Foo et al. Title. AJHG 2012; 21(3):1-5",
              "A", "T", SNP.AlleleOrientation.FORWARD,
              0.05 + 0.9 * rngMT.nextDouble(), -0.5 + rngMT.nextDouble());
          // Instantiate BinomialDistribution object for current SNP.
          genoDists[locusIdx] =
            new BinomialDistribution(null, 2,
              modelSNPs[locusIdx].getAllele2Freq());
        }
        /*
         * Initialize survivor function evaluation times and baseline survivor
         * function values for test case.
         */
        for (int timeIdx = 0; timeIdx < times.length; timeIdx++) {
          // Add integer time, implicitly recast to double.
          times[timeIdx] = timeIdx;
          if (timeIdx < 2) {
            /*
             * Baseline survivor function values of 0 and 1 are included for all
             * test cases.
             */
            baseSurv[timeIdx] = timeIdx;
          } else {
            // Remaining values are random and uniformly distributed in (0,1).
            baseSurv[timeIdx] = rngMT.nextDouble();
          }
        }
        /*
         * Put baseSurv array into non-increasing order by sorting into
         * ascending numerical order and reversing in-place.
         */
        Arrays.sort(baseSurv);
        ArrayUtils.reverse(baseSurv);
        // Number of multivariant genotypes to be enumerated.
        final int nGenos = ArithmeticUtils.pow(3, modelSNPs.length);
        // Marginal survivor function values, all initialized to 0.
        final double[] margSurv = new double[baseSurv.length];
        Arrays.fill(margSurv, 0.0D);
        /*
         * Linear predictor and ln probability holders for current multivariant
         * genotype.
         */
        double genoEtaVal, genoLnProb;
        /*
         * Current quotient (for ternary conversion; see below) and locus
         * genotype (ternary digit) holders.
         */
        int quotient, locusGeno;
        // For each multivariant genotype...
        for (int genoIdx = 0; genoIdx < nGenos; genoIdx++) {
          // Initialize quotient holder to decimal multivariant genotype index.
          quotient = genoIdx;
          // Initialize linear predictor and ln probability holders to 0.
          genoEtaVal = 0.0D;
          genoLnProb = 0.0D;
          /*
           * For each locus in the current multivariant genotype (in reverse
           * order for ternary conversion; see below)...
           */
          for (int locusIdx = modelSNPs.length - 1; locusIdx >= 0; locusIdx--) {
            /*
             * There are 3^k unique k-digit ternary numbers, and each one maps
             * to a unique decimal number in the representable range [0, 3^k -
             * 1]. Multivariant genotypes can thus be enumerated by converting
             * each decimal multivariant genotype index in the representable
             * range to its ternary representation (e.g., 21011 for 5 SNPs).
             * Ternary digit l then represents the genotype for SNP l. Ternary
             * digits (SNP genotypes), from right to left, are obtained by
             * calculating the current quotient mod 3 and then dividing the
             * quotient by 3 for calculation of the next digit (SNP genotype) to
             * the left. This serves as an implementation check of the
             * MultidimensionalCounter object used to enumerate multivariant
             * genotypes in the RiskModel.solveBaseSurv method.
             */
            locusGeno = quotient % 3;
            quotient /= 3;
            genoEtaVal += locusGeno * modelSNPs[locusIdx].getAllele2LnHR();
            /*
             * Use genoDists[locusIdx] BinomialDistribution object here to test
             * implementation of modelSNPs[locus].getLnProbGeno used in
             * RiskModel.solveBaseSurv method.
             */
            genoLnProb += genoDists[locusIdx].logProbability(locusGeno);
          }
          for (int timeIdx = 0; timeIdx < baseSurv.length; timeIdx++) {
            /*
             * Add current multivariant genotype's contribution to the marginal
             * survivor function at each evaluation time.
             */
            margSurv[timeIdx] +=
              FastMath.exp(FastMath.log(baseSurv[timeIdx]) *
                FastMath.exp(genoEtaVal) + genoLnProb);
          }
        }
        /*
         * After loop over multivariant genotypes is finished, check for
         * margSurv elements > 1 due to minor roundoff error. If an element
         * exceeds 1 by more than 1e-4, throw a RuntimeException. Otherwise,
         * reset corresponding margSurv element to exactly 1 to correct for
         * minor roundoff error.
         */
        for (int timeIdx = 0; timeIdx < margSurv.length; timeIdx++) {
          if (margSurv[timeIdx] > 1.0001D) {
            throw new RuntimeException(
              "SolveBaseSurvTest: Calculated marginal survivor function "
                + "value for test case exceeds 1 by more than 1e-4.");
          } else if (margSurv[timeIdx] > 1.0D) {
            margSurv[timeIdx] = 1.0D;
          }
        }
        // Add test case using exact genotype distribution.
        casesList
          .add(new Object[] { baseSurv, modelSNPs, times, margSurv, true });
        // Add test case using Monte Carlo approximation.
        casesList.add(new Object[] { baseSurv, modelSNPs, times, margSurv,
          false });
      }
      return casesList;
    }

    /*
     * Inject Object parameter array values into test class members for a
     * particular case.
     */
    /** {@code double} array containing test baseline survivor function values. */
    @Parameter(0)
    public double[] testBaseSurv;
    /** {@code SNP} array of test SNPs. */
    @Parameter(1)
    public SNP[] testSNPs;
    /**
     * {@code double} array containing test baseline survivor function
     * evaluation times.
     */
    @Parameter(2)
    public double[] testTimes;
    /**
     * {@code double} array containing test marginal survivor function values
     * calculated from test baseline survivor function values and test SNPs.
     */
    @Parameter(3)
    public double[] testMargSurv;
    /**
     * {@code boolean} indicator of whether exact genotype distribution (TRUE)
     * or Monte Carlo approximation (FALSE) should be used for this test.
     */
    @Parameter(4)
    public boolean testExactGenoDist;

    /**
     * Runs test on {@code SolveBaseSurvTest} object instantiated for a
     * particular test case.
     */
    @Test
    public final void test() {
      /*
       * Run test of baseline survivor function solving in RiskModel
       * constructor.
       */
      final RiskModel testRiskModel =
        new RiskModel(testSNPs, testTimes, testMargSurv, testExactGenoDist);
      /*
       * Print original and solved baseline survivor function values as well as
       * differences.
       */
      System.out.println("<<<=== RiskModelTest.SolveBaseSurvTest: " +
        testSNPs.length + " SNPs, " +
        (testExactGenoDist ? "Exact" : "Monte Carlo") + " ===>>>\n");
      System.out.println("Original baseline survivor function values:");
      System.out.println(Arrays.toString(testBaseSurv));
      System.out.println("Solved baseline survivor function values:");
      System.out.println(Arrays.toString(testRiskModel.getBaseSurv()));
      System.out.println("Differences:");
      System.out.println(Arrays.toString(MathArrays.ebeSubtract(testBaseSurv,
        testRiskModel.getBaseSurv())) + "\n");
      /*
       * Check that original and solved baseline survivor function values are
       * all are within absolute tolerance of each other (1e-8 for direct
       * enumeration, 6.16e-4 for Monte Carlo approximation). The tolerance of
       * 6.16e-4 for the Monte Carlo approximation was chosen based upon the
       * Monte Carlo sampling error bound for the estimate of the marginal
       * survivor function at a particular evaluation time.
       */
      Assert.assertArrayEquals("Baseline survivor function values not " +
        "successfully recovered by RiskModel constructor with " +
        (testExactGenoDist ? "direct enumeration."
          : "Monte Carlo approximation."), testBaseSurv,
        testRiskModel.getBaseSurv(), testExactGenoDist ? 1e-8 : 6.16e-4);
    }

  }

}
