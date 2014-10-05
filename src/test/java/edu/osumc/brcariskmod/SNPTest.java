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
import java.util.Collection;

import org.junit.Assert;
import org.junit.Test;
import org.junit.experimental.runners.Enclosed;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameter;
import org.junit.runners.Parameterized.Parameters;

/**
 * JUnit test class for {@code SNP} class.
 * 
 * @author Daniel Kinnamon
 * @version 2014-10-05
 * @since 2014-10-03
 */
@RunWith(Enclosed.class)
public final class SNPTest {

  /** Parameterized inner test class for {@code SNP.genoScore} method. */
  @RunWith(Parameterized.class)
  public static final class GenoScoreTest {

    /**
     * Produces a collection of {@code Object} arrays containing parameters
     * defining each test case.
     */
    @Parameters(name = "GenoScoreTest: Case {index}")
    public static final Collection<Object[]> cases() {
      // Declare ArrayList to hold test cases.
      final ArrayList<Object[]> casesList = new ArrayList<Object[]>();
      // Create two test case SNP objects.
      final SNP testSNP1 =
        new SNP("rs1", "Foo et al. Title. AJHG 2012; 21(3):1-5", "A", "g",
          SNP.AlleleOrientation.FORWARD, 0.2D, 0.5D);
      final SNP testSNP2 =
        new SNP("rs2", "Foo et al. Title. AJHG 2012; 21(3):1-5", "AtTacGcG",
          "-", SNP.AlleleOrientation.REVERSE, 0.5D, 0.25D);
      // IllegalArgumentException: Bad input allele.
      casesList.add(new Object[] { testSNP1, "A", "N",
        SNP.AlleleOrientation.FORWARD, Double.NaN,
        IllegalArgumentException.class, "invalid input allele" });
      // IllegalArgumentException: Bad input allele.
      casesList.add(new Object[] { testSNP1, "/", "G",
        SNP.AlleleOrientation.FORWARD, Double.NaN,
        IllegalArgumentException.class, "invalid input allele" });
      // IllegalArgumentException: One allele is missing, one is observed.
      casesList.add(new Object[] { testSNP1, "A", "0",
        SNP.AlleleOrientation.FORWARD, Double.NaN,
        IllegalArgumentException.class, "Neither or both" });
      // IllegalArgumentException: One allele is missing, one is observed.
      casesList.add(new Object[] { testSNP1, "0", "G",
        SNP.AlleleOrientation.FORWARD, Double.NaN,
        IllegalArgumentException.class, "Neither or both" });
      // Should return expected score for testSNP1 (0.2).
      casesList.add(new Object[] { testSNP1, "0", "0",
        SNP.AlleleOrientation.FORWARD, 0.2D, null, null });
      // Should return testSNP1 score for 0 copies of allele2.
      casesList.add(new Object[] { testSNP1, "a", "A",
        SNP.AlleleOrientation.FORWARD, 0.0D, null, null });
      // Should return testSNP1 score for 2 copies of allele2.
      casesList.add(new Object[] { testSNP1, "C", "c",
        SNP.AlleleOrientation.REVERSE, 1.0D, null, null });
      // Should return testSNP1 score for 1 copy of allele2.
      casesList.add(new Object[] { testSNP1, "t", "c",
        SNP.AlleleOrientation.REVERSE, 0.5D, null, null });
      /*
       * IllegalArgumentException: T and C are not possible alleles for testSNP1
       * if orientation is the same as in SNP object.
       */
      casesList.add(new Object[] { testSNP1, "t", "C",
        SNP.AlleleOrientation.FORWARD, Double.NaN,
        IllegalArgumentException.class, "possible population" });
      // Should return expected score for testSNP2 (0.25).
      casesList.add(new Object[] { testSNP2, "0", "0",
        SNP.AlleleOrientation.FORWARD, 0.25D, null, null });
      // Should return testSNP2 score for 2 copies of allele2.
      casesList.add(new Object[] { testSNP2, "-", "-",
        SNP.AlleleOrientation.FORWARD, 0.5D, null, null });
      // Should return testSNP2 score for 2 copies of allele2.
      casesList.add(new Object[] { testSNP2, "-", "-",
        SNP.AlleleOrientation.REVERSE, 0.5D, null, null });
      // Should return testSNP2 score for 1 copy of allele2.
      casesList.add(new Object[] { testSNP2, "attacgcg", "-",
        SNP.AlleleOrientation.REVERSE, 0.25D, null, null });
      // Should return testSNP2 score for 1 copy of allele2.
      casesList.add(new Object[] { testSNP2, "-", "TAATGCGC",
        SNP.AlleleOrientation.FORWARD, 0.25D, null, null });
      // Should return testSNP2 score for 0 copies of allele2.
      casesList.add(new Object[] { testSNP2, "taatgcgc", "TAATGCGC",
        SNP.AlleleOrientation.FORWARD, 0.0D, null, null });
      /*
       * IllegalArgumentException: TaaTGcGC is not a possible allele for
       * testSNP2 if orientation is the same as in SNP object.
       */
      casesList.add(new Object[] { testSNP2, "-", "TaaTGcGC",
        SNP.AlleleOrientation.REVERSE, Double.NaN,
        IllegalArgumentException.class, "possible population" });
      return casesList;
    }

    /*
     * Inject Object parameter array values into test class members for a
     * particular case.
     */
    /** {@code SNP} test SNP. */
    @Parameter(0)
    public SNP testSNP;
    /** {@code String} test allele 1. */
    @Parameter(1)
    public String testAllele1;
    /** {@code String} test allele 2. */
    @Parameter(2)
    public String testAllele2;
    /**
     * {@code SNP.AlleleOrientation} test allele orientation, relative to RefSNP
     * alleles.
     */
    @Parameter(3)
    public SNP.AlleleOrientation testOrientRs;
    /** {@code double} correct score for genotype given by test alleles 1/2. */
    @Parameter(4)
    public double correctScore;
    /**
     * {@code Class<Exception>} expected exception type (should be {@code null}
     * if no exception expected).
     */
    @Parameter(5)
    public Class<Exception> exceptType;
    /**
     * {@code String} expected exception message substring (should be
     * {@code null} if no exception expected).
     */
    @Parameter(6)
    public String exceptMessageSubstr;

    /**
     * Runs test on {@code GenoScoreTest} object instantiated for a particular
     * test case.
     */
    @Test
    public final void test() {
      // Use try-catch block for exception checking.
      try {
        double actualScore =
          testSNP.genoScore(testAllele1, testAllele2, testOrientRs);
        /*
         * If error-checking test case, code below this line should not be
         * executed.
         */
        if (exceptType != null || exceptMessageSubstr != null) {
          Assert.fail("Expected exception to be thrown.");
        }
        /*
         * Get actual score for test case and compare to expected score (with
         * absolute error of 1e-10).
         */
        Assert.assertEquals("Returned genotype score not within 1e-10 of "
          + "correct value.", correctScore, actualScore, 1e-10);
      } catch (Exception e) {
        // If not error-checking test case, no exception should occur.
        if (exceptType == null && exceptMessageSubstr == null) {
          Assert.fail("Unexpected exception was thrown.");
        }
        // Otherwise, check for expected exception.
        Assert.assertTrue("Expected " + exceptType.getSimpleName() + " to be " +
          "thrown.", e.getClass() == exceptType);
        Assert.assertTrue("Expected exception message to contain substring '" +
          exceptMessageSubstr + "'.",
          e.getMessage().contains(exceptMessageSubstr));
      }

    }

  }

}
