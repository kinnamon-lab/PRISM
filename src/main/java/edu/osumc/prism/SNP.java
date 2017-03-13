/*
 * This source file is part of the PRISM software package.
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

package edu.osumc.prism;

import java.io.Serializable;
import java.util.Locale;

import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.util.FastMath;

/**
 * Single nucleotide polymorphism (SNP) characteristics.
 * <p>
 * This class stores SNP information obtained from a source publication,
 * including the dbSNP refSNP identifier, the source publication reference,
 * alleles 1 and 2 from the source publication, the orientation of these alleles
 * relative to the dbSNP refSNP alleles, the allele 2 frequency given in the
 * source publication, and the Cox proportional hazards model ln hazard ratio
 * for each additional allele 2 obtained from the source publication. This class
 * is serializable.
 * 
 * @author Daniel Kinnamon
 * @author Chuhan Zhang
 * @version 2014-10-27
 * @since 2014-09-01
 */
public final class SNP implements Serializable {

  /**
   * The {@code serialVersionUID} should be automatically regenerated within the
   * IDE whenever the fields of the class change. It should not be changed if
   * only the methods change.
   */
  private static final long serialVersionUID = -4177085957451410411L;

  /** Allele orientation values for a SNP */
  public static enum AlleleOrientation {
    FORWARD, REVERSE
  }

  /** @serial {@code String} dbSNP refSNP identifier (i.e., rs number). */
  private final String rsID;
  /** @serial {@code String} source publication reference. */
  private final String sourcePub;
  /** @serial {@code String} uppercase allele 1 from source publication. */
  private final String allele1;
  /** @serial {@code String} uppercase allele 2 from source publication. */
  private final String allele2;
  /**
   * @serial {@code AlleleOrientation} orientation of source publication alleles
   *         relative to refSNP alleles.
   */
  private final AlleleOrientation orientRs;
  /** @serial {@code double} allele 2 frequency from source publication. */
  private final double allele2Freq;
  /**
   * @serial {@code double} Cox proportional hazards model ln hazard ratio for
   *         each additional allele 2 obtained from the source publication.
   */
  private final double allele2LnHR;

  /**
   * Constructs a SNP object.
   * 
   * @param rsID {@code String} dbSNP refSNP identifier (i.e., rs number)
   * @param sourcePub {@code String} source publication reference
   * @param allele1 {@code String} allele 1 from source publication
   * @param allele2 {@code String} allele 2 from source publication
   * @param orientRs {@code AlleleOrientation} orientation of source publication
   *          alleles relative to refSNP alleles
   * @param allele2Freq {@code double} allele 2 frequency from source
   *          publication
   * @param allele2LnHR {@code double} Cox proportional hazards model ln hazard
   *          ratio for each additional allele 2 obtained from the source
   *          publication
   */
  public SNP(final String rsID, final String sourcePub, final String allele1,
    final String allele2, final AlleleOrientation orientRs,
    final double allele2Freq, final double allele2LnHR) {
    this.rsID = rsID;
    // Check that stored rsID has valid form.
    if (!this.rsID.matches("^rs[0-9]+$")) {
      throw new IllegalArgumentException("SNP Constructor: " + this.rsID +
        " is an invalid rs ID.");
    }
    this.sourcePub = sourcePub;
    this.allele1 = allele1.toUpperCase(Locale.US);
    this.allele2 = allele2.toUpperCase(Locale.US);
    // Regular expression describing form of all valid allele values.
    final String validAlleleRegex = "^-$|^[ACGT]+$";
    // Check that stored alleles match this regex.
    if (!(this.allele1.matches(validAlleleRegex) && this.allele2
      .matches(validAlleleRegex))) {
      throw new IllegalArgumentException("SNP Constructor: There is an " +
        "invalid input allele for SNP " + this.rsID + ". Valid input alleles " +
        "can be '-' or any string containing only the characters " +
        "'A', 'C', 'G', and 'T'.");
    }
    this.orientRs = orientRs;
    this.allele2Freq = allele2Freq;
    // Check that stored allele 2 frequency is in (0,1).
    if (!((this.allele2Freq > 0.0D) && (this.allele2Freq < 1.0D))) {
      throw new IllegalArgumentException("SNP Constructor: Allele 2 " +
        "frequency for SNP " + this.rsID + " is not in (0,1).");
    }
    this.allele2LnHR = allele2LnHR;
  }

  /** Returns {@code String} dbSNP refSNP identifier (i.e., rs number). */
  public final String getRsID() {
    return rsID;
  }

  /** Returns {@code String} source publication reference. */
  public final String getSourcePub() {
    return sourcePub;
  }

  /** Returns {@code String} allele 1 from source publication. */
  public final String getAllele1() {
    return allele1;
  }

  /** Returns {@code String} allele 2 from source publication. */
  public final String getAllele2() {
    return allele2;
  }

  /**
   * Returns {@code AlleleOrientation} orientation of source publication alleles
   * relative to refSNP alleles.
   */
  public final AlleleOrientation getOrientRs() {
    return orientRs;
  }

  /** Returns {@code double} allele 2 frequency from source publication. */
  public final double getAllele2Freq() {
    return allele2Freq;
  }

  /**
   * Returns {@code double} Cox proportional hazards model ln hazard ratio for
   * each additional allele 2 obtained from the source publication.
   */
  public final double getAllele2LnHR() {
    return allele2LnHR;
  }

  /**
   * Returns {@code double} ln probability of a given SNP genotype under the
   * assumption of HWE.
   * 
   * @param allele2Num Genotype, as {@code int} number of 2 alleles
   */
  public final double getLnProbGeno(final int allele2Num) {
    // @formatter:off
    /*
     * ln(p_2) = ln(allele2Freq^2) = 2*ln(allele2Freq) 
     * ln(p_1) = ln(2*allele2Freq*(1-allele2Freq)) 
     *         = ln(2) + ln(allele2Freq) + ln(1-allele2Freq)
     * ln(p_0) = ln((1-allele2Freq)^2) = 2*ln(1-allele2Freq)
     */
    // @formatter:on
    return ((allele2Num == 1) ? FastMath.log(2) : 0.0D) + allele2Num *
      FastMath.log(allele2Freq) + (2 - allele2Num) *
      FastMath.log(1 - allele2Freq);
  }

  /**
   * Returns {@code int} random SNP genotype, as the number of 2 alleles, under
   * the assumption of HWE using a user-supplied random number generator.
   * <p>
   * The implementation uses the fact that, under HWE, the SNP genotype (number
   * of 2 alleles) is distributed Binomial(2, {@code allele2Freq}), which is
   * also the distribution of the number of successes in two independent
   * Bernoulli trials each with success probability {@code allele2Freq}.
   * 
   * @param rng {@code RandomGenerator} object that provides the underlying
   *          pseudorandom number stream
   */
  public final int getRandGeno(final RandomGenerator rng) {
    /*
     * The nextDouble() RandomGenerator method produces an independent random
     * double drawn from U(0,1), so the expression ((rng.nextDouble() <
     * allele2Freq) ? 1 : 0) takes the value 1 with probability allele2Freq and
     * 0 with probability 1 - allele2Freq. The sum of two of these expressions
     * is therefore the number of successes in two independent Bernoulli trials
     * each with success probability allele2Freq.
     */
    return ((rng.nextDouble() < allele2Freq) ? 1 : 0) +
      ((rng.nextDouble() < allele2Freq) ? 1 : 0);
  }

  /**
   * Returns the {@code double} # of times the SNP object {@code allele2}
   * appears in the input genotype given by {@code inAllele1/inAllele2}, times
   * the SNP object per-allele2 ln hazard ratio. If the input genotype is two
   * missing ("0") alleles, then the expected value over all possible genotypes
   * at this SNP is returned.
   * 
   * @param inAllele1 {@code String} input allele1
   * @param inAllele2 {@code String} input allele2
   * @param inOrientRs {@code AlleleOrientation} orientation of input alleles
   *          relative to refSNP alleles.
   * 
   * @return {@code double} genotype score
   */
  public final double genoScore(final String inAllele1, final String inAllele2,
    final AlleleOrientation inOrientRs) {
    /*
     * Stored alleles in SNP object are in English locale uppercase (ensured by
     * constructor). Must convert input alleles to uppercase using English
     * locale for proper matching.
     */
    final String upperInAllele1 = inAllele1.toUpperCase(Locale.US);
    final String upperInAllele2 = inAllele2.toUpperCase(Locale.US);
    // Regular expression describing form of all valid allele values.
    final String validAlleleRegex = "^-$|^0$|^[ACGT]+$";
    // Check that uppercase input alleles match this regex.
    if (!(upperInAllele1.matches(validAlleleRegex) && upperInAllele2
      .matches(validAlleleRegex))) {
      throw new IllegalArgumentException("SNP.genoScore: There is an " +
        "invalid input allele for SNP " + rsID + ". Valid input alleles " +
        "can be '-', '0', or any string containing only the characters " +
        "'A', 'C', 'G', and 'T'.");
    }
    // Check that either both uppercase input alleles are "0" or neither is "0".
    if (upperInAllele1.equals("0") ^ upperInAllele2.equals("0")) {
      throw new IllegalArgumentException("SNP.genoScore: Neither or both " +
        "of the two input alleles should be '0' for SNP " + rsID + ".");
    }
    // Double to hold result.
    double result;
    if (upperInAllele1.equals("0") && upperInAllele2.equals("0")) {
      /*
       * If both uppercase input alleles are "0", the input genotype is missing,
       * so the result is the expected X_j*Beta_j over all possible genotypes at
       * this SNP.
       */
      result =
        allele2LnHR * FastMath.exp(getLnProbGeno(1)) + 2 * allele2LnHR *
          FastMath.exp(getLnProbGeno(2));
    } else {
      /*
       * Otherwise, process uppercase input alleles into appropriate
       * orientation.
       */
      String orientedInAllele1;
      String orientedInAllele2;
      // Use != for safer enum comparison.
      if (inOrientRs != orientRs) {
        /*
         * If the uppercase input alleles are in the opposite orientation to the
         * SNP alleles, translate A -> T, C -> G, G -> C, and T -> A in each
         * uppercase input allele string.
         */
        orientedInAllele1 =
          StringUtils.replaceChars(upperInAllele1, "ACGT", "TGCA");
        orientedInAllele2 =
          StringUtils.replaceChars(upperInAllele2, "ACGT", "TGCA");
      } else {
        // Otherwise, copy uppercase input alleles.
        orientedInAllele1 = upperInAllele1;
        orientedInAllele2 = upperInAllele2;
      }
      /*
       * Create regex reflecting possible SNP alleles and check that oriented
       * uppercase input allele strings both match this regex.
       */
      final String possibleAllelesRegex = "^" + allele1 + "$|^" + allele2 + "$";
      if (!(orientedInAllele1.matches(possibleAllelesRegex) && orientedInAllele2
        .matches(possibleAllelesRegex))) {
        throw new IllegalArgumentException("One or both input alleles " +
          "differ from the possible population alleles for SNP " + rsID + ".");
      }
      /*
       * Obtain the number of times allele2 appears in the oriented uppercase
       * input alleles by counting the number of occurrences of the allele2
       * substring in the concatenation of the uppercase input allele strings.
       * The result is this number times the per-allele2 ln hazard ratio.
       */
      result =
        StringUtils
          .countMatches(orientedInAllele1 + orientedInAllele2, allele2) *
          allele2LnHR;
    }
    return result;
  }

  @Override
  public final String toString() {
    return rsID + "\t" + sourcePub + "\t" + allele1 + "/" + allele2 + "\t" +
      orientRs + "\t" + allele2Freq + "\t" + allele2LnHR;
  }

}
