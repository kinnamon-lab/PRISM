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
 * @version 2014-09-23
 * @since 2014-09-01
 */
public final class SNP implements Serializable {

  /** Allele orientation values for a SNP */
  public static enum AlleleOrientation {
    FORWARD, REVERSE
  }

  /**
   * The {@code serialVersionUID} should be incremented by 1 every time an
   * incompatible change (for serialization) is made to the class. If all
   * changes are compatible, then the {@code serialVersionUID} should not be
   * changed.
   */
  private static final long serialVersionUID = 1L;
  /** @serial {@code String} dbSNP refSNP identifier (i.e., rs number). */
  private final String rsID;
  /** @serial {@code String} source publication reference. */
  private final String sourcePub;
  /** @serial {@code String} allele 1 from source publication. */
  private final String allele1;
  /** @serial {@code String} allele 2 from source publication. */
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
    this.sourcePub = sourcePub;
    this.allele1 = allele1;
    this.allele2 = allele2;
    this.orientRs = orientRs;
    this.allele2Freq = allele2Freq;
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

  @Override
  public final String toString() {
    return rsID + "\t" + sourcePub + "\t" + allele1 + "/" + allele2 + "\t" +
      orientRs + "\t" + allele2Freq + "\t" + allele2LnHR;
  }
}
