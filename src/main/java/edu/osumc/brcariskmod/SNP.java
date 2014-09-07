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
 * @version 2014-09-07
 * @since 2014-09-01
 */
public final class SNP implements Serializable {

  /** Allele orientation values for a SNP */
  public static enum AlleleOrientation {
    FORWARD, REVERSE
  }

  /**
   * The SUID is given by the class modification date in YYYYMMDD format,
   * followed by committing author initials as 2-digit alphabet numbers in
   * FFMMLL format (e.g., initials DDK are 040411), followed by the 1-digit
   * revision number by that author on that date. Most of the time, this
   * revision number should be 1.
   */
  private static final long serialVersionUID = 201409070404111L;
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

  @Override
  public final String toString() {
    return rsID + "\t" + sourcePub + "\t" + allele1 + "/" + allele2 + "\t" +
      orientRs + "\t" + allele2Freq + "\t" + allele2LnHR;
  }
}
