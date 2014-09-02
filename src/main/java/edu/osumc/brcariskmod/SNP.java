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
 * @version 2014-09-01
 * @since 2014-09-01
 */
public class SNP implements Serializable {

  /** Allele orientation values for a SNP */
  public enum AlleleOrientation {
    FORWARD, REVERSE
  }

  /**
   * The SUID is given by the class modification date in YYYYMMDD format,
   * followed by committing author initials as 2-digit alphabet numbers in
   * FFMMLL format (e.g., initials DDK are 040411), followed by the 1-digit
   * revision number by that author on that date. Most of the time, this
   * revision number should be 1.
   */
  private static final long serialVersionUID = 201409010404111L;
  /** @serial String dbSNP refSNP identifier (i.e., rs number). */
  private String rsID;
  /** @serial String source publication reference. */
  private String sourcePub;
  /** @serial String allele 1 from source publication. */
  private String allele1;
  /** @serial String allele 2 from source publication. */
  private String allele2;
  /**
   * @serial AlleleOrientation orientation of source publication alleles
   *         relative to refSNP alleles.
   */
  private AlleleOrientation orientRs;
  /** @serial double allele 2 frequency from source publication. */
  private double allele2Freq;
  /**
   * @serial double Cox proportional hazards model ln hazard ratio for each
   *         additional allele 2 obtained from the source publication.
   */
  private double allele2LnHR;

  /**
   * Constructs a SNP object.
   * 
   * @param rsID String dbSNP refSNP identifier (i.e., rs number)
   * @param sourcePub String source publication reference
   * @param allele1 String allele 1 from source publication
   * @param allele2 String allele 2 from source publication
   * @param orientRs AlleleOrientation orientation of source publication alleles
   *          relative to refSNP alleles
   * @param allele2Freq double allele 2 frequency from source publication
   * @param allele2LnHR double Cox proportional hazards model ln hazard ratio
   *          for each additional allele 2 obtained from the source publication
   */
  public SNP(String rsID, String sourcePub, String allele1, String allele2,
    AlleleOrientation orientRs, double allele2Freq, double allele2LnHR) {
    this.rsID = rsID;
    this.sourcePub = sourcePub;
    this.allele1 = allele1;
    this.allele2 = allele2;
    this.orientRs = orientRs;
    this.allele2Freq = allele2Freq;
    this.allele2LnHR = allele2LnHR;
  }

  /** Returns String dbSNP refSNP identifier (i.e., rs number). */
  public String getRsID() {
    return rsID;
  }

  /** Returns String source publication reference. */
  public String getSourcePub() {
    return sourcePub;
  }

  /** Returns String allele 1 from source publication. */
  public String getAllele1() {
    return allele1;
  }

  /** Returns String allele 2 from source publication. */
  public String getAllele2() {
    return allele2;
  }

  /**
   * Returns AlleleOrientation orientation of source publication alleles
   * relative to refSNP alleles.
   */
  public AlleleOrientation getOrientRs() {
    return orientRs;
  }

  /** Returns double allele 2 frequency from source publication. */
  public double getAllele2Freq() {
    return allele2Freq;
  }

  /**
   * Returns double Cox proportional hazards model ln hazard ratio for each
   * additional allele 2 obtained from the source publication.
   */
  public double getAllele2LnHR() {
    return allele2LnHR;
  }

  @Override
  public String toString() {
    return rsID + "\t" + sourcePub + "\t" + allele1 + "/" + allele2 + "\t" +
      orientRs + "\t" + allele2Freq + "\t" + allele2LnHR;
  }
}
