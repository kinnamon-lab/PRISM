/*
 * This source file is part of the PRISM software package.
 *
 * Copyright 2014-2017 The Ohio State University Wexner Medical Center
 *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except
 * in compliance with the License. You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software distributed under the License
 * is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express
 * or implied. See the License for the specific language governing permissions and limitations under
 * the License.
 */

package edu.osumc.prism;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Locale;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.tuple.ImmutableTriple;

/**
 * Individual data wrapper class.
 */
public final class Individual {

  /**
   * Individual genotype data.
   * <p>
   * The constructor, setters, and getters are part of the API so that users can package individual
   * genotype data into a single {@code Genotype} object and pass it to the
   * {@link RiskModel#getRiskPrediction(Genotypes) RiskModel.getRiskPrediction} method to obtain
   * individual risk predictions.
   */
  public static final class Genotypes {

    /** {@code String} individual identifier. */
    private final String indivID;
    /**
     * {@code HashMap<String, ImmutableTriple<String, String, SNP.AlleleOrientation>>} containing
     * {@link ImmutableTriple} values storing individual genotype data (one string for each allele
     * and a {@link SNP.AlleleOrientation} constant) for each SNP identified by a {@code String}
     * dbSNP refSNP identifier key.
     */
    private final HashMap<String, ImmutableTriple<String, String, SNP.AlleleOrientation>> genotypeRsMap;

    /**
     * Constructs a {@code Genotypes} object ready to receive genotypes for an individual.
     *
     * @param indivID {@code String} individual identifier
     */
    public Genotypes(final String indivID) {
      this.indivID = indivID;
      genotypeRsMap = new HashMap<String, ImmutableTriple<String, String, SNP.AlleleOrientation>>();
    }

    /**
     * Adds genotype value for a given SNP.
     *
     * @param rsID {@code String} dbSNP refSNP identifier (i.e., rs number)
     * @param allele1 {@code String} allele 1 at SNP for individual
     * @param allele2 {@code String} allele 2 at SNP for individual
     * @param orientRsStr {@code String} allele orientation values (relative to RefSNP) for SNP
     */
    public final void addGenotype(final String rsID, final String allele1, final String allele2,
        final String orientRsStr) {
      // Check that input rsID has valid form.
      if (!rsID.matches("^rs[0-9]+$")) {
        throw new IllegalArgumentException(
            "Genotypes.addGenotype: " + rsID + " is an invalid rs ID.");
      }
      // Regular expression describing form of all valid allele values.
      final String validAlleleRegex = "^-$|^0$|^[AaCcGgTt]+$";
      // Check that input alleles match this regex.
      if (!(allele1.matches(validAlleleRegex) && allele2.matches(validAlleleRegex))) {
        throw new IllegalArgumentException("Genotypes.addGenotype: There "
            + "is an invalid input allele for SNP " + rsID + " in individual " + indivID
            + ". Valid input alleles can be '-', '0', or any string "
            + "containing only the characters 'A', 'C', 'G', and 'T' in upper " + "or lower case.");
      }
      // Check that either both input alleles are "0" or neither is "0".
      if (allele1.equals("0") ^ allele2.equals("0")) {
        throw new IllegalArgumentException("Genotypes.addGenotype: Neither "
            + "or both of the two input alleles should be '0' for SNP " + rsID + " in individual "
            + indivID + ".");
      }
      // Check that input allele orientation string is valid.
      if (!orientRsStr.toUpperCase(Locale.US).matches("^(FORWARD|REVERSE)$")) {
        throw new IllegalArgumentException("Genotypes.addGenotype: SNP " + rsID
            + " must have an allele orientation relative to RefSNP of "
            + "'Forward' or 'Reverse' provided.");
      }
      final SNP.AlleleOrientation orientRs = orientRsStr.toUpperCase(Locale.US).matches("FORWARD")
          ? SNP.AlleleOrientation.FORWARD : SNP.AlleleOrientation.REVERSE;
      genotypeRsMap.put(rsID,
          new ImmutableTriple<String, String, SNP.AlleleOrientation>(allele1, allele2, orientRs));
    }

    /** Returns {@code String} individual identifier. */
    public final String getIndivID() {
      return indivID;
    }

    /**
     * Returns {@code String} allele 1 for the SNP with identifier {@code rsID} if it is available
     * and the missing allele code "0" otherwise.
     *
     * @param rsID {@code String} dbSNP refSNP identifier (i.e., rs number)
     */
    public final String getAllele1(final String rsID) {
      return genotypeRsMap.containsKey(rsID) ? genotypeRsMap.get(rsID).getLeft() : "0";
    }

    /**
     * Returns {@code String} allele 2 for the SNP with identifier {@code rsID} if it is available
     * and the missing allele code "0" otherwise.
     *
     * @param rsID {@code String} dbSNP refSNP identifier (i.e., rs number)
     */
    public final String getAllele2(final String rsID) {
      return genotypeRsMap.containsKey(rsID) ? genotypeRsMap.get(rsID).getMiddle() : "0";
    }

    /**
     * Returns {@link SNP.AlleleOrientation} allele orientation (relative to RefSNP) for the SNP
     * with identifier {@code rsID} if it is available and {@code null} otherwise.
     *
     * @param rsID {@code String} dbSNP refSNP identifier (i.e., rs number)
     */
    public final SNP.AlleleOrientation getOrientRs(final String rsID) {
      return genotypeRsMap.containsKey(rsID) ? genotypeRsMap.get(rsID).getRight() : null;
    }

  }

  /**
   * Individual risk prediction from a particular risk model.
   * <p>
   * Getters are part of the API, but the constructor and setters can only be called from within
   * this package. This ensures that {@code RiskPrediction} object state can only be modified within
   * methods of other package classes, namely the {@link RiskModel#getRiskPrediction(Genotypes)
   * RiskModel.getRiskPrediction} method. This obviates the need for constructor/setter argument
   * checking.
   */
  public static final class RiskPrediction {

    /** {@code String} individual identifier. */
    private final String indivID;
    /** {@code String} risk model name. */
    private final String modelName;
    /**
     * {@code LinkedHashMap<String, String>} containing {@code String} genotype values (in the form
     * "allele1/allele2") used to generate the risk prediction stored in this object. Used SNP
     * genotypes are identified by their {@code String} dbSNP refSNP identifier keys and are stored
     * in a {@code LinkedHashMap} so that the iterator always extracts SNPs in the same order that
     * they were added.
     */
    private final LinkedHashMap<String, String> usedGenotypesRsMap;
    /**
     * {@code double} prognostic index for an individual.
     */
    private double prognosticIndex;
    /**
     * {@code double} prognostic index percentile.
     */
    private double prognosticIndexPctl;
    /**
     * {@code ArrayList<Double>} list of of times at which the the predicted cumulative risk
     * function is evaluated (abscissae).
     */
    private final ArrayList<Double> timesList;
    /**
     * {@code ArrayList<Double>} list of predicted cumulative risks for an individual, where
     * {@code predCumRiskList[i]} is the individual's predicted cumulative risk at
     * {@code timesList[i]}.
     */
    private final ArrayList<Double> predCumRiskList;

    /**
     * Constructs a {@code RiskPrediction} object ready to receive predictions for a given
     * individual from a particular risk model.
     *
     * @param indivID {@code String} individual identifier
     * @param modelName {@code String} risk model name
     */
    RiskPrediction(final String indivID, final String modelName) {
      this.indivID = indivID;
      this.modelName = modelName;
      usedGenotypesRsMap = new LinkedHashMap<String, String>();
      timesList = new ArrayList<Double>();
      predCumRiskList = new ArrayList<Double>();
    }

    /**
     * Adds a used genotype to the end of the {@code LinkedHashMap} object member.
     *
     * @param rsID {@code String} dbSNP refSNP identifier (i.e., rs number)
     * @param allele1 {@code String} used allele 1 at SNP for individual
     * @param allele2 {@code String} used allele 2 at SNP for individual
     */
    final void addUsedGenotype(final String rsID, final String allele1, final String allele2) {
      usedGenotypesRsMap.put(rsID, allele1 + "/" + allele2);
    }

    /**
     * Sets prognostic index field values.
     *
     * @param prognosticIndex {@code double} prognostic index for an individual
     * @param prognosticIndexPctl {@code double} prognostic index percentile
     */
    final void setPrognosticIndex(final double prognosticIndex, final double prognosticIndexPctl) {
      this.prognosticIndex = prognosticIndex;
      this.prognosticIndexPctl = prognosticIndexPctl;
    }

    /**
     * Adds a time to the end of the {@code timesList} object member and the predicted cumulative
     * risk at this time to the end of the {@code predCumRiskList} object member.
     *
     * @param t {@code double} time at which the the predicted cumulative risk function is evaluated
     *        (abscissa)
     * @param predCumRiskT {@code double} predicted cumulative risk for an individual at time
     *        {@code t}
     */
    final void addPredCumRiskT(final double t, final double predCumRiskT) {
      timesList.add(t);
      predCumRiskList.add(predCumRiskT);
    }

    /** Returns {@code String} individual identifier. */
    public final String getIndivID() {
      return indivID;
    }

    /** Returns {@code String} risk model name. */
    public final String getModelName() {
      return modelName;
    }

    /**
     * Returns {@code LinkedHashMap<String, String>} containing {@code String} genotype values (in
     * the form "allele1/allele2") used to generate the risk prediction stored in this object. Used
     * SNP genotypes are identified by their {@code String} dbSNP refSNP identifier keys and are
     * stored in a {@code LinkedHashMap} so that the iterator always extracts SNPs in the same order
     * that they were added.
     */
    public final LinkedHashMap<String, String> getUsedGenotypesRsMap() {
      return usedGenotypesRsMap;
    }

    /**
     * Returns {@code double} prognostic index for an individual.
     */
    public final double getPrognosticIndex() {
      return prognosticIndex;
    }

    /**
     * Returns {@code double} prognostic index percentile.
     */
    public final double getPrognosticIndexPctl() {
      return prognosticIndexPctl;
    }

    /**
     * Returns {@code double[]} array of of times at which the the predicted cumulative risk
     * function is evaluated (abscissae).
     */
    public final double[] getTimes() {
      return ArrayUtils.toPrimitive(timesList.toArray(new Double[0]));
    }

    /**
     * Returns {@code double[]} array of predicted cumulative risks for an individual, where
     * {@code getPredCumRisk()[i]} is the individual's predicted cumulative risk at
     * {@code getTimes()[i]}.
     */
    public final double[] getPredCumRisk() {
      return ArrayUtils.toPrimitive(predCumRiskList.toArray(new Double[0]));
    }

  }

  /** Prevents construction of {@code Individual} class instances. */
  private Individual() {}

}
