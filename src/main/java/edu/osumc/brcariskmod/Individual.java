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
import java.util.HashMap;
import java.util.LinkedHashMap;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.tuple.ImmutablePair;

/**
 * Individual data wrapper class.
 * 
 * @author Daniel Kinnamon
 * @version 2014-10-26
 * @since 2014-10-26
 */
public final class Individual {

  /** Individual genotype data. */
  public static final class Genotypes {

    /** {@code String} individual identifier. */
    private final String indivID;
    /**
     * {@code HashMap<String, SNP.AlleleOrientation>} containing
     * {@link SNP.AlleleOrientation} allele orientation values (relative to
     * RefSNP) for SNPs identified by {@code String} dbSNP refSNP identifier
     * keys.
     */
    private final HashMap<String, SNP.AlleleOrientation> orientRsMap;
    /**
     * {@code HashMap<String, MutablePair<String, String>>} containing
     * {@code MutablePair<String, String>} genotype values (one string for each
     * allele) for SNPs identified by {@code String} dbSNP refSNP identifier
     * keys.
     */
    private final HashMap<String, ImmutablePair<String, String>> genotypeRsMap;

    /**
     * Constructs a {@code Genotypes} object ready to receive genotypes for an
     * individual.
     * 
     * @param indivID {@code String} individual identifier
     * @param orientRsMap {@code HashMap<String, SNP.AlleleOrientation>}
     *          containing {@link SNP.AlleleOrientation} allele orientation
     *          values (relative to RefSNP) for SNPs identified by
     *          {@code String} dbSNP refSNP identifier keys
     */
    public Genotypes(final String indivID,
      final HashMap<String, SNP.AlleleOrientation> orientRsMap) {
      this.indivID = indivID;
      this.orientRsMap = orientRsMap;
      genotypeRsMap = new HashMap<String, ImmutablePair<String, String>>();
    }

    /**
     * Sets genotype value for a given SNP.
     * 
     * @param rsID {@code String} dbSNP refSNP identifier (i.e., rs number)
     * @param allele1 {@code String} allele 1 at SNP for individual
     * @param allele2 {@code String} allele 2 at SNP for individual
     */
    public final void setGenotype(final String rsID, final String allele1,
      final String allele2) {
      genotypeRsMap.put(rsID, new ImmutablePair<String, String>(allele1,
        allele2));
    }

    /** Returns {@code String} individual identifier. */
    public final String getIndivID() {
      return indivID;
    }

    /**
     * Returns {@code String} allele 1 for the SNP with identifier {@code rsID}
     * if it is available and the missing allele code "0" otherwise.
     * 
     * @param rsID {@code String} dbSNP refSNP identifier (i.e., rs number)
     */
    public final String getAllele1(final String rsID) {
      return genotypeRsMap.containsKey(rsID) ? genotypeRsMap.get(rsID)
        .getLeft() : "0";
    }

    /**
     * Returns {@code String} allele 2 for the SNP with identifier {@code rsID}
     * if it is available and the missing allele code "0" otherwise.
     * 
     * @param rsID {@code String} dbSNP refSNP identifier (i.e., rs number)
     */
    public final String getAllele2(final String rsID) {
      return genotypeRsMap.containsKey(rsID) ? genotypeRsMap.get(rsID)
        .getRight() : "0";
    }

    /**
     * Returns {@code SNP.AlleleOrientation} allele orientation for the SNP with
     * identifier {@code rsID} if it is available and {@code null} otherwise.
     * 
     * @param rsID {@code String} dbSNP refSNP identifier (i.e., rs number)
     */
    public final SNP.AlleleOrientation getOrientRs(final String rsID) {
      return orientRsMap.containsKey(rsID) ? orientRsMap.get(rsID) : null;
    }

  }

  /** Individual risk prediction from a particular risk model. */
  public static final class RiskPrediction {

    /** {@code String} individual identifier. */
    private final String indivID;
    /** {@code String} risk model name. */
    private final String modelName;
    /**
     * {@code LinkedHashMap<String, String>} containing {@code String} genotype
     * values (in the form "allele1/allele2") used to generate the risk
     * prediction stored in this object. Used SNP genotypes are identified by
     * their {@code String} dbSNP refSNP identifier keys and are stored in a
     * {@code LinkedHashMap} so that the iterator always extracts SNPs in the
     * same order that they were added.
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
     * {@code ArrayList<Double>} list of of times at which the the predicted
     * cumulative risk function is evaluated (abscissae).
     */
    private final ArrayList<Double> timesList;
    /**
     * {@code ArrayList<Double>} list of predicted cumulative risks for an
     * individual, where {@code predCumRiskList[i]} is the individual's
     * predicted cumulative risk at {@code timesList[i]}.
     */
    private final ArrayList<Double> predCumRiskList;

    /**
     * Constructs a {@code RiskPrediction} object ready to receive predictions
     * for a given individual from a particular risk model.
     * 
     * @param indivID {@code String} individual identifier
     * @param modelName {@code String} risk model name
     */
    public RiskPrediction(final String indivID, final String modelName) {
      this.indivID = indivID;
      this.modelName = modelName;
      usedGenotypesRsMap = new LinkedHashMap<String, String>();
      timesList = new ArrayList<Double>();
      predCumRiskList = new ArrayList<Double>();
    }

    /**
     * Adds a used genotype to the end of the {@code LinkedHashMap} object
     * member.
     * 
     * @param rsID {@code String} dbSNP refSNP identifier (i.e., rs number)
     * @param allele1 {@code String} used allele 1 at SNP for individual
     * @param allele2 {@code String} used allele 2 at SNP for individual
     */
    public final void addUsedGenotype(final String rsID, final String allele1,
      final String allele2) {
      usedGenotypesRsMap.put(rsID, allele1 + "/" + allele2);
    }

    /**
     * Sets prognostic index field values.
     * 
     * @param prognosticIndex {@code double} prognostic index for an individual
     * @param prognosticIndexPctl {@code double} prognostic index percentile
     */
    public final void setPrognosticIndex(final double prognosticIndex,
      final double prognosticIndexPctl) {
      this.prognosticIndex = prognosticIndex;
      this.prognosticIndexPctl = prognosticIndexPctl;
    }

    /**
     * Adds a time to the end of the {@code timesList} object member and the
     * predicted cumulative risk at this time to the end of the
     * {@code predCumRiskList} object member.
     * 
     * @param t {@code double} time at which the the predicted cumulative risk
     *          function is evaluated (abscissa)
     * @param predCumRiskT {@code double} predicted cumulative risk for an
     *          individual at time {@code t}
     */
    public final void
      addPredCumRiskT(final double t, final double predCumRiskT) {
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
     * Returns {@code LinkedHashMap<String, String>} containing {@code String}
     * genotype values (in the form "allele1/allele2") used to generate the risk
     * prediction stored in this object. Used SNP genotypes are identified by
     * their {@code String} dbSNP refSNP identifier keys and are stored in a
     * {@code LinkedHashMap} so that the iterator always extracts SNPs in the
     * same order that they were added.
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
     * Returns {@code double[]} array of of times at which the the predicted
     * cumulative risk function is evaluated (abscissae).
     */
    public final double[] getTimes() {
      return ArrayUtils.toPrimitive(timesList.toArray(new Double[0]));
    }

    /**
     * Returns {@code double[]} array of predicted cumulative risks for an
     * individual, where {@code getPredCumRisk()[i]} is the individual's
     * predicted cumulative risk at {@code getTimes()[i]}.
     */
    public final double[] getPredCumRisk() {
      return ArrayUtils.toPrimitive(predCumRiskList.toArray(new Double[0]));
    }

  }

  /** Prevents construction of {@code Individual} class instances. */
  private Individual() {
  }

}
