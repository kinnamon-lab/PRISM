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

import java.io.BufferedReader;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Scanner;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.Precision;

/**
 * Constructs and serializes {@link RiskModel} objects during the build process.
 * <p>
 * This package-private class contains a private nested helper class, several
 * private static methods, and a {@code main} method called during the build
 * process. Because it is used only in the build process, it is neither part of
 * the package API nor included in the final artifact.
 * 
 * @author Daniel Kinnamon
 * @version 2014-10-30
 * @since 2014-10-30
 */
final class RiskModelBuilder {

  /** Holds raw data obtained from input files. */
  private static final class RiskModelRawData {

    /** {@code String} risk model name. */
    private final String modelName;
    /** {@code ArrayList<SNP>} list of model SNPs. */
    private final ArrayList<SNP> modelSNPsList;
    /**
     * {@code ArrayList<Integer>} list of of consecutive annual ages.
     */
    private final ArrayList<Integer> ageYrsList;
    /**
     * {@code ArrayList<Double>} list of annual incidences, where
     * {@code annIncList[i]} is the annual incidence in the year preceding age
     * {@code ageYrsList[i]}.
     */
    private final ArrayList<Double> annIncList;

    /**
     * Constructs a {@code RiskModelRawData} object ready to receive data.
     * 
     * @param modelName {@code String} risk model name
     */
    private RiskModelRawData(final String modelName) {
      this.modelName = modelName;
      modelSNPsList = new ArrayList<SNP>();
      ageYrsList = new ArrayList<Integer>();
      annIncList = new ArrayList<Double>();
    }

    /**
     * Adds SNP to the end of the {@code modelSNPList} object member. SNPs can
     * be added in any order for a given risk model.
     */
    private final void addSNP(final SNP newSNP) {
      modelSNPsList.add(newSNP);
    }

    /**
     * Adds a consecutive yearly age to the end of the {@code ageYrsList} object
     * member and the annual incidence in the year preceding this age to the end
     * of the {@code annIncList} object member. Performs checks to make sure
     * that annual incidences for consecutive yearly ages are being added in
     * order and that {@code ageYrsList} and {@code annIncList} are growing at
     * the same rate.
     * 
     * @param ageYrs {@code int} age in years
     * @param annInc {@code double} annual incidence in the year preceding
     *          {@code ageYrs}
     */
    private final void addAnnInc(final int ageYrs, final double annInc) {
      // Check for an annual incidence of 0 at age 0.
      if (ageYrs == 0 && !Precision.equals(annInc, 0.0D)) {
        throw new IllegalArgumentException("Annual incidence in the year "
          + "preceding birth (age 0) must be 0 so that S(0) = 1. Please "
          + "check the annual incidence input file");
      }
      // Check that all annual incidences are non-negative.
      if (annInc < 0.0D) {
        throw new IllegalArgumentException("Annual incidences must be "
          + "non-negative. Please check the annual incidence input file");
      }
      ageYrsList.add(ageYrs);
      annIncList.add(annInc);
      /*
       * Check to make sure that ageYrsList[ageYrs] = ageYrs and
       * annIncList[ageYrs] = annInc.
       */
      if (ageYrsList.get(ageYrs) != ageYrs ||
        !Precision.equals(annIncList.get(ageYrs), annInc)) {
        throw new IllegalArgumentException("Annual incidences for "
          + "consecutive yearly ages starting at 0 must be supplied in "
          + "order. Please check the annual incidence input file.");
      }
      // Check that lists are growing at the same rate.
      if (ageYrsList.size() != annIncList.size()) {
        throw new IllegalArgumentException("ageYrsList and annIncList are "
          + "growing at different rates. Please check the annual incidence "
          + " input file.");
      }
    }

    /** Returns {@code String} risk model name. */
    private final String getModelName() {
      return modelName;
    }

    /** Returns {@code SNP[]} array of model SNPs. */
    private final SNP[] getModelSNPs() {
      return modelSNPsList.toArray(new SNP[0]);
    }

    /**
     * Returns {@code int[]} array of of ages at which the the marginal survivor
     * function is evaluated (abscissae).
     */
    private final int[] getAgeYrs() {
      return ArrayUtils.toPrimitive(ageYrsList.toArray(new Integer[0]));
    }

    /**
     * Returns {@code double[]} array of marginal survivor function values,
     * where {@code getMargSurv()[i]} is the marginal survivor function at
     * {@code getAgeYrs()[i]}.
     */
    private final double[] getMargSurv() {
      ArrayList<Double> margSurvList = new ArrayList<Double>();
      /*
       * Use the annual incidences to calculate the continuous-time survivor
       * function assuming a piecewise constant yearly hazard (i.e., annual
       * incidence). Because of checks in the addAnnInc method, annIncList[0] =
       * 0 and the indices in ageYrsList and annIncList both correspond to
       * annual ages. The cumulative hazard from age 0 to age i, H(i), is
       * therefore given by sum_{j=0}^{i} annIncList[j], or the cumulative sum
       * over elements 0 to i of annIncList. Calculating the cumulative sum at
       * each index i therefore gives the cumulative hazard to age i, and the
       * relationship S(i) = exp(-H(i)) is used to obtain the survivor function
       * at age i, which is stored in margSurvList[i].
       */
      double cumHazI = 0.0D;
      for (int i = 0; i < annIncList.size(); i++) {
        cumHazI += annIncList.get(i);
        final double survI = FastMath.exp(-cumHazI);
        margSurvList.add(survI);
        /*
         * Check that S(i) was added to margSurvList[i] so that S(ageYrsList[i])
         * = S(i) = margSurvList[i].
         */
        if (!Precision.equals(margSurvList.get(i), survI)) {
          throw new RuntimeException("Problem with calculating survivor "
            + "function from annual incidences.");
        }
      }
      return ArrayUtils.toPrimitive(margSurvList.toArray(new Double[0]));
    }
  }

  /** Prevents construction of {@code RiskModelBuilder} class instances. */
  private RiskModelBuilder() {
  }

  /**
   * Parses model SNPs and annual incidence files and returns a {@link HashMap}
   * linking {@code String} model ID keys to {@code RiskModelRawData} values
   * that can be used to build a {@link RiskModel} object.
   * 
   * @param modelSNPsFile {@code String} path to properly formatted model SNPs
   *          file
   * @param annIncFile {@code String} path to properly formatted annual
   *          incidence file
   */
  private static final HashMap<String, RiskModelRawData> parseInputFiles(
    final String modelSNPsFile, final String annIncFile) throws IOException {
    HashMap<String, RiskModelRawData> modelRawDataMap =
      new HashMap<String, RiskModelRawData>();
    // Parse model SNPs file.
    try (
      BufferedReader fileReader =
        Files.newBufferedReader(Paths.get(modelSNPsFile),
          StandardCharsets.UTF_8)) {
      String line;
      int linesRead = 0;
      // Read lines from file...
      while ((line = fileReader.readLine()) != null) {
        // For each line, try to create a scanner using tab delimiter.
        try (
          Scanner lineScanner =
            new Scanner(line).useLocale(Locale.US).useDelimiter("\t")
              .useRadix(10)) {
          // Parse tab-delimited tokens from the current line.
          if (linesRead == 0) {
            // If first line, check that file column headers are correct.
            List<String> expColHeaderList =
              Arrays.asList("BRCAType", "cancerType", "rsID", "sourcePub",
                "allele1", "allele2", "orientRs", "allele2Freq", "allele2HR");
            for (String expColHeader : expColHeaderList) {
              if (!expColHeader.equals(lineScanner.next())) {
                throw new IOException("Model SNPs file column headers "
                  + "are incorrect.");
              }
            }
          } else {
            // Otherwise, read SNP data.
            final short brcaType = lineScanner.nextShort();
            final String cancerType = lineScanner.next();
            final String rsID = lineScanner.next();
            final String sourcePub = lineScanner.next();
            final String allele1 = lineScanner.next();
            final String allele2 = lineScanner.next();
            final String orientRsStr =
              lineScanner.next().toUpperCase(Locale.US);
            final double allele2Freq = lineScanner.nextDouble();
            final double allele2lnHR = FastMath.log(lineScanner.nextDouble());
            // Create modelID.
            final String modelID =
              "brca" + brcaType + cancerType.toLowerCase(Locale.US);
            /*
             * Check that string input is parsed properly (alleles are checked
             * by the SNP constructor).
             */
            if (!modelID.matches("brca[12](breast|ovarian)")) {
              throw new IOException("Model SNPs file BRCA mutation type must "
                + "be '1' or '2' and cancer type must be 'Breast' or "
                + "'Ovarian'.");
            }
            if (!orientRsStr.matches("FORWARD|REVERSE")) {
              throw new IOException("Model SNPs file allele orientation "
                + "relative to RefSNP must be 'Forward' or 'Reverse'.");
            }
            // If new modelID, construct map entry.
            if (!modelRawDataMap.containsKey(modelID)) {
              modelRawDataMap.put(modelID, new RiskModelRawData("BRCA " +
                brcaType + " " + cancerType + " Cancer"));
            }
            // Add SNP to RiskModelRawData object for modelID.
            modelRawDataMap.get(modelID).addSNP(
              new SNP(rsID, sourcePub, allele1, allele2, orientRsStr
                .matches("FORWARD") ? SNP.AlleleOrientation.FORWARD
                : SNP.AlleleOrientation.REVERSE, allele2Freq, allele2lnHR));
          }
          // Check that there are no input tokens other than the ones expected.
          if (lineScanner.hasNext()) {
            throw new IOException("Input tokens available past expected "
              + "end of line in model SNPs file.");
          }
        } // Scanner automatically closed here.
        // Increment number of model SNPs file lines read.
        linesRead++;
      }
    } // Model SNPs file BufferedReader automatically closed here.
    // Parse annual Incidence file.
    try (
      BufferedReader fileReader =
        Files.newBufferedReader(Paths.get(annIncFile), StandardCharsets.UTF_8)) {
      String line;
      int linesRead = 0;
      // Read lines from file...
      while ((line = fileReader.readLine()) != null) {
        // For each line, try to create a scanner using tab delimiter.
        try (
          Scanner lineScanner =
            new Scanner(line).useLocale(Locale.US).useDelimiter("\t")
              .useRadix(10)) {
          // Parse tab-delimited tokens from the current line.
          if (linesRead == 0) {
            // If first line, check that column headers are correct.
            List<String> expColHeaderList =
              Arrays.asList("BRCAType", "cancerType", "ageYrs", "annInc");
            for (String expColHeader : expColHeaderList) {
              if (!expColHeader.equals(lineScanner.next())) {
                throw new IOException("Annual incidence file column headers "
                  + "are incorrect.");
              }
            }
          } else {
            // Otherwise, read annual incidence data.
            final short brcaType = lineScanner.nextShort();
            final String cancerType = lineScanner.next();
            int ageYrs = lineScanner.nextInt();
            double annInc = lineScanner.nextDouble();
            // Get modelID.
            final String modelID =
              "brca" + brcaType + cancerType.toLowerCase(Locale.US);
            /*
             * If the model ID for the current line corresponds to a model for
             * which we have SNPs, then process the current line.
             */
            if (modelRawDataMap.containsKey(modelID)) {
              /*
               * Annual incidence file rows should be in order of consecutive
               * annual ages within each modelID. This is checked by the
               * addAnnInc method.
               */
              modelRawDataMap.get(modelID).addAnnInc(ageYrs, annInc);
            }
          }
          // Check that there are no input tokens other than the ones expected.
          if (lineScanner.hasNext()) {
            throw new IOException("Input tokens available past expected "
              + "end of line in annual incidence file.");
          }
        } // Scanner automatically closed here.
        // Increment number of annual incidence file lines read.
        linesRead++;
      }
    } // Annual incidence file BufferedReader automatically closed here.
    // After files are processed, return modelRawDataMap.
    return modelRawDataMap;
  }

  /**
   * Constructs a {@link RiskModel} object for each model with a
   * {@code RiskModelRawData} entry in the input {@link HashMap} and saves
   * serialized version (*.rmo file) to the specified path.
   * 
   * @param savePath {@code String} path in which to save *.rmo files
   * @param modelRawDataMap {@code HashMap<String,
   *          RiskModelRawData>} with each entry containing a {@code String}
   *          model ID key and a {@code RiskModelRawData} value containing input
   *          data parsed from input files
   */
  private static final void buildRiskModels(final String savePath,
    final HashMap<String, RiskModelRawData> modelRawDataMap) throws IOException {
    // For each modelRawDataMap entry...
    for (Map.Entry<String, RiskModelRawData> modelRawDataEntry : modelRawDataMap
      .entrySet()) {
      final String modelID = modelRawDataEntry.getKey();
      final RiskModelRawData modelRawData = modelRawDataEntry.getValue();
      // Convert ageYrs int array to double array.
      final int[] ageYrsIntArr = modelRawData.getAgeYrs();
      final double[] ageYrsArr = new double[ageYrsIntArr.length];
      for (int ageIdx = 0; ageIdx < ageYrsIntArr.length; ageIdx++) {
        ageYrsArr[ageIdx] = ageYrsIntArr[ageIdx];
      }
      // Construct RiskModel object for current modelID.
      final RiskModel modelRiskModel =
        new RiskModel(modelRawData.getModelName(), modelRawData.getModelSNPs(),
          ageYrsArr, modelRawData.getMargSurv());
      // Create path to [modelID].rmo file to store serialized RiskModel object.
      final Path modelOutputPath =
        Paths.get(savePath).resolve(modelID + ".rmo");
      // Try to open an ObjectOutputStream that writes to this file.
      try (
        ObjectOutputStream rmObjOut =
          new ObjectOutputStream(Files.newOutputStream(modelOutputPath))) {
        // Write current modelID's RiskModel object to [modelID].rmo file.
        rmObjOut.writeObject(modelRiskModel);
      } // ObjectOutputStream automatically closed here.
    }
  }

  /**
   * Main method of {@code RiskModelBuilder} class.
   * 
   * @param args {@code String[]} providing required command line arguments.
   *          {@code args[0]} should contain the path to the model SNPs input
   *          file, {@code args[1]} should contain the path to the annual
   *          incidence input file, and {@code args[2]} should contain the path
   *          to the directory where the serialized {@code RiskModel} objects
   *          (*.rmo files) should be saved.
   */
  public static void main(String[] args) {
    // Set exit code to 0 (success).
    int exitCode = 0;
    // Parse arguments.
    final String modelSNPsFile = args[0];
    final String annIncFile = args[1];
    final String savePath = args[2];
    // Try parsing input files and building RiskModel objects.
    try {
      HashMap<String, RiskModelRawData> modelRawDataMap =
        parseInputFiles(modelSNPsFile, annIncFile);
      buildRiskModels(savePath, modelRawDataMap);
    } catch (Exception e) {
      /*
       * If any exception is thrown during the process, a message will be
       * printed to standard output and details provided. The exit code will be
       * set to 1 (failure) to alert the build system to the failure.
       */
      System.out.println("Problem building or writing risk model object "
        + "files.");
      System.out.println(e.toString());
      exitCode = 1;
    }
    // Exit with the appropriate exit code.
    System.exit(exitCode);
  }

}