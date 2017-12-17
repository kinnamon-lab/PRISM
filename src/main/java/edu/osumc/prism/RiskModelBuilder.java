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

import java.io.BufferedReader;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.NoSuchFileException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Scanner;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.text.WordUtils;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.Precision;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * Constructs and serializes {@link RiskModel} objects.
 * <p>
 * This package-private class contains a private nested helper class, several private static
 * methods, and a {@code main} method for building serialized risk model objects.
 */
final class RiskModelBuilder {

  // Initialize logger
  private static final Logger LOGGER = LogManager.getLogger();

  /** Holds raw data obtained from input files. */
  private static final class RiskModelRawData {

    /** {@code String} risk model name. */
    private final String modelID;
    /** {@code ArrayList<SNP>} list of model SNPs. */
    private final ArrayList<SNP> modelSNPsList;
    /**
     * {@code ArrayList<Integer>} list of of consecutive annual ages.
     */
    private final ArrayList<Integer> ageYrsList;
    /**
     * {@code ArrayList<Double>} list of annual incidences, where {@code annIncList[i]} is the
     * annual incidence in the year preceding age {@code ageYrsList[i]}.
     */
    private final ArrayList<Double> annIncList;

    /**
     * Constructs a {@code RiskModelRawData} object ready to receive data.
     *
     * @param modelID {@code String} risk model name
     */
    private RiskModelRawData(final String modelID) {
      this.modelID = modelID;
      modelSNPsList = new ArrayList<SNP>();
      ageYrsList = new ArrayList<Integer>();
      annIncList = new ArrayList<Double>();
    }

    /**
     * Adds SNP to the end of the {@code modelSNPList} object member. SNPs can be added in any order
     * for a given risk model.
     */
    private final void addSNP(final SNP newSNP) {
      modelSNPsList.add(newSNP);
    }

    /**
     * Adds a consecutive yearly age to the end of the {@code ageYrsList} object member and the
     * annual incidence in the year preceding this age to the end of the {@code annIncList} object
     * member. Performs checks to make sure that annual incidences for consecutive yearly ages are
     * being added in order and that {@code ageYrsList} and {@code annIncList} are growing at the
     * same rate.
     *
     * @param ageYrs {@code int} age in years
     * @param annInc {@code double} annual incidence in the year preceding {@code ageYrs}
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
       * Check to make sure that ageYrsList[ageYrs] = ageYrs and annIncList[ageYrs] = annInc.
       */
      if (ageYrsList.size() - 1 < ageYrs || ageYrsList.get(ageYrs) != ageYrs
          || !Precision.equals(annIncList.get(ageYrs), annInc)) {
        throw new IllegalArgumentException(
            "Annual incidences for " + "consecutive yearly ages starting at 0 must be supplied in "
                + "order. Please check the annual incidence input file.");
      }
      // Check that lists are growing at the same rate.
      if (ageYrsList.size() != annIncList.size()) {
        throw new IllegalArgumentException("ageYrsList and annIncList are "
            + "growing at different rates. Please check the annual incidence " + " input file.");
      }
    }

    /** Returns {@code String} risk model name. */
    private final String getModelName() {
      return modelID;
    }

    /** Returns {@code SNP[]} array of model SNPs. */
    private final SNP[] getModelSNPs() {
      return modelSNPsList.toArray(new SNP[0]);
    }

    /**
     * Returns {@code int[]} array of of ages at which the the marginal survivor function is
     * evaluated (abscissae).
     */
    private final int[] getAgeYrs() {
      return ArrayUtils.toPrimitive(ageYrsList.toArray(new Integer[0]));
    }

    /**
     * Returns {@code double[]} array of marginal survivor function values, where
     * {@code getMargSurv()[i]} is the marginal survivor function at {@code getAgeYrs()[i]}.
     */
    private final double[] getMargSurv() {
      final ArrayList<Double> margSurvList = new ArrayList<Double>();
      /*
       * Use the annual incidences to calculate the continuous-time survivor function assuming a
       * piecewise constant yearly hazard (i.e., annual incidence). Because of checks in the
       * addAnnInc method, annIncList[0] = 0 and the indices in ageYrsList and annIncList both
       * correspond to annual ages. The cumulative hazard from age 0 to age i, H(i), is therefore
       * given by sum_{j=0}^{i} annIncList[j], or the cumulative sum over elements 0 to i of
       * annIncList. Calculating the cumulative sum at each index i therefore gives the cumulative
       * hazard to age i, and the relationship S(i) = exp(-H(i)) is used to obtain the survivor
       * function at age i, which is stored in margSurvList[i].
       */
      double cumHazI = 0.0D;
      for (int i = 0; i < annIncList.size(); i++) {
        cumHazI += annIncList.get(i);
        final double survI = FastMath.exp(-cumHazI);
        margSurvList.add(survI);
        /*
         * Check that S(i) was added to margSurvList[i] so that S(ageYrsList[i]) = S(i) =
         * margSurvList[i].
         */
        if (!Precision.equals(margSurvList.get(i), survI)) {
          throw new RuntimeException(
              "Problem with calculating survivor function from annual incidences.");
        }
      }
      return ArrayUtils.toPrimitive(margSurvList.toArray(new Double[0]));
    }
  }

  /** Prevents construction of {@code RiskModelBuilder} class instances. */
  private RiskModelBuilder() {}

  /**
   * Constructs a {@link DirectoryStream.Filter<Path>} object for use in creating a new directory
   * stream containing only SNP files.
   */
  private static final DirectoryStream.Filter<Path> filterSNPs =
      new DirectoryStream.Filter<Path>() {
        @Override
        public boolean accept(final Path file) throws IOException {
          return (file.getFileName().toString().endsWith("_SNPs.dat"));
        }
      };

  /**
   * Constructs a {@link DirectoryStream.Filter<Path>} object for use in creating a new directory
   * stream containing only annInc files.
   */
  private static final DirectoryStream.Filter<Path> filterannInc =
      new DirectoryStream.Filter<Path>() {
        @Override
        public boolean accept(final Path file) throws IOException {
          return (file.getFileName().toString().endsWith("_annInc.dat"));
        }
      };

  /**
   * Parses model SNPs and annual incidence files and returns a {@link HashMap} linking
   * {@code String} model ID keys to {@code RiskModelRawData} values that can be used to build a
   * {@link RiskModel} object.
   *
   * @param targetModelID {@code String} Specific modelID to build, or left blank to build all
   *        models present in the source directory
   * @param sourceDirectory {@code String} path to directory containing all properly formatted SNP
   *        and Annual Incidence files, or the files belonging to the specified targetModelID
   */
  private static final HashMap<String, RiskModelRawData> parseInputFiles(final String targetModelID,
      final String sourceDirectory) throws IOException {
    // Try to run the application.
    final HashMap<String, RiskModelRawData> modelRawDataMap =
        new HashMap<String, RiskModelRawData>();
    // Parse model SNPs file(s).

    // Start by attempting to read all SNPs.dat files in the source directory
    LOGGER.debug("Now attempting to load SNPs.dat files from directory: " + sourceDirectory);
    final Path sourceDirPath = Paths.get(sourceDirectory);
    LOGGER.debug("Path resolves to: " + sourceDirPath.toString());
    try (DirectoryStream<Path> sourceFolder = Files.newDirectoryStream(sourceDirPath, filterSNPs)) {
      LOGGER.debug("Directory stream for SNPs.dat files open..");

      // For each file in the directory, compare to the targetModelID. If it matches, process it
      for (final Path snpFile : sourceFolder) {
        if (snpFile.endsWith(targetModelID + "_SNPs.dat") || targetModelID.isEmpty()) {
          LOGGER.debug("Now attempting to read SNPs file: " + snpFile.toString());
          try (BufferedReader fileReader =
              Files.newBufferedReader(snpFile, StandardCharsets.UTF_8)) {
            String line;
            int linesRead = 0;
            // Read lines from file...
            while ((line = fileReader.readLine()) != null) {
              // For each line, try to create a scanner using tab delimiter.
              try (Scanner lineScanner =
                  new Scanner(line).useLocale(Locale.US).useDelimiter("\t").useRadix(10)) {
                // Parse tab-delimited tokens from the current line.
                if (linesRead == 0) {
                  // If first line, check that file column headers are correct.
                  final List<String> expColHeaderList = Arrays.asList("modelID", "rsID",
                      "sourcePub", "allele1", "allele2", "orientRs", "allele2Freq", "allele2lnHR");
                  for (final String expColHeader : expColHeaderList) {
                    if (!expColHeader.equals(lineScanner.next())) {
                      throw new IOException("Model SNPs file column headers are incorrect.");
                    }
                  }
                } else {
                  // Otherwise, read SNP data.
                  final String modelID = lineScanner.next();
                  final String rsID = lineScanner.next();
                  final String sourcePub = lineScanner.next();
                  final String allele1 = lineScanner.next();
                  final String allele2 = lineScanner.next();
                  final String orientRsStr = lineScanner.next().toUpperCase(Locale.US);
                  final double allele2Freq = lineScanner.nextDouble();
                  final double allele2lnHR = lineScanner.nextDouble();

                  if (!orientRsStr.matches("^(FORWARD|REVERSE)$")) {
                    throw new IOException("Model SNPs file allele orientation "
                        + "relative to RefSNP must be 'Forward' or 'Reverse'.");
                  }
                  // If new modelID, construct map entry.
                  if (!modelRawDataMap.containsKey(modelID)) {
                    modelRawDataMap.put(modelID, new RiskModelRawData(modelID));
                  }
                  // Add SNP to RiskModelRawData object for modelID.
                  modelRawDataMap.get(modelID)
                      .addSNP(new SNP(rsID, sourcePub,
                          allele1, allele2, orientRsStr.matches("FORWARD")
                              ? SNP.AlleleOrientation.FORWARD : SNP.AlleleOrientation.REVERSE,
                          allele2Freq, allele2lnHR));
                }
                // Check that there are no input tokens other than the ones expected.
                if (lineScanner.hasNext()) {
                  throw new IOException(
                      "Input tokens available past expected end of line in model SNPs file.");
                }
              } // Scanner automatically closed here.
              // Increment number of model SNPs file lines read.
              linesRead++;
            }
          }
        } // Model SNPs file BufferedReader automatically closed here.
      }
    } catch (final Exception e) {
      LOGGER.error("Exception thrown. Details:", e);
    }

    // Parse annual Incidence file(s).
    // Start by attempting to read all annInc.dat files in the source directory
    LOGGER.debug("Now attempting to load annInc.dat files...");


    // Start by attempting to read all SNPs.dat files in the source directory
    try (DirectoryStream<Path> sourceFolder =
        Files.newDirectoryStream(sourceDirPath, filterannInc)) {
      // For each file in the directory, compare to the targetModelID. If it matches, process it
      for (final Path annIncFile : sourceFolder) {
        if (annIncFile.endsWith(targetModelID + "_annInc.dat") || targetModelID.isEmpty()) {
          LOGGER.debug("Now attempting to read annInc file: " + annIncFile.toString());
          try (BufferedReader fileReader =
              Files.newBufferedReader(annIncFile, StandardCharsets.UTF_8)) {
            String line;
            int linesRead = 0;
            // Read lines from file...
            while ((line = fileReader.readLine()) != null) {
              // For each line, try to create a scanner using tab delimiter.
              try (Scanner lineScanner =
                  new Scanner(line).useLocale(Locale.US).useDelimiter("\t").useRadix(10)) {
                // Parse tab-delimited tokens from the current line.
                if (linesRead == 0) {
                  // If first line, check that column headers are correct.
                  final List<String> expColHeaderList =
                      Arrays.asList("modelID", "ageYrs", "annInc");
                  for (final String expColHeader : expColHeaderList) {
                    if (!expColHeader.equals(lineScanner.next())) {
                      throw new IOException("Annual incidence file column headers are incorrect.");
                    }
                  }
                } else {
                  // Otherwise, read annual incidence data.
                  final String modelID = lineScanner.next();
                  final int ageYrs = lineScanner.nextInt();
                  final double annInc = lineScanner.nextDouble();

                  /*
                   * If the model ID for the current line corresponds to a model for which we have
                   * SNPs, then process the current line.
                   */
                  if (modelRawDataMap.containsKey(modelID)) {
                    /*
                     * Annual incidence file rows should be in order of consecutive annual ages
                     * within each modelID. This is checked by the addAnnInc method.
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
        }
      } // End of for loop for each annInc file
    }
    // After files are processed, return modelRawDataMap.
    return modelRawDataMap;
  }


  /**
   * Constructs a {@link RiskModel} object for each model with a {@code RiskModelRawData} entry in
   * the input {@link HashMap} and saves serialized version (*.rmo file) to the specified path.
   *
   * @param savePath {@code String} path in which to save *.rmo files
   * @param modelRawDataMap {@code HashMap<String,
   *          RiskModelRawData>} with each entry containing a {@code String} model ID key and a
   *        {@code RiskModelRawData} value containing input data parsed from input files
   */
  private static final void buildRiskModels(final String savePath,
      final HashMap<String, RiskModelRawData> modelRawDataMap) throws IOException {
    // For each modelRawDataMap entry...
    for (final Map.Entry<String, RiskModelRawData> modelRawDataEntry : modelRawDataMap.entrySet()) {
      final String modelID = modelRawDataEntry.getKey();
      final RiskModelRawData modelRawData = modelRawDataEntry.getValue();
      // Convert ageYrs int array to double array.
      final int[] ageYrsIntArr = modelRawData.getAgeYrs();
      final double[] ageYrsArr = new double[ageYrsIntArr.length];
      for (int ageIdx = 0; ageIdx < ageYrsIntArr.length; ageIdx++) {
        ageYrsArr[ageIdx] = ageYrsIntArr[ageIdx];
      }
      // Construct RiskModel object for current modelID.
      final RiskModel modelRiskModel = new RiskModel(modelRawData.getModelName(),
          modelRawData.getModelSNPs(), ageYrsArr, modelRawData.getMargSurv());
      // Create path to [modelID].rmo file to store serialized RiskModel object.
      final Path modelOutputPath = Paths.get(savePath).resolve(modelID + ".rmo");
      // Try to open an ObjectOutputStream that writes to this file.
      try (ObjectOutputStream rmObjOut =
          new ObjectOutputStream(Files.newOutputStream(modelOutputPath))) {
        // Write current modelID's RiskModel object to [modelID].rmo file.
        rmObjOut.writeObject(modelRiskModel);
      } // ObjectOutputStream automatically closed here.
    }
  }

  /**
   * Main method of {@code RiskModelBuilder} class.
   *
   * @param args {@code String[]} providing required command line arguments
   */
  public static void main(final String[] args) {
    // Set exit code to 0 (success).
    int exitCode = 0;

    /*
     * Wrap standard output in PrintWriter (necessary for some Apache Commons CLI methods).
     */
    final PrintWriter stdOutPw =
        new PrintWriter(new OutputStreamWriter(System.out, StandardCharsets.UTF_8), true);

    LOGGER.debug("Welcome to RiskModelBuilder.");
    LOGGER.debug("Now attempting to build options.");

    // Define command line options.
    final Option help = new Option("h", "help", false, "Print this usage information and exit.");

    final Option optModelID = new Option("m", "model_id", true,
        "The builder will build a risk model using the supplied annual incidences "
            + "and risk allele hazard ratios in given by <modelID>_annInc.dat / "
            + "<modelID>_SNPs.dat and save the resulting model as a risk model object "
            + "<modelID>.rmo. If left blank, the builder will attempt to build models "
            + "for all source files in the source directory <sourceDir> and parse "
            + "the modelIDs from these file names.");
    optModelID.setArgName("modelID");

    final Option sources = new Option("s", "src_dir", true,
        "The directory in which the builder will look for model source files. If "
            + "left unset, the builder will look in the current working directory.");
    sources.setArgName("sourceDir");

    final Option destination = new Option("o", "out_dir", true,
        "The directory in which the risk model object files will be saved. If left "
            + "unset, the builder will save the risk model object files in the current "
            + "working directory.");
    destination.setArgName("outputDir");

    // Set command-line options.
    final Options options = new Options();
    options.addOption(help);
    options.addOption(optModelID);
    options.addOption(sources);
    options.addOption(destination);

    LOGGER.debug("Building of options complete. Now attempting to parse them.");
    try {
      // Print masthead to standard output.
      Util.printMasthead(RiskModelBuilder.class.getSimpleName(), stdOutPw);
      // Parse command line and echo arguments to standard output.
      final CommandLine cmd = new PosixParser().parse(options, args);
      Util.printArgs(cmd, stdOutPw);

      if (cmd.hasOption("h")) {
        // If in "h" mode, print usage message to standard output.
        Util.printHelp(options, RiskModelBuilder.class.getSimpleName(), stdOutPw);
      } else {
        LOGGER.debug("Build mode selected");
        /*
         * Build specified risk models, and produce output files containing the saved models to the
         * specified directory
         */
        final String modelID = cmd.getOptionValue("model_id", "");
        final String sourceDir = cmd.getOptionValue("src_dir", System.getProperty("user.dir"));
        final String outputDir = cmd.getOptionValue("out_dir", System.getProperty("user.dir"));
        LOGGER.debug("Confirming, read modelID: " + modelID);
        LOGGER.debug("Confirming, read sourceDir: " + sourceDir);
        LOGGER.debug("Confirming, read outputDir: " + outputDir);

        // Try parsing input files and building RiskModel objects.

        LOGGER.debug("Beginning attempt to parse input files..");
        final HashMap<String, RiskModelRawData> modelRawDataMap =
            parseInputFiles(modelID, sourceDir);

        LOGGER.debug("Parsing complete. Beginning attempt to build risk models...");
        LOGGER.debug(modelRawDataMap.toString());
        buildRiskModels(outputDir, modelRawDataMap);
      }
    } catch (final ParseException e) {
      LOGGER.error("Parse exception!");
      /*
       * If a ParseException is thrown to indicate a problem with the command line, print the
       * ParseException message and the usage message and set the exit code to 1 (failure).
       */
      stdOutPw.println();
      stdOutPw.println(WordUtils.wrap(e.getMessage(), 80));
      stdOutPw.println();
      Util.printHelp(options, RiskModelBuilder.class.getSimpleName(), stdOutPw);
      exitCode = 1;
    } catch (final NoSuchFileException e) {
      stdOutPw.println("ERROR! File not found!");
      stdOutPw.println(
          "The following file could not be found. Please check your command and try again.");
      stdOutPw.println(WordUtils.wrap(e.getMessage(), 80));
      exitCode = 1;
    } catch (final Exception e) {
      /*
       * If any other exception is thrown during the process, a message will be printed to standard
       * output and details provided. The exit code will be set to 1 (failure).
       */
      LOGGER.error("Exception thrown. Details:", e);
      System.out.println(e.toString());
      exitCode = 1;
    }
    // Flush and close the standard output PrintWriter.
    stdOutPw.flush();
    stdOutPw.close();
    // Exit with the appropriate exit code.
    System.exit(exitCode);
  }

}
