/*
 * This source file is part of the PRISM software package.
 * 
 * Copyright 2014 The Ohio State University Wexner Medical Center
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
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionGroup;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.PosixParser;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.Precision;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * Constructs and serializes {@link RiskModel} objects during the build process.
 * <p>
 * This package-private class contains a private nested helper class, several private static
 * methods, and a {@code main} method called during the build process. Because it is used only in
 * the build process, it is neither part of the package API nor included in the final artifact.
 * 
 * @author Daniel Kinnamon
 * @version 2014-10-31
 * @since 2014-10-30
 */
final class RiskModelBuilder {

  // Initialize logger
  private static final Logger LOGGER = LogManager.getLogger();

  /**
   * Prints a command-line help message to a {@code PrintWriter}.
   * 
   * @param options {@link Options} object containing possible options for the command line
   * @param pw {@code PrintWriter} for output
   */
  private static final void printHelp(final Options options, final PrintWriter pw) {
    final String version = Package.getPackage("edu.osumc.prism").getImplementationVersion();
    new HelpFormatter().printHelp(pw, HelpFormatter.DEFAULT_WIDTH, "java -jar prism-" + version
        + ".jar", "", options, HelpFormatter.DEFAULT_LEFT_PAD, HelpFormatter.DEFAULT_DESC_PAD, "",
        true);
  }

  /**
   * Prints the arguments with which the program was invoked to a {@code PrintWriter}.
   * 
   * @param cmd {@link CommandLine} object containing parsed command line
   * @param pw {@code PrintWriter} for output
   */
  private static final void printArgs(final CommandLine cmd, final PrintWriter pw) {
    pw.println();
    pw.println("Invoked with arguments:");
    for (Option option : cmd.getOptions()) {
      pw.println(HelpFormatter.DEFAULT_LONG_OPT_PREFIX + option.getLongOpt()
          + (option.hasArg() ? " " + option.getValue() : ""));
    }
    pw.println();
  }


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
      if (ageYrsList.size()-1 < ageYrs || ageYrsList.get(ageYrs) != ageYrs || !Precision.equals(annIncList.get(ageYrs), annInc)) {
        throw new IllegalArgumentException("Annual incidences for "
            + "consecutive yearly ages starting at 0 must be supplied in "
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
      ArrayList<Double> margSurvList = new ArrayList<Double>();
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
          throw new RuntimeException("Problem with calculating survivor "
              + "function from annual incidences.");
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
        public boolean accept(Path file) throws IOException {
          return (file.getFileName().toString().endsWith("_SNPs.dat"));
        }
      };

  /**
   * Constructs a {@link DirectoryStream.Filter<Path>} object for use in creating a new directory
   * stream containing only annInc files.
   */
  private static final DirectoryStream.Filter<Path> filterannInc =
      new DirectoryStream.Filter<Path>() {
        public boolean accept(Path file) throws IOException {
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
  private static final HashMap<String, RiskModelRawData> parseInputFiles(
      final String targetModelID, final String sourceDirectory) throws IOException {
    // Try to run the application.
    HashMap<String, RiskModelRawData> modelRawDataMap = new HashMap<String, RiskModelRawData>();
    // Parse model SNPs file(s).

    // Start by attempting to read all SNPs.dat files in the source directory
    LOGGER.debug("Now attempting to load SNPs.dat files from directory: " + sourceDirectory);
    Path sourceDirPath = Paths.get(sourceDirectory);
    LOGGER.debug("Path resolves to: " + sourceDirPath.toString());
    try (DirectoryStream<Path> sourceFolder =
        Files.newDirectoryStream(sourceDirPath, filterSNPs)) {
      LOGGER.debug("Directory stream for SNPs.dat files open..");
      
      // For each file in the directory, compare to the targetModelID. If it matches, process it
      for (Path snpFile : sourceFolder) {
        if (snpFile.endsWith(targetModelID + "_SNPs.dat") || targetModelID.isEmpty()) {
          LOGGER.debug("Now attempting to read SNPs file: " + snpFile.toString());
          try (BufferedReader fileReader = Files.newBufferedReader(snpFile, StandardCharsets.UTF_8)) {
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
                  List<String> expColHeaderList =
                      Arrays.asList("modelID", "rsID", "sourcePub", "allele1", "allele2",
                          "orientRs", "allele2Freq", "allele2lnHR");
                  for (String expColHeader : expColHeaderList) {
                    if (!expColHeader.equals(lineScanner.next())) {
                      throw new IOException("Model SNPs file column headers " + "are incorrect.");
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
                  modelRawDataMap.get(modelID).addSNP(
                      new SNP(rsID, sourcePub, allele1, allele2,
                          orientRsStr.matches("FORWARD") ? SNP.AlleleOrientation.FORWARD
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
          }
        } // Model SNPs file BufferedReader automatically closed here.
      }
    } catch (Exception e){
      LOGGER.error("Exception thrown. Details:", e);
    }

    // Parse annual Incidence file(s).
    // Start by attempting to read all annInc.dat files in the source directory
    LOGGER.debug("Now attempting to load annInc.dat files...");
    
    
    // Start by attempting to read all SNPs.dat files in the source directory
    try (DirectoryStream<Path> sourceFolder =
        Files.newDirectoryStream(sourceDirPath, filterannInc)) {
      // For each file in the directory, compare to the targetModelID. If it matches, process it
      for (Path annIncFile : sourceFolder) {
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
                  List<String> expColHeaderList = Arrays.asList("modelID", "ageYrs", "annInc");
                  for (String expColHeader : expColHeaderList) {
                    if (!expColHeader.equals(lineScanner.next())) {
                      throw new IOException("Annual incidence file column headers "
                          + "are incorrect.");
                    }
                  }
                } else {
                  // Otherwise, read annual incidence data.
                  final String modelID = lineScanner.next();
                  int ageYrs = lineScanner.nextInt();
                  double annInc = lineScanner.nextDouble();

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
    for (Map.Entry<String, RiskModelRawData> modelRawDataEntry : modelRawDataMap.entrySet()) {
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
          new RiskModel(modelRawData.getModelName(), modelRawData.getModelSNPs(), ageYrsArr,
              modelRawData.getMargSurv());
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
   * @param args {@code String[]} providing required command line arguments. {@code args[0]} should
   *        contain the path to the model SNPs input file, {@code args[1]} should contain the path
   *        to the annual incidence input file, and {@code args[2]} should contain the path to the
   *        directory where the serialized {@code RiskModel} objects (*.rmo files) should be saved.
   */
  public static void main(String[] args) {
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
    final Option help = new Option("h", "help", false, "Print this usage information.");
    final Option build =
        new Option("b", "buildModel", false,
            "Generate new risk predictions model using the supplied annual incidences "
                + "and risk allele hazard ratios in given by <modelID>_annInc.dat / "
                + "<modelID>_SNPs.dat and save the resulting model as a risk model object"
                + "(<modelID>.rmo).");
    
    final Option optModelID = new Option("m", "modelID", true, "Set the modelID that the "
        + "builder will attempt to build. If left blank, the builder will attempt to build"
        + "all models in the source directory (<SourceDir>).");
    optModelID.setArgName("modelID");

    final Option sources =
        new Option("s", "SourceDirectory", true, "Set the directory in "
            + "which the builder will look for model source files (<SourceDir>). If "
            + "left unset, the builder will look in the directory containing itself.");
    sources.setArgName("sourceDir");

    final Option destination =
        new Option("d", "DestinationDirectory", true, "Set the"
            + " directory in which the output model will be saved. If left "
            + "unset, the builder will save the model in the directory " + "containing itself.");
    destination.setArgName("saveDir");


    /*
     * Make options a mutually exclusive group so only one can be selected at a time. Make the group
     * required so that at least one of these options must be selected.
     */
    final OptionGroup modes = new OptionGroup();
    modes.addOption(help);
    modes.addOption(build);
    modes.setRequired(true);
    // Set command-line options.
    final Options options = new Options();
    options.addOptionGroup(modes);
    options.addOption(optModelID);
    options.addOption(sources);
    options.addOption(destination);

    LOGGER.debug("Building of options complete. Now attempting to parse them.");
    try {
      // Parse command line and echo arguments to standard output.
      final CommandLine cmd = new PosixParser().parse(options, args);
      printArgs(cmd, stdOutPw);

      switch (modes.getSelected()) {
        case "h":
          // If in "h" mode, print usage message to standard output.
          printHelp(options, stdOutPw);
          break;
        case "b":
          LOGGER.debug("Build mode selected");
          /*
           * If in "b" mode, build specified risk models, and produce output files containing the
           * saved models to the specified directory
           */
          String modelID = cmd.getOptionValue("modelID", "");
          String sourceDir = cmd.getOptionValue("sourceDir", "");
          String outputDir = cmd.getOptionValue("destinationDir", "");
          LOGGER.debug("Confirming, read modelID: " + modelID);
          LOGGER.debug("Confirming, read sourceDir: " + sourceDir);
          LOGGER.debug("Confirming, read outputDir: " + outputDir);
          
          // Try parsing input files and building RiskModel objects.

          LOGGER.debug("Beginning attempt to parse input files..");
          HashMap<String, RiskModelRawData> modelRawDataMap = parseInputFiles(modelID, sourceDir);
          
          LOGGER.debug("Parsing complete. Beginning attempt to build risk models...");
          LOGGER.debug(modelRawDataMap.toString());
          buildRiskModels(outputDir, modelRawDataMap);

          break;
      }



    } catch (Exception e) {
      /*
       * If any exception is thrown during the process, a message will be printed to standard output
       * and details provided. The exit code will be set to 1 (failure) to alert the build system to
       * the failure.
       */
      System.out.println("Problem building or writing risk model object " + "files.");
      LOGGER.error("Exception thrown. Details:", e);
      System.out.println(e.toString());
      exitCode = 1;
    }
    // Exit with the appropriate exit code.
    System.exit(exitCode);
  }

}
