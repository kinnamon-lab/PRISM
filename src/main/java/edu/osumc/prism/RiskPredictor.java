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
import java.io.InputStream;
import java.io.ObjectInputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.net.JarURLConnection;
import java.net.URISyntaxException;
import java.nio.charset.StandardCharsets;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.NoSuchFileException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Locale;
import java.util.Map;
import java.util.Scanner;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionGroup;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import org.apache.commons.lang3.text.WordUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * Command-line interface for obtaining risk predictions from the PRISM package.
 * <p>
 * This package-private class contains several private static methods and a {@code main} method to
 * provide an application entry point in the package JAR file.
 */
final class RiskPredictor {

  // Initialize logger
  private static final Logger LOGGER = LogManager.getLogger();

  /** Prevents construction of {@code RiskPredictor} class instances. */
  private RiskPredictor() {}

  /**
   * Constructs a {@link DirectoryStream.Filter<Path>} object for use in creating a new directory
   * stream containing only SNP files.
   */
  private static final DirectoryStream.Filter<Path> filterRMOs =
      new DirectoryStream.Filter<Path>() {
        @Override
        public boolean accept(final Path file) throws IOException {
          return (file.getFileName().toString().endsWith(".rmo"));
        }
      };

  /**
   * Searches the {@link RiskModel} objects (*.rmo files) contained in the modelDirectory for
   * targetModelID and returns a {@code LinkedHashMap} of the model identified by its modelID.
   */
  private static final LinkedHashMap<String, RiskModel> loadRiskModel(final String modelDirectory,
      final String targetModelID) throws URISyntaxException, IOException, ClassNotFoundException {
    final LinkedHashMap<String, RiskModel> riskModelMap = new LinkedHashMap<String, RiskModel>();

    LOGGER.debug("Attempting to load file at: " + modelDirectory + targetModelID + ".rmo");
    final Path modelPath = Paths.get(modelDirectory, targetModelID + ".rmo");
    try (InputStream modelStream = Files.newInputStream(modelPath)) {
      try (ObjectInputStream riskModelObjInStream = new ObjectInputStream(modelStream)) {
        /*
         * Deserialize RiskModel object entry to RiskModel object stored in map with key modelID.
         */
        riskModelMap.put(targetModelID, (RiskModel) riskModelObjInStream.readObject());
        LOGGER.debug(targetModelID + " model successfully loaded");
      } // ObjectInputStream automatically closed here.
    }
    if (riskModelMap.isEmpty()) {
      throw new NoSuchFileException(
          targetModelID + ".rmo not found in " + Paths.get(modelDirectory).toString());
    }
    return riskModelMap;
  }

  /**
   * Loads all {@link RiskModel} objects (*.rmo files) contained in the modelDirectory and returns a
   * {@code LinkedHashMap} of these objects identified by filename root/modelID keys.
   */
  private static final LinkedHashMap<String, RiskModel> loadRiskModels(final String modelDirectory)
      throws URISyntaxException, IOException, ClassNotFoundException {
    final LinkedHashMap<String, RiskModel> riskModelMap = new LinkedHashMap<String, RiskModel>();
    // Build list of paths to Risk Model Objects and then load them into
    // the riskModelMap
    try (DirectoryStream<Path> modelsFolder =
        Files.newDirectoryStream(Paths.get(modelDirectory), filterRMOs)) {
      LOGGER.debug("Directory stream for .rmo files open.");

      for (final Path modelPath : modelsFolder) {
        LOGGER.debug("Now looking at " + modelPath.getFileName().toString());
        // Slightly redundant due to filter, but better safe than sorry
        if (modelPath.getFileName().toString().endsWith(".rmo")) {
          // If RiskModel object entry, extract modelID from filename.
          final String modelID = modelPath.getFileName().toString().replaceFirst("\\.rmo$", "");
          // Create new ObjectInputStream to read RiskModel object entry.
          try (ObjectInputStream riskModelObjInStream =
              new ObjectInputStream(Files.newInputStream(modelPath))) {
            /*
             * Deserialize RiskModel object entry to RiskModel object stored in map with key
             * modelID.
             */
            riskModelMap.put(modelID, (RiskModel) riskModelObjInStream.readObject());
            LOGGER.debug(modelID + " successfully loaded.");
          } // ObjectInputStream automatically closed here.
        }
      }
    } // DirectoryStream automatically closed here.

    if (riskModelMap.isEmpty()) {
      throw new NoSuchFileException(
          "No .rmo files found in " + Paths.get(modelDirectory).toString());
    }
    return riskModelMap;
  }

  /**
   * Prints loaded {@link RiskModel} information to a {@code PrintWriter}.
   *
   * @param riskModelMap {@code LinkedHashMap<String, RiskModel>} containing modelID-RiskModel
   *        key-value pairs
   * @param pw {@code PrintWriter} for output
   */
  private static final void printRiskModels(final LinkedHashMap<String, RiskModel> riskModelMap,
      final PrintWriter pw) {
    pw.println();
    for (final RiskModel riskModel : riskModelMap.values()) {
      riskModel.print(pw);
      pw.println();
    }
  }

  /**
   * Reads individual genotypes from a PED/MAP file pair and outputs individual risk predictions
   * from each loaded {@link RiskModel} object to a separate file, logging the process.
   *
   * @param riskModelMap {@code LinkedHashMap<String, RiskModel>} containing modelID-RiskModel
   *        key-value pairs
   * @param genoFilesRoot {@code String} root path to PED/MAP file pair (i.e., without file
   *        extensions)
   * @param pw {@code PrintWriter} for logging
   */
  private static final void outputRiskPredictions(
      final LinkedHashMap<String, RiskModel> riskModelMap, final String genoFilesRoot,
      final PrintWriter pw) throws IOException {
    // Parse MAP File.
    final LinkedHashMap<String, String> orientRsMap = new LinkedHashMap<String, String>();
    final Path mapFilePath = Paths.get(genoFilesRoot + ".map");
    try (BufferedReader mapBufferedReader =
        Files.newBufferedReader(mapFilePath, StandardCharsets.UTF_8)) {
      String line;
      // Read lines from file...
      while ((line = mapBufferedReader.readLine()) != null) {
        /*
         * For each line, try to create a scanner using the default whitespace delimiter.
         */
        try (Scanner lineScanner = new Scanner(line).useLocale(Locale.US)) {
          // String rs number.
          final String rsID = lineScanner.next();
          // String allele orientation.
          final String orientRsStr = lineScanner.next();
          // Add SNP to orientRsMap (entries maintained in order of addition).
          orientRsMap.put(rsID, orientRsStr);
          // Check that there are no input tokens other than the ones expected.
          if (lineScanner.hasNext()) {
            throw new IOException("Input tokens available past expected end of line in MAP file.");
          }
        } // Scanner automatically closed here.
      }
    } // MAP file BufferedReader automatically closed here.

    // Parse PED file.
    final ArrayList<Individual.Genotypes> inputGenoList = new ArrayList<Individual.Genotypes>();
    final Path pedFilePath = Paths.get(genoFilesRoot + ".ped");
    try (BufferedReader pedBufferedReader =
        Files.newBufferedReader(pedFilePath, StandardCharsets.UTF_8)) {
      String line;
      // Read lines from file...
      while ((line = pedBufferedReader.readLine()) != null) {
        /*
         * For each line, try to create a scanner using the default whitespace delimiter.
         */
        try (Scanner lineScanner = new Scanner(line).useLocale(Locale.US)) {
          // Individual ID.
          final String indivID = lineScanner.next();
          // Construct Genotypes object for this individual.
          final Individual.Genotypes inputGenos = new Individual.Genotypes(indivID);
          /*
           * Load the individual's genotypes for each SNP into the Genotypes object. Pairs of SNP
           * allele columns appear in the same order from left to right as the rsIDs appeared in the
           * rows of the MAP file. Because orientRsMap is a LinkedHashMap and rsID keys were added
           * in the order that they appeared in the rows of the MAP file, we can iterate over the
           * rsID keys in orientRsMap and pairs of allele columns in the PED file to get the
           * individual's genotype for each rsID key.
           */
          for (final String rsID : orientRsMap.keySet()) {
            /*
             * Get data from next pair of allele columns corresponding to current rsID.
             */
            final String allele1 = lineScanner.next();
            final String allele2 = lineScanner.next();
            // Add individual's data for current SNP to Genotypes object.
            inputGenos.addGenotype(rsID, allele1, allele2, orientRsMap.get(rsID));
          }
          // Check that there are no input tokens other than the ones expected.
          if (lineScanner.hasNext()) {
            throw new IOException("Input tokens available past expected end of line in PED file.");
          }
          // Add individual's Genotype object to the end of the list.
          inputGenoList.add(inputGenos);
        } // Scanner automatically closed here.
      }
    } // PED file BufferedReader automatically closed here.

    /*
     * Output risk predictions for all individuals from each model to separate files.
     */

    pw.println("Outputting risk predictions to: ");
    for (final Map.Entry<String, RiskModel> riskModelEntry : riskModelMap.entrySet()) {
      final ArrayList<Individual.RiskPrediction> outputRiskPredList =
          new ArrayList<Individual.RiskPrediction>();
      /*
       * Add each individual's risk prediction from the current model to a results ArrayList.
       */
      for (final Individual.Genotypes inputGenos : inputGenoList) {
        outputRiskPredList.add(riskModelEntry.getValue().getRiskPrediction(inputGenos));
      }
      /*
       * Output the risk prediction for each individual to a single row of a tab-delimited text file
       * named "[PED/MAP root path]-[modelID].prd".
       */
      final Path modelOutFilePath =
          Paths.get(genoFilesRoot + "-" + riskModelEntry.getKey() + ".prd");
      // Log output file path for current risk model.
      pw.println(riskModelEntry.getValue().getModelName() + " --> " + modelOutFilePath.toString());
      try (PrintWriter outPw =
          new PrintWriter(Files.newBufferedWriter(modelOutFilePath, StandardCharsets.UTF_8))) {
        // Output file header row.
        outPw.print("IndivID\tModel");
        for (final String rsID : outputRiskPredList.get(0).getUsedGenotypesRsMap().keySet()) {
          outPw.print("\t" + rsID);
        }
        outPw.print("\tPI\tPIPctl");
        /*
         * Built risk models have getTimes()[i] = i, where i is age in years, so getPredCumRisk()[i]
         * is the predicted cumulative risk at age i.
         */
        for (int ageYrs = 0; ageYrs < outputRiskPredList.get(0).getPredCumRisk().length; ageYrs++) {
          outPw.format("\tPredCumRiskAge%1$d", ageYrs);
        }
        outPw.println();
        // Output individual risk predictions.
        for (final Individual.RiskPrediction riskPrediction : outputRiskPredList) {
          outPw.print(riskPrediction.getIndivID());
          outPw.print("\t" + riskPrediction.getModelName());
          for (final String geno : riskPrediction.getUsedGenotypesRsMap().values()) {
            outPw.print("\t" + geno);
          }
          outPw.print("\t" + riskPrediction.getPrognosticIndex());
          outPw.print("\t" + riskPrediction.getPrognosticIndexPctl());
          final double[] predCumRisk = riskPrediction.getPredCumRisk();
          for (int ageYrs = 0; ageYrs < predCumRisk.length; ageYrs++) {
            outPw.print("\t" + predCumRisk[ageYrs]);
          }
          outPw.println();
        }
        // Flush output file PrintWriter.
        outPw.flush();
      } // Output file PrintWriter automatically closed here.
    }
  }

  /**
   * Main method of {@code RiskPredictor} class.
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

    LOGGER.debug("Welcome to RiskPredictor.");
    LOGGER.debug("Now attempting to build options.");

    // Define command line options.
    final Option help = new Option("h", "help", false, "Print this usage information and exit.");
    final Option listModels = new Option("l", "list_models", false,
        "Print available risk model " + "information to console.");
    final Option predict = new Option("p", "predict", false,
        "Generate risk predictions for individual genotype data supplied "
            + "in the PED/MAP file pair given by <InputBase>.ped/<InputBase>.map. "
            + "Predictions will be output to <InputBase>-<modelID>.prd files, "
            + "one for each risk model specified by <modelID> located in <ModelSaveDir>"
            + "and logging information will be printed to the console.");


    final Option optModelID = new Option("m", "model_id", true,
        "Declare which model to use for generating risk predictions. Leave blank to use all "
            + "available models in <ModelSaveDir>.");
    optModelID.setArgName("modelID");

    final Option sources = new Option("s", "input_base", true,
        "Sets <InputBase> described in the --predict option. "
            + "If left unset, the builder will look in the current working directory "
            + "for a PED/MAP pair named subjects_genotypes.ped/subjects_genotypes.map.");
    sources.setArgName("InputBase");

    final Option destination = new Option("d", "model_dir", true,
        "Defines the directory in which risk model objects created by previous runs of PRISM "
            + "RiskModelBuilder are stored. If left unset, PRISM will look for risk model objects "
            + "in the current working directory.");
    destination.setArgName("ModelSaveDir");


    /*
     * Make options a mutually exclusive group so only one can be selected at a time. Make the group
     * required so that at least one of these options must be selected.
     */
    final OptionGroup modes = new OptionGroup();
    modes.addOption(help);
    modes.addOption(listModels);
    modes.addOption(predict);
    modes.setRequired(true);
    // Set command-line options.
    final Options options = new Options();
    options.addOptionGroup(modes);
    options.addOption(optModelID);
    options.addOption(sources);
    options.addOption(destination);



    LOGGER.debug("After declaring options...");
    // Try to run the application.
    try {
      // Print masthead to standard output.
      Util.printMasthead(RiskPredictor.class.getSimpleName(), stdOutPw);
      // Parse command line and echo arguments to standard output.
      final CommandLine cmd = new PosixParser().parse(options, args);
      Util.printArgs(cmd, stdOutPw);
      // Declare LinkedHashMap to load RiskModel objects.
      final LinkedHashMap<String, RiskModel> riskModelMap;

      final String modelID = cmd.getOptionValue("model_id", "");
      final String inputBase = cmd.getOptionValue("input_base", "subjects_genotypes");
      final String modelDir = cmd.getOptionValue("model_dir", System.getProperty("user.dir"));

      LOGGER.debug("Confirming, read modelID: " + modelID);
      LOGGER.debug("Confirming, read SourceDir: " + inputBase);
      LOGGER.debug("Confirming, read ModelSaveDir: " + modelDir);

      switch (modes.getSelected()) {
        case "h":
          // If in "h" mode, print usage message to standard output.
          Util.printHelp(options, RiskPredictor.class.getSimpleName(), stdOutPw);
          break;
        case "l":
          /*
           * If in "l" mode, load available risk models and print information to standard output.
           */
          riskModelMap = loadRiskModels(modelDir);
          printRiskModels(riskModelMap, stdOutPw);
          break;
        case "p":
          /*
           * If in "p" mode, load available risk models, print information to standard output, and
           * produce output files containing risk predictions for individuals in input genotype
           * files under each risk model.
           */
          riskModelMap = loadRiskModels(modelDir);

          // If a model ID has been specified, just load that model
          if (!modelID.isEmpty()) {
            LOGGER.debug("modelID specified, now loading a single risk model.");
            outputRiskPredictions(loadRiskModel(modelDir, modelID), inputBase, stdOutPw);
          } else {
            // If not, load all models
            printRiskModels(riskModelMap, stdOutPw);
            outputRiskPredictions(riskModelMap, inputBase, stdOutPw);
          }

          break;
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
      Util.printHelp(options, RiskPredictor.class.getSimpleName(), stdOutPw);
      exitCode = 1;
    } catch (final NoSuchFileException e) {
      stdOutPw.println("ERROR! File not found!");
      stdOutPw.println(
          "The following file could not be found. Please check your command and try again.");
      stdOutPw.println(WordUtils.wrap(e.getMessage(), 80));
      exitCode = 1;
    } catch (final Exception e) {
      LOGGER.error("Not a parse exception!");
      LOGGER.error("Exception thrown. Details:", e);
      /*
       * For all other exceptions, print the exception string and set the exit code to 1 (failure).
       */
      stdOutPw.println(e.toString());
      exitCode = 1;
    }
    // Flush and close the standard output PrintWriter.
    stdOutPw.flush();
    stdOutPw.close();
    // Exit with the appropriate exit code.
    System.exit(exitCode);
  }
}
