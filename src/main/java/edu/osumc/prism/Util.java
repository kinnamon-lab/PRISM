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

import java.io.PrintWriter;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.text.WordUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * PRISM utility functions.
 * <p>
 * This package-private class contains several static methods used in multiple other classes.
 */
final class Util {

  // Initialize logger
  private static final Logger LOGGER = LogManager.getLogger();

  /** Prevents construction of {@code Util} class instances. */
  private Util() {};
  
  /**
   * Prints program masthead for a specific tool to a {@code PrintWriter}.
   *
   * @param tool {@code String} identifying the tool being used (i.e., a class name)
   * @param pw {@code PrintWriter} for output
   */
  protected static final void printMasthead(final String tool, final PrintWriter pw) {
    final String version = Package.getPackage("edu.osumc.prism").getImplementationVersion();
    pw.println(StringUtils.repeat("=", 79));
    pw.println(StringUtils.center("PRISM - Polygenic RIsk Score Modeler", 79));
    pw.println(StringUtils.center("* Version " + version + " *", 79));
    pw.println(StringUtils.center(tool + " Tool", 79));
    pw.println(StringUtils.center(
        "-- Copyright 2014-2017 The Ohio State " + "University Wexner Medical Center --", 79));
    pw.println(StringUtils.center("Authors: Daniel D. Kinnamon, PhD; Carl A. Starkey, PhD", 79));
    pw.println(StringUtils.center("Licensed under the Apache License, Version 2.0", 79));
    pw.println(StringUtils.repeat("=", 79));
    pw.println(WordUtils.wrap("NOTE: This software is distributed on "
        + "an \"AS IS\" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, "
        + "either express or implied.", 79));
    pw.println(StringUtils.repeat("-", 79));
  }

  /**
   * Prints a command-line help message to a {@code PrintWriter}.
   *
   * @param options {@link Options} object containing possible options for the command line
   * @param tool {@code String} identifying the tool being used (i.e., a class name)
   * @param pw {@code PrintWriter} for output
   */
  protected static final void printHelp(final Options options, final String tool,
      final PrintWriter pw) {
    final String version = Package.getPackage("edu.osumc.prism").getImplementationVersion();
    new HelpFormatter().printHelp(pw, HelpFormatter.DEFAULT_WIDTH,
        "java -cp prism-" + version + ".jar edu.osumc.prism." + tool, "", options,
        HelpFormatter.DEFAULT_LEFT_PAD, HelpFormatter.DEFAULT_DESC_PAD, "", true);
  }

  /**
   * Prints the arguments with which the program was invoked to a {@code PrintWriter}.
   *
   * @param cmd {@link CommandLine} object containing parsed command line
   * @param pw {@code PrintWriter} for output
   */
  protected static final void printArgs(final CommandLine cmd, final PrintWriter pw) {
    pw.println();
    pw.println("Invoked with arguments:");
    for (final Option option : cmd.getOptions()) {
      pw.println(HelpFormatter.DEFAULT_LONG_OPT_PREFIX + option.getLongOpt()
          + (option.hasArg() ? " " + option.getValue() : ""));
    }
    pw.println();
  }

}
