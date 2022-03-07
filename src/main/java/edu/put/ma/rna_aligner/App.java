package edu.put.ma.rna_aligner;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.lang.Thread;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang3.StringUtils;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.io.PDBFileReader;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class App {
  private static final Logger LOGGER = LoggerFactory.getLogger(App.class);

  public static void main(final String[] args) {
    CommandLine parsed_args = validateInput(args);

    AlignerConfig config = new AlignerConfig();
    try {
      config.rmsdLimit = Double.parseDouble(parsed_args.getOptionValue("rmsd", "3.5"));
      config.returnTime = Integer.parseInt(parsed_args.getOptionValue("time-limit", "300"));
      config.threads = Integer.parseInt(
          parsed_args.getOptionValue("threads", Integer.toString(Thread.activeCount())));
      // Geometric specific
      config.pairRmsdLimit = Double.parseDouble(parsed_args.getOptionValue("pair-rmsd", "0.65"));
      config.tripleRmsdLimit =
          Double.parseDouble(parsed_args.getOptionValue("triple-rmsd", "1.00"));
      // Genetic specific
      config.geometricPopulation = parsed_args.hasOption("geometric-pop");
      config.populationSize = Integer.parseInt(parsed_args.getOptionValue("pop-size", "200"));
    } catch (NumberFormatException e) {
      LOGGER.error(e.getMessage());
      System.exit(2);
    }

    final File outputDirectory = getOutputDirectory(parsed_args);

    // Parse PDB and create default coarse-grained structures.
    final PdbParser parser = new PdbParser();
    final ArrayList<Nucleotide> referenceStructure =
        parser.ParsePdbToCoarseGrained(parsed_args.getOptionValue("reference"));
    final ArrayList<Nucleotide> targetStructure =
        parser.ParsePdbToCoarseGrained(parsed_args.getOptionValue("target"));

    // Aligner function input: config, coarse-grained RNA structures.

    final boolean isSequenceDependent =
        StringUtils.equalsIgnoreCase(parsed_args.getOptionValue("mode", "seq-indep"), "seq-dep");

    AlignerOutput output = null;
    if (StringUtils.equalsIgnoreCase(
            parsed_args.getOptionValue("method", "geometric"), "geometric")) {
      GeometricAligner aligner = new GeometricAligner(
          config, referenceStructure, targetStructure, isSequenceDependent, config.rmsdLimit);
      output = aligner.calculate();
    } else {
      GeneticAligner aligner = new GeneticAligner(
          config, referenceStructure, targetStructure, isSequenceDependent, config.rmsdLimit);
      output = aligner.calculate();
    }

    final String modelNameWithoutExtension =
        FilenameUtils.removeExtension(FilenameUtils.getName(parsed_args.getOptionValue("target")));

    final StringBuilder outputStringBuilder = new StringBuilder();
    outputStringBuilder.append(String.format("Aligning mode: %s\n",
        ((isSequenceDependent) ? "sequence-dependent" : "sequence-independent")));
    outputStringBuilder.append(String.format("Maximal RMSD threshold: %.2f\n", config.rmsdLimit));
    outputStringBuilder.append(
        String.format("Residues number of reference structure: %d\n", referenceStructure.size()));
    outputStringBuilder.append(
        String.format("Residues number of model: %d\n", targetStructure.size()));

    if ((output != null) && (output.aligned > 0)) {
      outputStringBuilder.append(
          String.format("Number of aligned nucleotides: %d\n", output.aligned));
      outputStringBuilder.append(String.format("RMSD score: %.3f\n", output.rmsd));
      outputStringBuilder.append(
          String.format("Processing time [ms]: %d\n", output.processingTime));
    } else {
      outputStringBuilder.append(String.format("Alignment is not found.\n"));
    }

    saveDataToFile(new StringBuilder(outputDirectory.getAbsolutePath())
                       .append(File.separator)
                       .append(modelNameWithoutExtension)
                       .append("-output.txt")
                       .toString(),
        outputStringBuilder.toString());

    if ((output != null) && (output.aligned > 0)) {
      final StringBuilder alignedResiduesMappingStringBuilder =
          new StringBuilder("REF  \t<->\tMODEL\n");
      for (int i = 0; i < output.referenceIndexes.size(); ++i) {
        final int targetIndex = output.targetMapping.get(i).intValue();
        alignedResiduesMappingStringBuilder.append(String.format("%s\t<->\t%s\n",
            referenceStructure.get(output.referenceIndexes.get(i)).getKey(),
            (targetIndex > -1) ? targetStructure.get(targetIndex).getKey() : "-"));
      }
      saveDataToFile(new StringBuilder(outputDirectory.getAbsolutePath())
                         .append(File.separator)
                         .append(modelNameWithoutExtension)
                         .append("-residue-mapping.txt")
                         .toString(),
          alignedResiduesMappingStringBuilder.toString());

      saveDataToFile(new StringBuilder(outputDirectory.getAbsolutePath())
                         .append(File.separator)
                         .append(modelNameWithoutExtension)
                         .append("-sequence-alignment.txt")
                         .toString(),
          getAlignment(referenceStructure, targetStructure, output));

      final Structure model = readStructure(parsed_args.getOptionValue("target"));
      rotateAndShift(model, output.superimposer);
      saveDataToFile(new StringBuilder(outputDirectory.getAbsolutePath())
                         .append(File.separator)
                         .append(modelNameWithoutExtension)
                         .append("-superimposed.pdb")
                         .toString(),
          model.toPDB());

      final StringBuilder refModelSuperimpositionStringBuilder = new StringBuilder();
      refModelSuperimpositionStringBuilder.append(
          "MODEL        1                                              \n");
      final Structure reference = readStructure(parsed_args.getOptionValue("reference"));
      refModelSuperimpositionStringBuilder.append(reference.toPDB());
      refModelSuperimpositionStringBuilder.append(
          "ENDMDL                                                      \n");
      refModelSuperimpositionStringBuilder.append(
          "MODEL        2                                              \n");
      refModelSuperimpositionStringBuilder.append(model.toPDB());
      refModelSuperimpositionStringBuilder.append(
          "ENDMDL                                                      \n");
      refModelSuperimpositionStringBuilder.append(
          "END                                                         ");

      final String referenceNameWithoutExtension = FilenameUtils.removeExtension(
          FilenameUtils.getName(parsed_args.getOptionValue("reference")));
      saveDataToFile(new StringBuilder(outputDirectory.getAbsolutePath())
                         .append(File.separator)
                         .append(referenceNameWithoutExtension)
                         .append("-")
                         .append(modelNameWithoutExtension)
                         .append("-superimposed.pdb")
                         .toString(),
          refModelSuperimpositionStringBuilder.toString());
    }
  }

  private static final File getOutputDirectory(final CommandLine args) {
    final File path =
        Paths.get(args.getOptionValue("output", args.getOptionValue("target"))).toFile();

    if (path.exists()) {
      if (path.isDirectory()) {
        // Should be 'output' flag.
        return path;
      } else {
        final File targetDir = path.getAbsoluteFile().getParentFile();
        if (targetDir.exists() && targetDir.isDirectory()) {
          // Should be 'target' flag. Extract directory from it.
          return targetDir;
        }
      }
    }

    LOGGER.error("Output path directory is not a proper folder path");
    System.exit(3);
    return null;
  }

  private static final CommandLine validateInput(final String[] args) {
    Options options = new Options();

    options.addOption(OptionBuilder.withLongOpt("reference")
                          .withDescription("Reference structure in .pdb format.")
                          .hasArg()
                          .withArgName("reference.pdb")
                          .isRequired()
                          .create('r'));

    options.addOption(OptionBuilder.withLongOpt("target")
                          .withDescription("Target structure in .pdb format.")
                          .hasArg()
                          .withArgName("target.pdb")
                          .isRequired()
                          .create('t'));

    options.addOption(
        OptionBuilder.withLongOpt("output")
            .withDescription("(optional) Output directory for all results and alignements.\n"
                + "Default: use directory of target structure")
            .hasArg()
            .withArgName("path")
            .create('o'));

    options.addOption(OptionBuilder.withLongOpt("method")
                          .withDescription("(optional) Method used to align input structures.\n"
                              + "Available: geometric, genetic\n"
                              + "Default: geometric")
                          .hasArg()
                          .withArgName("method")
                          .create('m'));

    options.addOption(
        OptionBuilder.withLongOpt("mode")
            .withDescription("(optional) Aligning mode used with each method. Can be either "
                + "sequence independent or sequence dependent\n"
                + "Available: seq-indep, seq-dep\n"
                + "Default: seq-indep")
            .hasArg()
            .withArgName("aligning-mode")
            .create());

    options.addOption(
        OptionBuilder.withLongOpt("rmsd")
            .withDescription("(optional) Maximum RMSD (in Ångström) that the aligned fragment will "
                + "not exceed.\n"
                + "Default: 3.5")
            .hasArg()
            .withType(Float.class)
            .withArgName("rmsd")
            .create());

    options.addOption(
        OptionBuilder.withLongOpt("time-limit")
            .withDescription("(optional) Maximum execution time in seconds before program returns. "
                + "Execution time can be well below set value if no improvements are found for long "
                + "time.\n"
                + "Default: 300")
            .hasArg()
            .withType(Integer.class)
            .withArgName("seconds")
            .create());

    options.addOption(
        OptionBuilder.withLongOpt("threads")
            .withDescription("(optional) Number of threads used by algoritm. Easy way to speedup "
                + "the processing.\n"
                + "Default: all system threads")
            .hasArg()
            .withType(Integer.class)
            .withArgName("threads")
            .create());

    options.addOption(
        OptionBuilder.withLongOpt("pair-rmsd")
            .withDescription("(optional) Maximum RMSD (in Ångström) that cores with 2 nucleotides "
                + "will not exceed. Increase leads to wider and longer search. Must be lower or "
                + "equal to triple-rmsd.\n"
                + "Default: 0.65")
            .hasArg()
            .withType(Float.class)
            .withArgName("rmsd")
            .create());

    options.addOption(
        OptionBuilder.withLongOpt("triple-rmsd")
            .withDescription("(optional) Maximum RMSD (in Ångström) that cores with 3 nucleotides "
                + "will not exceed. Increase leads to wider and longer search. Must be higher or "
                + "equal to pair-rmsd.\n"
                + "Default: 1.0")
            .hasArg()
            .withType(Float.class)
            .withArgName("rmsd")
            .create());

    options.addOption(
        OptionBuilder.withLongOpt("pop-size")
            .withDescription("(optional) Population size for each generation and thread\n"
                + "Default: 200")
            .hasArg()
            .withType(Integer.class)
            .withArgName("size")
            .create());

    options.addOption(
        OptionBuilder.withLongOpt("geometric-pop")
            .withDescription("(optional) Generate initial population using first results obrained "
                + "from the geometric algorithm.\n")
            .create());

    CommandLineParser parser = new DefaultParser();
    HelpFormatter formatter = new HelpFormatter();
    CommandLine cmd = null;

    try {
      cmd = parser.parse(options, args);
    } catch (ParseException e) {
      LOGGER.error(e.getMessage());
      formatter.printHelp(100,
          "java -jar rna-hugs.jar -r <reference.pdb> -t <target.pdb> [OPTIONS]", "", options, "");
      System.exit(1);
    }

    return cmd;
  }

  private static final Structure readStructure(final String path) {
    final PDBFileReader reader = new PDBFileReader();
    Structure structure = null;
    try (FileInputStream fis = new FileInputStream(path)) {
      structure = reader.getStructure(fis);
    } catch (IOException e) {
      LOGGER.error(e.getMessage(), e);
    }
    return structure;
  }

  private static final void saveDataToFile(final String outputPath, final String data) {
    try {
      Files.write(Paths.get(outputPath), data.getBytes());
    } catch (IOException e) {
      LOGGER.error(e.getMessage(), e);
    }
  }

  private static final void rotateAndShift(
      final Structure structure, final Superimposer superimposer) {
    final List<Chain> model = structure.getModel(0);
    for (Chain chain : model) {
      for (Group residue : chain.getAtomGroups()) {
        for (org.biojava.nbio.structure.Atom atom : residue.getAtoms()) {
          rotateAndShiftAtom(atom, superimposer);
        }
      }
    }
  }

  private static final void rotateAndShiftAtom(final Atom atom, final Superimposer superimposer) {
    Calc.rotate(atom, superimposer.getRotation());
    shift(atom, superimposer.translation);
  }

  public static final void shift(final Atom a, final Coordinates shift) {
    a.setX(a.getX() + shift.x);
    a.setY(a.getY() + shift.y);
    a.setZ(a.getZ() + shift.z);
  }

  private static final String getAlignment(final ArrayList<Nucleotide> referenceStructure,
      final ArrayList<Nucleotide> targetStructure, final AlignerOutput output) {
    final StringBuilder alignment = new StringBuilder();
    StringBuilder referenceStringBuilder = new StringBuilder();
    StringBuilder targetStringBuilder = new StringBuilder();
    StringBuilder linkerStringBuilder = new StringBuilder();
    int count = 0;
    for (int i = 0; i < referenceStructure.size(); i++) {
      referenceStringBuilder.append(referenceStructure.get(i).getCode());
      final int referenceIndex = output.referenceIndexes.indexOf(Integer.valueOf(i));
      if (referenceIndex > -1) {
        final int targetIndex = output.targetMapping.get(referenceIndex);
        if (targetIndex > -1) {
          targetStringBuilder.append(targetStructure.get(targetIndex).getCode());
          linkerStringBuilder.append('|');
        } else {
          targetStringBuilder.append('-');
          linkerStringBuilder.append(' ');
        }
      } else {
        targetStringBuilder.append('-');
        linkerStringBuilder.append(' ');
      }
      count++;
      if (count % 80 == 0) {
        alignment.append("REF:   ").append(referenceStringBuilder.toString()).append('\n');
        alignment.append("       ").append(linkerStringBuilder.toString()).append('\n');
        alignment.append("MODEL: ").append(targetStringBuilder.toString()).append('\n');
        alignment.append('\n');
        referenceStringBuilder = new StringBuilder();
        targetStringBuilder = new StringBuilder();
        linkerStringBuilder = new StringBuilder();
        count = 0;
      }
    }
    alignment.append("REF:   ").append(referenceStringBuilder.toString()).append('\n');
    alignment.append("       ").append(linkerStringBuilder.toString()).append('\n');
    alignment.append("MODEL: ").append(targetStringBuilder.toString()).append('\n');
    return alignment.toString();
  }
}
