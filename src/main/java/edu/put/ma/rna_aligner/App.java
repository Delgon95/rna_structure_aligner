package edu.put.ma.rna_aligner;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.lang.Thread;
import java.nio.charset.StandardCharsets;
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
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang3.StringUtils;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.io.PDBFileReader;
import org.biojava.nbio.structure.jama.Matrix;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import pl.poznan.put.pdb.PdbAtomLine;
import pl.poznan.put.pdb.analysis.CifModel;
import pl.poznan.put.pdb.analysis.CifParser;
import pl.poznan.put.pdb.analysis.DefaultPdbModel;
import pl.poznan.put.pdb.analysis.ImmutableDefaultCifModel;
import pl.poznan.put.pdb.analysis.ImmutableDefaultPdbModel;
import pl.poznan.put.pdb.analysis.PdbModel;
import pl.poznan.put.pdb.analysis.PdbParser;
import pl.poznan.put.pdb.analysis.PdbResidue;
import pl.poznan.put.pdb.analysis.ResidueCollection;

public class App {
  private static final Logger LOGGER = LoggerFactory.getLogger(App.class);

  public static void main(final String[] args) {
    CommandLine parsed_args = validateInput(args);

    AlignerConfig config = new AlignerConfig();
    try {
      config.rmsdLimit = Double.parseDouble(parsed_args.getOptionValue("rmsd", "3.5"));
      config.returnTime = Integer.parseInt(parsed_args.getOptionValue("time-limit", "300"));
      config.threads = Integer.parseInt(parsed_args.getOptionValue(
          "threads", Integer.toString(Runtime.getRuntime().availableProcessors())));
      // Geometric specific
      config.pairRmsdLimit = Double.parseDouble(parsed_args.getOptionValue("pair-rmsd", "0.65"));
      config.tripleRmsdLimit = Double.parseDouble(parsed_args.getOptionValue("triple-rmsd", "1.0"));
      // Genetic specific
      config.geometricPopulation = parsed_args.hasOption("geometric-pop");
      config.populationSize = Integer.parseInt(parsed_args.getOptionValue("pop-size", "200"));
    } catch (NumberFormatException e) {
      LOGGER.error(e.getMessage());
      System.exit(2);
    }


    final File outputDirectory = getOutputDirectory(parsed_args);

    // Parse PDB and create default coarse-grained structures.
    final StructureParser parser = new StructureParser();
    final ArrayList<Nucleotide> referenceStructure = parser.StructureToCoarseGrained(
        parsed_args.getOptionValue("reference"), parsed_args.getOptionValue("pop-size", "auto"),
        parsed_args.hasOption("allow-incomplete"));
    final ArrayList<Nucleotide> targetStructure = parser.StructureToCoarseGrained(
        getModelArg(parsed_args) , parsed_args.getOptionValue("pop-size", "auto"),
        parsed_args.hasOption("allow-incomplete"));

    // Aligner function input: config, coarse-grained RNA structures.

    final boolean isSequenceDependent =
        StringUtils.equalsIgnoreCase(parsed_args.getOptionValue("mode", "seq-indep"), "seq-dep");
    final boolean respectOrder = parsed_args.hasOption("respect-order");

    AlignerOutput output = null;
    if (StringUtils.equalsIgnoreCase(
            parsed_args.getOptionValue("method", "geometric"), "geometric")) {
      GeometricAligner aligner = new GeometricAligner(
          config, referenceStructure, targetStructure, isSequenceDependent, config.rmsdLimit, respectOrder);
      output = aligner.calculate();
    } else {
      GeneticAligner aligner = new GeneticAligner(
          config, referenceStructure, targetStructure, isSequenceDependent, config.rmsdLimit, respectOrder);
      output = aligner.calculate();
    }

    final String modelNameWithoutExtension =
        FilenameUtils.removeExtension(FilenameUtils.getName(getModelArg(parsed_args)));

    final StringBuilder outputStringBuilder = new StringBuilder();
    outputStringBuilder.append(String.format("Alignment mode: %s\n",
        ((isSequenceDependent) ? "sequence-dependent" : "sequence-independent")));
    outputStringBuilder.append(String.format("RMSD threshold [Å]: %.2f\n", config.rmsdLimit));
    outputStringBuilder.append(
        String.format("Reference structure size [nts]: %d\n", referenceStructure.size()));
    outputStringBuilder.append(
        String.format("Aligned structure size [nts]: %d\n", targetStructure.size()));


    if ((output != null) && (output.aligned > 0)) {
      outputStringBuilder.append(
          String.format("Number of aligned residues: %d\n", output.aligned));
      float min_size = Math.min(referenceStructure.size(), targetStructure.size());
      int percentage = 0;
      // Should always happen.
      if (min_size > 0) {
        percentage = new Double(Math.floor(Float.valueOf(output.aligned) / min_size * 100.0f)).intValue();
      }
      outputStringBuilder.append(
          String.format("Percentage of aligned residues [%%]: %d\n", percentage));
      outputStringBuilder.append(String.format("RMSD of aligned fragments [Å]: %.3f\n", output.rmsd));
      outputStringBuilder.append(
          String.format("Processing time [sec]: %d.%03d\n", output.processingTime / 1000, output.processingTime % 1000));
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
            referenceStructure.get(output.referenceIndexes.get(i)).getDelimitedChainNumber(":"),
            (targetIndex > -1) ? targetStructure.get(targetIndex).getDelimitedChainNumber(":") : "-"));
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

      ImmutableDefaultPdbModel model = getRotatedStructure(getModelArg(parsed_args),
          parsed_args.getOptionValue("pop-size", "auto"), output.superimposer);
      saveDataToFile(new StringBuilder(outputDirectory.getAbsolutePath())
                         .append(File.separator)
                         .append(modelNameWithoutExtension)
                         .append("-superimposed.pdb")
                         .toString(),
          model.toPdb());

      saveDataToFile(new StringBuilder(outputDirectory.getAbsolutePath())
                         .append(File.separator)
                         .append(modelNameWithoutExtension)
                         .append("-superimposed.cif")
                         .toString(),
          model.toCif());
    }
  }

  private static final String getModelArg(final CommandLine args) {
    return args.getOptionValue("model", args.getOptionValue("target"));
  }

  private static final File getOutputDirectory(final CommandLine args) {
    final File path =
        Paths.get(args.getOptionValue("output", getModelArg(args))).toFile();

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

    options.addOption(
        OptionBuilder.withLongOpt("reference")
            .withDescription(
                "Reference structure in .pdb/.cif format. Can force format with --input-format")
            .hasArg()
            .withArgName("reference.pdb")
            .isRequired()
            .create('r'));

    options.addOption(
        OptionBuilder.withLongOpt("target")
            .withDescription(
                "Same as --model (Deprecated)")
            .hasArg()
            .withArgName("model.pdb")
            .create('t'));
    options.addOption(
        OptionBuilder.withLongOpt("model")
            .withDescription(
                "Model structure in .pdb/.cif format. Can force format with --input-format")
            .hasArg()
            .withArgName("model.pdb")
            .create('m'));
    options.addOption(
        OptionBuilder.withLongOpt("input-format")
            .withDescription(
                "(optional) Format type of both input structures. Auto allows for different "
                + "formats between model and reference\n"
                + "Available: auto, pdb, cif\n"
                + "Default: auto")
            .hasArg()
            .withArgName("format")
            .create());

    options.addOption(
        OptionBuilder.withLongOpt("output")
            .withDescription("(optional) Output directory for all results and alignements.\n"
                + "Default: use directory of model structure")
            .hasArg()
            .withArgName("path")
            .create('o'));

    options.addOption(OptionBuilder.withLongOpt("method")
                          .withDescription("(optional) Method used to align input structures.\n"
                              + "Available: geometric, genetic\n"
                              + "Default: geometric")
                          .hasArg()
                          .withArgName("method")
                          .create());

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

    options.addOption(
        OptionBuilder.withLongOpt("respect-order")
            .withDescription("(optional) The resulting alignemnt will respect nucleotide order "
                + "in the structure. This might make shorter but more biologicaly accurate "
                + "alignments.\n")
            .create());


    options.addOption(
        OptionBuilder.withLongOpt("allow-incomplete")
            .withDescription("(optional) Allow usage of incomplete atoms in coarse-grained "
                + "structure creation. By default, all of the atoms specified in the code "
                + "are required to include molecule in calculations.\n")
            .create());

    CommandLineParser parser = new DefaultParser();
    HelpFormatter formatter = new HelpFormatter();
    CommandLine cmd = null;

    try {
      cmd = parser.parse(options, args);
    } catch (ParseException e) {
      LOGGER.error(e.getMessage());
      formatter.printHelp(100,
          "java -jar rna-hugs.jar -r <reference.pdb> -m <model.pdb> [OPTIONS]", "", options, "");
      System.exit(1);
    }

    if (! cmd.hasOption("model") && ! cmd.hasOption("target")) {
      LOGGER.error("Model structure is required");
      formatter.printHelp(100,
          "java -jar rna-hugs.jar -r <reference.pdb> -m <model.pdb> [OPTIONS]", "", options, "");
      System.exit(1);
    }

    return cmd;
  }

  private static final ImmutableDefaultPdbModel getRotatedStructure(
      final String path, final String inputType, final Superimposer superimposer) {
    final ResidueCollection originalModel = readStructure(path, inputType);
    final List<PdbAtomLine> shiftedAtoms = getRotatedAndShifted(originalModel, superimposer);
    return ImmutableDefaultPdbModel.copyOf((DefaultPdbModel) originalModel).withAtoms(shiftedAtoms);
  }

  private static final ResidueCollection readStructure(final String path, final String inputType) {
    boolean isPdb = true;

    final File file = new File(path);
    try {
      final String structureContent = FileUtils.readFileToString(file, StandardCharsets.UTF_8);

      if (inputType == "auto") {
        // Try Detect if the file is mmCIF or PDB.
        // Search for _loop in the file. If it is present we can assume it is CIF format.
        if (structureContent.indexOf("_loop") != -1 || structureContent.indexOf("loop_") != -1) {
          isPdb = false;
        } // else remains isPdb = true;
      } else if (inputType == "cif") {
        isPdb = false;
      } // else remains isPdb = true;

      // parse the data
      final PdbParser pdbParser = new PdbParser(false);
      final List<PdbModel> pdbModels = (isPdb) ? pdbParser.parse(structureContent) : null;
      final List<CifModel> cifModels = (!isPdb) ? CifParser.parse(structureContent) : null;

      return (isPdb) ? pdbModels.get(0) : pdbParser.parse(cifModels.get(0).toPdb()).get(0);
    } catch (IOException e) {
      e.printStackTrace();
    }
    return null;
  }

  private static final void saveDataToFile(final String outputPath, final String data) {
    try {
      Files.write(Paths.get(outputPath), data.getBytes());
    } catch (IOException e) {
      LOGGER.error(e.getMessage(), e);
    }
  }

  private static List<PdbAtomLine> getRotatedAndShifted(
      ResidueCollection structure, final Superimposer superimposer) {
    List<PdbAtomLine> shiftedAtoms = new ArrayList<PdbAtomLine>();

    for (PdbResidue residue : structure.residues()) {
      for (PdbAtomLine atom : residue.atoms()) {
        shiftedAtoms.add(rotateAndShiftAtom(atom, superimposer));
      }
    }

    return shiftedAtoms;
  }

  private static final PdbAtomLine rotateAndShiftAtom(
      PdbAtomLine atom, final Superimposer superimposer) {
    org.biojava.nbio.structure.Atom bioJavaAtom = atom.toBioJavaAtom();
    Calc.rotate(bioJavaAtom, superimposer.getRotation());
    shift(bioJavaAtom, superimposer.translation);
    return PdbAtomLine.fromBioJavaAtom(bioJavaAtom);
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
