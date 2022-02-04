package edu.put.ma.rna_aligner;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
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
    final AlignerConfig config = new AlignerConfig();

    validateInput(args, config);

    final File outputDirectory = getOutputDirectory(args);

    // Parse PDB and create default coarse-grained structures.
    final PdbParser parser = new PdbParser();
    final ArrayList<Nucleotide> referenceStructure = parser.ParsePdbToCoarseGrained(args[1]);
    final ArrayList<Nucleotide> targetStructure = parser.ParsePdbToCoarseGrained(args[2]);

    // Aligner function input: config, coarse-grained RNA structures.

    boolean isSequenceDependent = true;
    if (StringUtils.equalsIgnoreCase("seq-indep", args[3]))
      isSequenceDependent = false;

    final double rmsdLimit = Double.parseDouble(args[4]);

    AlignerOutput output = null;
    if (args[0].equals("genetic")) {
      GeneticAligner aligner = new GeneticAligner(
          config, referenceStructure, targetStructure, isSequenceDependent, rmsdLimit);
      output = aligner.calculate();
    } else {
      GeometricAligner aligner = new GeometricAligner(
          config, referenceStructure, targetStructure, isSequenceDependent, rmsdLimit);
      output = aligner.calculate();
    }

    final String modelNameWithoutExtension =
        FilenameUtils.removeExtension(FilenameUtils.getName(args[2]));

    final StringBuilder outputStringBuilder = new StringBuilder();
    outputStringBuilder.append(String.format("Aligning mode: %s\n",
        ((isSequenceDependent) ? "sequence-dependent" : "sequence-independent")));
    outputStringBuilder.append(String.format("Maximal RMSD threshold: %.2f\n", rmsdLimit));
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
    } else
      outputStringBuilder.append(String.format("Alignment is not found.\n"));

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

      final Structure model = readStructure(args[2]);
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
      final Structure reference = readStructure(args[1]);
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

      final String referenceNameWithoutExtension =
          FilenameUtils.removeExtension(FilenameUtils.getName(args[1]));
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

  private static final File getOutputDirectory(final String[] args) {
    File outputDirectory = null;
    if (args.length == 6) {
      final File userOutputFile = Paths.get(args[5]).toFile();
      if ((userOutputFile.exists()) && (userOutputFile.isDirectory())) {
        outputDirectory = userOutputFile;
      }
    }
    if (outputDirectory == null) {
      final File modelFile = Paths.get(args[2]).toFile();
      final File modelDirectory = modelFile.getParentFile();
      if ((modelDirectory.exists()) && (modelDirectory.isDirectory())) {
        outputDirectory = modelDirectory;
      }
    }

    if (outputDirectory == null) {
      LOGGER.error("Output directory must be defined!");
    }
    return outputDirectory;
  }

  private static final void validateInput(final String[] args, final AlignerConfig _config) {
    if (args.length < 5 || args.length > 6) {
      LOGGER.info(
          "Usage: <method> <reference.pdb> <target.pdb> <aligning-mode> <rmsd-threshold> <output-directory> (optional)");
      LOGGER.info("Provided methods: geometric, genetic");
      LOGGER.info(
          "Provided aligning modes: sequence-dependent (seq-dep), sequence-independent (seq-indep). By default sequence-dependent is used.");
      LOGGER.info(
          "RMSD threshold: maximal RMSD score of resultant alignment. By default 3.5A is used.");
      System.exit(1);
    }

    if ((StringUtils.equalsIgnoreCase(args[0], "geometric"))
        || (StringUtils.equalsIgnoreCase(args[0], "genetic"))) {
      if (!StringUtils.isAllLowerCase(args[0])) {
        args[0] = StringUtils.lowerCase(args[0]);
      }
    } else {
      LOGGER.warn(String.format("Invalid method: %s", args[0]));
      LOGGER.warn("By default geometric is used.");
      args[0] = "geometric";
    }

    final File referenceFile = Paths.get(args[1]).toFile();
    if (!((referenceFile.exists()) && (referenceFile.isFile()))) {
      LOGGER.error("Reference structure file must be defined!");
      System.exit(1);
    }

    final File modelFile = Paths.get(args[2]).toFile();
    if (!((modelFile.exists()) && (modelFile.isFile()))) {
      LOGGER.error("Model file must be defined!");
      System.exit(1);
    }

    if ((StringUtils.equalsIgnoreCase(args[3], "seq-dep"))
        || (StringUtils.equalsIgnoreCase(args[3], "seq-indep"))) {
      if (!StringUtils.isAllLowerCase(args[3])) {
        args[3] = StringUtils.lowerCase(args[3]);
      }
    } else {
      LOGGER.warn(String.format("Invalid aligning mode: %s", args[3]));
      LOGGER.warn("By default seq-dep is used.");
      args[3] = "seq-dep";
    }

    try {
      Double.parseDouble(StringUtils.trim(args[4]));
    } catch (Exception e) {
      LOGGER.warn(String.format("Invalid RMSD threshold: %s", args[4]));
      LOGGER.warn("By default 3.5A is used.");
      args[4] = String.valueOf(_config.rmsdLimit);
    }
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
