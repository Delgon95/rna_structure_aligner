package edu.put.ma.rna_aligner;

import static java.util.Map.entry;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.StringUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import pl.poznan.put.pdb.PdbAtomLine;
import pl.poznan.put.pdb.analysis.MoleculeType;
import pl.poznan.put.pdb.analysis.CifModel;
import pl.poznan.put.pdb.analysis.CifParser;
import pl.poznan.put.pdb.analysis.PdbModel;
import pl.poznan.put.pdb.analysis.PdbParser;

public class StructureParser {
  private static final Logger LOGGER = LoggerFactory.getLogger(PdbParser.class);

  private static final ArrayList<String> rest_atoms =
      new ArrayList<String>(Arrays.asList("P", "O5'", "C5'"));
  private static final ArrayList<String> ryboze_atoms =
      new ArrayList<String>(Arrays.asList("C4'", "C3'", "C2'", "C1'", "O4'"));
  private static final ArrayList<String> U_atoms =
      new ArrayList<String>(Arrays.asList("N1", "C6", "C5", "C4", "O4", "N3", "C2"));
  private static final ArrayList<String> C_atoms =
      new ArrayList<String>(Arrays.asList("N1", "C6", "C5", "C4", "N3", "C2", "O2"));
  private static final ArrayList<String> A_atoms =
      new ArrayList<String>(Arrays.asList("N9", "C8", "N7", "C5", "C6", "N1", "C2", "N3", "C4"));
  private static final ArrayList<String> G_atoms =
      new ArrayList<String>(Arrays.asList("N9", "C8", "N7", "C5", "C6", "N1", "C2", "N3", "C4"));

  private static final Map<String, String> modified_3to1 = Map.ofEntries(entry("A23", "A"),
      entry("A2L", "A"), entry("A2M", "A"), entry("A39", "A"), entry("A3P", "A"), entry("A44", "A"),
      entry("A5O", "A"), entry("A6A", "A"), entry("A7E", "A"), entry("A9Z", "A"), entry("ADI", "A"),
      entry("ADP", "A"), entry("AET", "A"), entry("AMD", "A"), entry("AMO", "A"), entry("AP7", "A"),
      entry("AVC", "A"), entry("MA6", "A"), entry("MAD", "A"), entry("MGQ", "A"), entry("MIA", "A"),
      entry("MTU", "A"), entry("M7A", "A"), entry("26A", "A"), entry("2MA", "A"), entry("6IA", "A"),
      entry("6MA", "A"), entry("6MC", "A"), entry("6MP", "A"), entry("6MT", "A"), entry("6MZ", "A"),
      entry("6NW", "A"), entry("F3N", "A"), entry("N79", "A"), entry("RIA", "A"), entry("V3L", "A"),
      entry("ZAD", "A"), entry("31H", "A"), entry("31M", "A"), entry("7AT", "A"), entry("O2Z", "A"),
      entry("SRA", "A"), entry("00A", "A"), entry("45A", "A"), entry("8AN", "A"), entry("LCA", "A"),
      entry("P5P", "A"), entry("PPU", "A"), entry("PR5", "A"), entry("PU",  "A"), entry("T6A", "A"),
      entry("TBN", "A"), entry("TXD", "A"), entry("TXP", "A"), entry("12A", "A"), entry("1MA", "A"),
      entry("5FA", "A"), entry("A6G", "G"), entry("E6G", "G"), entry("E7G", "G"), entry("EQ4", "G"),
      entry("IG",  "G"), entry("IMP", "G"), entry("M2G", "G"), entry("MGT", "G"), entry("MGV", "G"),
      entry("MHG", "G"), entry("QUO", "G"), entry("YG",  "G"), entry("YYG", "G"), entry("23G", "G"),
      entry("2EG", "G"), entry("2MG", "G"), entry("2SG", "G"), entry("B8K", "G"), entry("B8W", "G"),
      entry("B9B", "G"), entry("BGH", "G"), entry("N6G", "G"), entry("RFJ", "G"), entry("ZGU", "G"),
      entry("7MG", "G"), entry("CG1", "G"), entry("G1G", "G"), entry("G25", "G"), entry("G2L", "G"),
      entry("G46", "G"), entry("G48", "G"), entry("G7M", "G"), entry("GAO", "G"), entry("GDO", "G"),
      entry("GDP", "G"), entry("GH3", "G"), entry("GNG", "G"), entry("GOM", "G"), entry("GRB", "G"),
      entry("GTP", "G"), entry("KAG", "G"), entry("KAK", "G"), entry("O2G", "G"), entry("OMG", "G"),
      entry("8AA", "G"), entry("8OS", "G"), entry("LG",  "G"), entry("PGP", "G"), entry("P7G", "G"),
      entry("TPG", "G"), entry("TG",  "G"), entry("XTS", "G"), entry("102", "G"), entry("18M", "G"),
      entry("1MG", "G"), entry("A5M", "C"), entry("A6C", "C"), entry("E3C", "C"), entry("IC",  "C"),
      entry("M4C", "C"), entry("M5M", "C"), entry("6OO", "C"), entry("B8Q", "C"), entry("B8T", "C"),
      entry("B9H", "C"), entry("JMH", "C"), entry("N5M", "C"), entry("RPC", "C"), entry("RSP", "C"),
      entry("RSQ", "C"), entry("ZBC", "C"), entry("ZCY", "C"), entry("73W", "C"), entry("C25", "C"),
      entry("C2L", "C"), entry("C31", "C"), entry("C43", "C"), entry("C5L", "C"), entry("CBV", "C"),
      entry("CCC", "C"), entry("CH",  "C"), entry("CSF", "C"), entry("OMC", "C"), entry("S4C", "C"),
      entry("4OC", "C"), entry("LC",  "C"), entry("LHH", "C"), entry("LV2", "C"), entry("PMT", "C"),
      entry("TC",  "C"), entry("10C", "C"), entry("1SC", "C"), entry("5HM", "C"), entry("5IC", "C"),
      entry("5MC", "C"), entry("A6U", "U"), entry("IU",  "U"), entry("I4U", "U"), entry("MEP", "U"),
      entry("MNU", "U"), entry("U25", "U"), entry("U2L", "U"), entry("U2P", "U"), entry("U31", "U"),
      entry("U34", "U"), entry("U36", "U"), entry("U37", "U"), entry("U8U", "U"), entry("UAR", "U"),
      entry("UBB", "U"), entry("UBD", "U"), entry("UD5", "U"), entry("UPV", "U"), entry("UR3", "U"),
      entry("URD", "U"), entry("US5", "U"), entry("UZR", "U"), entry("UMO", "U"), entry("U23", "U"),
      entry("2AU", "U"), entry("2MU", "U"), entry("2OM", "U"), entry("B8H", "U"), entry("FHU", "U"),
      entry("FNU", "U"), entry("F2T", "U"), entry("RUS", "U"), entry("ZBU", "U"), entry("3AU", "U"),
      entry("3ME", "U"), entry("3MU", "U"), entry("3TD", "U"), entry("70U", "U"), entry("75B", "U"),
      entry("CNU", "U"), entry("OMU", "U"), entry("ONE", "U"), entry("S4U", "U"), entry("SSU", "U"),
      entry("SUR", "U"), entry("4SU", "U"), entry("85Y", "U"), entry("DHU", "U"), entry("H2U", "U"),
      entry("LHU", "U"), entry("PSU", "U"), entry("PYO", "U"), entry("P4U", "U"), entry("T31", "U"),
      entry("125", "U"), entry("126", "U"), entry("127", "U"), entry("1RN", "U"), entry("5BU", "U"),
      entry("5FU", "U"), entry("5MU", "U"), entry("9QV", "U"), entry("5GP", "G"));

  public ArrayList<Nucleotide> StructureToCoarseGrained(
      final String filename, final String inputType, final boolean allowIncomplete) {
    ArrayList<Nucleotide> result = new ArrayList<Nucleotide>();
    boolean isPdb = true;

    final File file = new File(filename);
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

      // Only consider first model
      //for (int i = model_num - 1; i < size; ++i) {
        final PdbModel model = (isPdb) ? pdbModels.get(0).filteredNewInstance(MoleculeType.RNA) : cifModels.get(0).filteredNewInstance(MoleculeType.RNA);
        double bsgcCounter = 0;
        double rbgcCounter = 0;
        double restCounter = 0;
        Coordinates bsgcAtom = new Coordinates();
        Coordinates rbgcAtom = new Coordinates();
        Coordinates restAtom = new Coordinates();

        String currentResidueName = null;
        String currentAtomKey = null;

        for (final PdbAtomLine atom : model.atoms()) {
          final String atomName = atom.atomName();
          final String key = (StringUtils.isEmpty(atom.insertionCode()))
              ? new StringBuilder(atom.chainIdentifier()).append(atom.residueNumber()).toString()
              : new StringBuilder(atom.chainIdentifier())
                    .append(atom.residueNumber())
                    .append(atom.insertionCode())
                    .toString();

          final String residueName = getOneLetterCode(atom.residueName());
          if (residueName.isBlank()) {
            LOGGER.warn(
                String.format("Found incompatible atom residue name: %s", atom.residueName()));
            continue;
          }
          if (StringUtils.isBlank(currentAtomKey)) {
            currentAtomKey = key;
          }
          if (StringUtils.isBlank(currentResidueName)) {
            currentResidueName = residueName;
          }

          if (!StringUtils.equals(currentAtomKey, key)) {
            if (bsgcCounter > 0) {
              bsgcAtom.x = bsgcAtom.x / bsgcCounter;
              bsgcAtom.y = bsgcAtom.y / bsgcCounter;
              bsgcAtom.z = bsgcAtom.z / bsgcCounter;
            }
            if (rbgcCounter > 0) {
              rbgcAtom.x = rbgcAtom.x / rbgcCounter;
              rbgcAtom.y = rbgcAtom.y / rbgcCounter;
              rbgcAtom.z = rbgcAtom.z / rbgcCounter;
            }
            if (restCounter > 0) {
              restAtom.x = restAtom.x / restCounter;
              restAtom.y = restAtom.y / restCounter;
              restAtom.z = restAtom.z / restCounter;
            }
//            if ((restCounter == rest_atoms.size()) && (rbgcCounter == ryboze_atoms.size())
//                && (bsgcCounter == getBaseAtomsCount(currentResidueName))) {
//          Allow incomplete but must be at least one from each
            if ((restCounter >= 1) && (rbgcCounter >= 1)
                && (bsgcCounter >= 1 )) {
              ArrayList<Coordinates> grained =
                  new ArrayList<Coordinates>(Arrays.asList(bsgcAtom, rbgcAtom, restAtom));
              result.add(new Nucleotide(grained, currentResidueName, currentAtomKey));
            } else {
              if (allowIncomplete && bsgcCounter >= 1. && rbgcCounter >= 1. && restCounter >= 1.) {
                final ArrayList<Coordinates> grained =
                    new ArrayList<Coordinates>(Arrays.asList(bsgcAtom, rbgcAtom, restAtom));
                result.add(new Nucleotide(grained, currentResidueName, currentAtomKey));
                LOGGER.warn(String.format("%s %s.", restCounter, rest_atoms.size()));
                LOGGER.warn(String.format("%s %s.", rbgcCounter, ryboze_atoms.size()));
                LOGGER.warn(
                    String.format("%s %s.", bsgcCounter, getBaseAtomsCount(currentResidueName)));
                LOGGER.warn(String.format("Incomplete residue %s.", currentAtomKey));
              }
            }
            // Reset
            currentAtomKey = key;
            currentResidueName = residueName;
            bsgcAtom = new Coordinates();
            rbgcAtom = new Coordinates();
            restAtom = new Coordinates();
            bsgcCounter = 0;
            rbgcCounter = 0;
            restCounter = 0;
          }

          for (String name : rest_atoms) {
            if (name.equals(atomName)) {
              restAtom.x = restAtom.x + atom.x();
              restAtom.y = restAtom.y + atom.y();
              restAtom.z = restAtom.z + atom.z();
              restCounter++;
              break;
            }
          }
          for (String name : ryboze_atoms) {
            if (name.equals(atomName)) {
              rbgcAtom.x = rbgcAtom.x + atom.x();
              rbgcAtom.y = rbgcAtom.y + atom.y();
              rbgcAtom.z = rbgcAtom.z + atom.z();
              rbgcCounter++;
              break;
            }
          }

          final ArrayList<String> base_atoms;
          if (residueName.equals("U")) {
            base_atoms = U_atoms;
          } else if (residueName.equals("C")) {
            base_atoms = C_atoms;
          } else if (residueName.equals("A")) {
            base_atoms = A_atoms;
          } else { // G
            base_atoms = G_atoms;
          }

          for (String name : base_atoms) {
            if (name.equals(atomName)) {
              bsgcAtom.x = bsgcAtom.x + atom.x();
              bsgcAtom.y = bsgcAtom.y + atom.y();
              bsgcAtom.z = bsgcAtom.z + atom.z();
              bsgcCounter++;
              break;
            }
          }
        }

        // Average and add last residue.
        if (bsgcCounter >= 1.) {
          bsgcAtom.x = bsgcAtom.x / bsgcCounter;
          bsgcAtom.y = bsgcAtom.y / bsgcCounter;
          bsgcAtom.z = bsgcAtom.z / bsgcCounter;
        }
        if (rbgcCounter >= 1.) {
          rbgcAtom.x = rbgcAtom.x / rbgcCounter;
          rbgcAtom.y = rbgcAtom.y / rbgcCounter;
          rbgcAtom.z = rbgcAtom.z / rbgcCounter;
        }
        if (restCounter >= 1.) {
          restAtom.x = restAtom.x / restCounter;
          restAtom.y = restAtom.y / restCounter;
          restAtom.z = restAtom.z / restCounter;
        }
//            if ((restCounter == rest_atoms.size()) && (rbgcCounter == ryboze_atoms.size())
//                && (bsgcCounter == getBaseAtomsCount(currentResidueName))) {
            if ((restCounter >= 1) && (rbgcCounter >= 1)
                && (bsgcCounter >= 1 )) {
          final ArrayList<Coordinates> grained =
              new ArrayList<Coordinates>(Arrays.asList(bsgcAtom, rbgcAtom, restAtom));
          result.add(new Nucleotide(grained, currentResidueName, currentAtomKey));
        } else {
          if (allowIncomplete && bsgcCounter >= 1. && rbgcCounter >= 1. && restCounter >= 1.) {
            final ArrayList<Coordinates> grained =
                new ArrayList<Coordinates>(Arrays.asList(bsgcAtom, rbgcAtom, restAtom));
            result.add(new Nucleotide(grained, currentResidueName, currentAtomKey));
            LOGGER.warn(String.format("%s %s.", restCounter, rest_atoms.size()));
            LOGGER.warn(String.format("%s %s.", rbgcCounter, ryboze_atoms.size()));
            LOGGER.warn(
                String.format("%s %s.", bsgcCounter, getBaseAtomsCount(currentResidueName)));
            LOGGER.warn(String.format("Incomplete residue %s.", currentAtomKey));
          }
        }
      //} // for
      // Only if no exception was thrown
      // It is possible that auto detect did not work properly.
      // We can try and run the parser again with different parse type.
      if (result.size() == 0 && inputType == "auto") {
        LOGGER.warn(
            "Parsed structure without any atoms to analyze. Trying another type pdb/cif instead.");
        return (isPdb) ? StructureToCoarseGrained(filename, "cif", allowIncomplete)
                       : StructureToCoarseGrained(filename, "pdb", allowIncomplete);
      }
    } catch (IOException|IllegalArgumentException e) {
      e.printStackTrace();
    }

    return result;
  }

  private static final int getBaseAtomsCount(final String resName) {
    final ArrayList<String> base_atoms;
    if (resName.equals("U")) {
      base_atoms = U_atoms;
    } else if (resName.equals("C")) {
      base_atoms = C_atoms;
    } else if (resName.equals("A")) {
      base_atoms = A_atoms;
    } else { // G
      base_atoms = G_atoms;
    }
    return base_atoms.size();
  }

  private static String getOneLetterCode(final String resName) {
    if ((StringUtils.equalsIgnoreCase("URI", resName))
        || (StringUtils.equalsIgnoreCase("URA", resName))) {
      return "U";
    } else if (StringUtils.equalsIgnoreCase("ADE", resName)) {
      return "A";
    } else if (StringUtils.equalsIgnoreCase("GUA", resName)) {
      return "G";
    } else if (StringUtils.equalsIgnoreCase("CYT", resName)) {
      return "C";
    } else if (StringUtils.equalsIgnoreCase(resName, "U")) {
      return "U";
    } else if (StringUtils.equalsIgnoreCase(resName, "C")) {
      return "C";
    } else if (StringUtils.equalsIgnoreCase(resName, "G")) {
      return "G";
    } else if (StringUtils.equalsIgnoreCase(resName, "A")) {
      return "A";
    } else if (modified_3to1.containsKey(resName.toUpperCase())) {
      return modified_3to1.get(resName.toUpperCase());
    } else if (StringUtils.equalsIgnoreCase(resName, "DT")) {
      // Convert DT into U (as it should work be similar)
      return "U";
    } else if (StringUtils.equalsIgnoreCase(resName, "DC")) {
      return "C";
    } else if (StringUtils.equalsIgnoreCase(resName, "DG")) {
      return "G";
    } else if (StringUtils.equalsIgnoreCase(resName, "DA")) {
      return "A";
    } else if (StringUtils.endsWithIgnoreCase(resName, "U")) {
      return "U";
    } else if (StringUtils.endsWithIgnoreCase(resName, "C")) {
      return "C";
    } else if (StringUtils.endsWithIgnoreCase(resName, "G")) {
      return "G";
    } else if (StringUtils.endsWithIgnoreCase(resName, "A")) {
      return "A";
    }
    return "";
    // throw new IllegalArgumentException(String.format("Only ADE, GUA, CYT, URI|URA are supported.
    // Provided %s", resName));
  }
}
