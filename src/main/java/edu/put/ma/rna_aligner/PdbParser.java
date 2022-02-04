package edu.put.ma.rna_aligner;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import org.apache.commons.lang3.StringUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class PdbParser {
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

  @SuppressWarnings("resource")
  public ArrayList<Nucleotide> ParsePdbToCoarseGrained(final String pdb_filename) {
    ArrayList<Nucleotide> result = new ArrayList<Nucleotide>();

    String line;
    try {
      FileReader fileReader = new FileReader(pdb_filename);
      BufferedReader bufferedReader = new BufferedReader(fileReader);
      String currentKey = null, currentResName = null;
      double bsgcCounter = 0;
      double rbgcCounter = 0;
      double restCounter = 0;
      Coordinates bsgcAtom = new Coordinates();
      Coordinates rbgcAtom = new Coordinates();
      Coordinates restAtom = new Coordinates();
      Coordinates atom = new Coordinates();
      ArrayList<Coordinates> atoms = new ArrayList<Coordinates>();
      while ((line = bufferedReader.readLine()) != null) {
        double x, y, z;
        int residue;
        String chain, iCode, resName, key;
        String[] split_line = line.split("\\s+");
        if (!split_line[0].equals("ATOM")) {
          continue;
        }
        if (split_line.length != 12) {
          x = Double.parseDouble(line.substring(31, 38).trim());
          y = Double.parseDouble(line.substring(38, 46).trim());
          z = Double.parseDouble(line.substring(46, 54).trim());
          chain = line.substring(21, 22).trim();
          residue = Integer.parseInt(line.substring(23, 26).trim());
          iCode = line.substring(26, 27).trim();
          resName = line.substring(17, 20).trim();
        } else {
          x = Double.parseDouble(split_line[6]);
          y = Double.parseDouble(split_line[7]);
          z = Double.parseDouble(split_line[8]);
          atom.x = x;
          atom.y = y;
          atom.z = z;
          atoms.add(atom);
          if (StringUtils.isNumeric(split_line[5])) {
            residue = Integer.parseInt(split_line[5]);
            iCode = "";
          } else {
            int residueSize = StringUtils.length(split_line[5]);
            residue = Integer.parseInt(StringUtils.substring(split_line[5], 0, residueSize - 1));
            iCode = StringUtils.substring(split_line[5], residueSize - 1, residueSize);
          }
          chain = split_line[4];
          resName = split_line[3];
        }
        if (StringUtils.isEmpty(chain)) {
          chain = "A";
        }
        key = (StringUtils.isEmpty(iCode))
            ? new StringBuilder(chain).append(String.valueOf(residue)).toString()
            : new StringBuilder(chain).append(String.valueOf(residue)).append(iCode).toString();
        if (StringUtils.isBlank(currentKey)) {
          currentKey = key;
        }
        if (StringUtils.length(resName) == 3) {
          resName = getOneLetterCode(resName);
        }
        if (StringUtils.isBlank(currentResName)) {
          currentResName = resName;
        }

        atom = new Coordinates();
        /*BSGC, RBGC, C5'
         * Field: 0         1  2     3 4   5     6       7       8        9     10             11
         * Line:  ATOM      1  P     U A   1     -12.818 -47.368 -15.305  1.00  0.00           P
         *
         * Fields:
         * 1 - Atom number
         * 2 - Atom type
         * 3 - Type (U, G, A, C)
         * 5 - Residue number
         * 6 - X
         * 7 - Y
         * 8 - Z
         */

        if (!StringUtils.equals(currentKey, key)) {
          bsgcAtom.x = bsgcAtom.x / bsgcCounter;
          bsgcAtom.y = bsgcAtom.y / bsgcCounter;
          bsgcAtom.z = bsgcAtom.z / bsgcCounter;
          rbgcAtom.x = rbgcAtom.x / rbgcCounter;
          rbgcAtom.y = rbgcAtom.y / rbgcCounter;
          rbgcAtom.z = rbgcAtom.z / rbgcCounter;
          restAtom.x = restAtom.x / restCounter;
          restAtom.y = restAtom.y / restCounter;
          restAtom.z = restAtom.z / restCounter;
          if ((restCounter == rest_atoms.size()) && (rbgcCounter == ryboze_atoms.size())
              && (bsgcCounter == getBaseAtomsCount(currentResName))) {
            ArrayList<Coordinates> grained =
                new ArrayList<Coordinates>(Arrays.asList(bsgcAtom, rbgcAtom, restAtom));
            result.add(new Nucleotide(grained, currentResName, currentKey));
          } else {
            ArrayList<Coordinates> grained =
                new ArrayList<Coordinates>(Arrays.asList(bsgcAtom, rbgcAtom, restAtom));
            result.add(new Nucleotide(grained, currentResName, currentKey));
            LOGGER.warn(String.format("%s %s.", restCounter, rest_atoms.size()));
            LOGGER.warn(String.format("%s %s.", rbgcCounter, ryboze_atoms.size()));
            LOGGER.warn(String.format("%s %s.", bsgcCounter, getBaseAtomsCount(currentResName)));
            LOGGER.warn(String.format("Incomplete residue %s.", currentKey));
          }
          // Reset
          currentKey = key;
          currentResName = resName;
          bsgcAtom = new Coordinates();
          rbgcAtom = new Coordinates();
          restAtom = new Coordinates();
          atoms = new ArrayList<Coordinates>();
          bsgcCounter = 0;
          rbgcCounter = 0;
          restCounter = 0;
        }
        for (String name : rest_atoms) {
          if (name.equals(split_line[2])) {
            restAtom.x = restAtom.x + x;
            restAtom.y = restAtom.y + y;
            restAtom.z = restAtom.z + z;
            restCounter++;
            break;
          }
        }
        for (String name : ryboze_atoms) {
          if (name.equals(split_line[2])) {
            rbgcAtom.x = rbgcAtom.x + x;
            rbgcAtom.y = rbgcAtom.y + y;
            rbgcAtom.z = rbgcAtom.z + z;
            rbgcCounter++;
            break;
          }
        }

        final ArrayList<String> base_atoms;
        if (split_line[3].equals("U")) {
          base_atoms = U_atoms;
        } else if (split_line[3].equals("C")) {
          base_atoms = C_atoms;
        } else if (split_line[3].equals("A")) {
          base_atoms = A_atoms;
        } else { // G
          base_atoms = G_atoms;
        }

        for (String name : base_atoms) {
          if (name.equals(split_line[2])) {
            bsgcAtom.x = bsgcAtom.x + x;
            bsgcAtom.y = bsgcAtom.y + y;
            bsgcAtom.z = bsgcAtom.z + z;
            bsgcCounter++;
            break;
          }
        }
      }
      // Average and add last residue.

      bsgcAtom.x = bsgcAtom.x / bsgcCounter;
      bsgcAtom.y = bsgcAtom.y / bsgcCounter;
      bsgcAtom.z = bsgcAtom.z / bsgcCounter;
      rbgcAtom.x = rbgcAtom.x / rbgcCounter;
      rbgcAtom.y = rbgcAtom.y / rbgcCounter;
      rbgcAtom.z = rbgcAtom.z / rbgcCounter;
      restAtom.x = restAtom.x / restCounter;
      restAtom.y = restAtom.y / restCounter;
      restAtom.z = restAtom.z / restCounter;
      if ((restCounter == rest_atoms.size()) && (rbgcCounter == ryboze_atoms.size())
          && (bsgcCounter == getBaseAtomsCount(currentResName))) {
        final ArrayList<Coordinates> grained =
            new ArrayList<Coordinates>(Arrays.asList(bsgcAtom, rbgcAtom, restAtom));
        result.add(new Nucleotide(grained, currentResName, currentKey));
      } else {
        final ArrayList<Coordinates> grained =
            new ArrayList<Coordinates>(Arrays.asList(bsgcAtom, rbgcAtom, restAtom));
        result.add(new Nucleotide(grained, currentResName, currentKey));
        LOGGER.warn(String.format("%s %s.", restCounter, rest_atoms.size()));
        LOGGER.warn(String.format("%s %s.", rbgcCounter, ryboze_atoms.size()));
        LOGGER.warn(String.format("%s %s.", bsgcCounter, getBaseAtomsCount(currentResName)));
        LOGGER.warn(String.format("Incomplete residue %s.", currentKey));
      }

    } catch (IOException e) {
      LOGGER.error(e.getMessage(), e);
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
        || (StringUtils.equalsIgnoreCase("URA", resName)))
      return "U";
    else if (StringUtils.equalsIgnoreCase("ADE", resName))
      return "A";
    else if (StringUtils.equalsIgnoreCase("GUA", resName))
      return "G";
    else if (StringUtils.equalsIgnoreCase("CYT", resName))
      return "C";
    else if (StringUtils.endsWithIgnoreCase(resName, "U"))
      return "U";
    else if (StringUtils.endsWithIgnoreCase(resName, "C"))
      return "C";
    else if (StringUtils.endsWithIgnoreCase(resName, "G"))
      return "G";
    else if (StringUtils.endsWithIgnoreCase(resName, "A"))
      return "A";
    throw new IllegalArgumentException("Only ADE, GUA, CYT, URI|URA are supported.");
  }
}
