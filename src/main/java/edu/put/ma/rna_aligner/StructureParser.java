package edu.put.ma.rna_aligner;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.StringUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import pl.poznan.put.pdb.PdbAtomLine;
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

  
  public ArrayList<Nucleotide> StructureToCoarseGrained(final String filename, final String inputType) {
	ArrayList<Nucleotide> result = new ArrayList<Nucleotide>();  
	boolean isPdb = true;
	
    final File file =  new File(filename);
    try {
		final String structureContent = FileUtils.readFileToString(file, StandardCharsets.UTF_8);
		
		if (inputType == "auto") {
			// Try Detect if the file is mmCIF or PDB.
			// Search for _loop in the file. If it is present we can assume it is CIF format.
			if (structureContent.indexOf("_loop") != -1) {
				isPdb = false;
			} // else remains isPdb = true;
		} else if (inputType == "cif") {
			isPdb = false;
		} // else remains isPdb = true;

		
		    // parse the data
		    final PdbParser pdbParser = new PdbParser(false);
		    final List<PdbModel> pdbModels = (isPdb) ? pdbParser.parse(structureContent) : null;
		    final List<CifModel> cifModels = (!isPdb) ? CifParser.parse(structureContent) : null;
		    
		    //final int size = (isPdb) ? pdbModels.size() : cifModels.size();
		    // Only consider first model
		    final int size = 1;
		    
		    for (int i = 0; i < size; ++i) {
		    	final PdbModel model = (isPdb) ? pdbModels.get(i) : cifModels.get(i);
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
		                    : new StringBuilder(atom.chainIdentifier()).append(atom.residueNumber()).append(atom.insertionCode()).toString();

		    	    
		    		final String residueName = getOneLetterCode(atom.residueName());
		    		if (residueName.isBlank()) {
		    			LOGGER.warn(String.format("Found incompatible atom residue name: %s", atom.residueName()));
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
		                if ((restCounter == rest_atoms.size()) && (rbgcCounter == ryboze_atoms.size())
		                    && (bsgcCounter == getBaseAtomsCount(currentResidueName))) {
		                  ArrayList<Coordinates> grained =
		                      new ArrayList<Coordinates>(Arrays.asList(bsgcAtom, rbgcAtom, restAtom));
		                  result.add(new Nucleotide(grained, currentResidueName, currentAtomKey));
		                } else {
		  		          if (bsgcCounter >= 1. && rbgcCounter >= 1. && restCounter >= 1.) {
		  			          final ArrayList<Coordinates> grained =
		  			              new ArrayList<Coordinates>(Arrays.asList(bsgcAtom, rbgcAtom, restAtom));
		  			          result.add(new Nucleotide(grained, currentResidueName, currentAtomKey));
		  			          }
		                  LOGGER.warn(String.format("%s %s.", restCounter, rest_atoms.size()));
		                  LOGGER.warn(String.format("%s %s.", rbgcCounter, ryboze_atoms.size()));
		                  LOGGER.warn(String.format("%s %s.", bsgcCounter, getBaseAtomsCount(currentResidueName)));
		                  LOGGER.warn(String.format("Incomplete residue %s.", currentAtomKey));
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
		        if ((restCounter == rest_atoms.size()) && (rbgcCounter == ryboze_atoms.size())
		            && (bsgcCounter == getBaseAtomsCount(currentResidueName))) {
		          final ArrayList<Coordinates> grained =
		              new ArrayList<Coordinates>(Arrays.asList(bsgcAtom, rbgcAtom, restAtom));
		          result.add(new Nucleotide(grained, currentResidueName, currentAtomKey));
		        } else {
		          if (bsgcCounter >= 1. && rbgcCounter >= 1. && restCounter >= 1.) {
		          final ArrayList<Coordinates> grained =
		              new ArrayList<Coordinates>(Arrays.asList(bsgcAtom, rbgcAtom, restAtom));
		          result.add(new Nucleotide(grained, currentResidueName, currentAtomKey));
		          }
		          LOGGER.warn(String.format("%s %s.", restCounter, rest_atoms.size()));
		          LOGGER.warn(String.format("%s %s.", rbgcCounter, ryboze_atoms.size()));
		          LOGGER.warn(String.format("%s %s.", bsgcCounter, getBaseAtomsCount(currentResidueName)));
		          LOGGER.warn(String.format("Incomplete residue %s.", currentAtomKey));
		        }
		    }
	} catch (IOException e) {
		e.printStackTrace();
	}
    // It is possible that auto detect did not work properly.
    // We can try and run the parser again with different parse type.
    if (result.size() == 0 && inputType == "auto") {
      LOGGER.warn("Parsed structure without any atoms to analyze. Trying another type pdb/cif instead.");
      return (isPdb) ? StructureToCoarseGrained(filename, "cif") : StructureToCoarseGrained(filename, "pdb");
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
    //throw new IllegalArgumentException(String.format("Only ADE, GUA, CYT, URI|URA are supported. Provided %s", resName));
  }
}
