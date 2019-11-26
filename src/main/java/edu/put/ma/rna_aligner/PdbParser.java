package edu.put.ma.rna_aligner;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

import edu.put.ma.rna_aligner.Nucleotide;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.AtomImpl;

public class PdbParser {
	private static final ArrayList<String> rest_atoms = new ArrayList<String>(Arrays.asList("C5'"));
	private static final ArrayList<String> ryboze_atoms = new ArrayList<String>(Arrays.asList("C2'", "C1'", "O4'"));
	private static final ArrayList<String> U_atoms = new ArrayList<String>(Arrays.asList("N1", "C6", "C5", "C4", "O4", "N3", "C2", "O2"));
	private static final ArrayList<String> C_atoms = new ArrayList<String>(Arrays.asList("N1", "C6", "C5", "C4", "N4", "N3", "C2", "O2"));
	private static final ArrayList<String> A_atoms = new ArrayList<String>(Arrays.asList("N9", "C8", "N7", "C5", "C6", "N6", "N1", "C2", "N3", "C4"));
	private static final ArrayList<String> G_atoms = new ArrayList<String>(Arrays.asList("N9", "C8", "N7", "C5", "C6", "O6", "N1", "C2", "N2", "N3", "C4"));
	
	@SuppressWarnings("resource")
	public ArrayList<Nucleotide> ParsePdbToCoarseGrained(final String pdb_filename) {
		ArrayList<Nucleotide> result = new ArrayList<Nucleotide>();
		
		String line;
		try {
			FileReader fileReader = new FileReader(pdb_filename);
			BufferedReader bufferedReader = new BufferedReader(fileReader);
			int current_residue = 1;
			
			double bsgcCounter = 0;
			double rbgcCounter = 0;
			double restCounter = 0;
			Atom bsgcAtom = new AtomImpl();
			Atom rbgcAtom = new AtomImpl();
			Atom restAtom = new AtomImpl();
			Atom atom = new AtomImpl();
			ArrayList<Atom> atoms = new ArrayList<Atom>();
			
			while ((line = bufferedReader.readLine()) != null) {
				double x, y, z;
				int residue;
				String[] split_line = line.split("\\s+");
				if (!split_line[0].equals("ATOM")) {
					continue;
				}
				if (split_line.length != 12) {
					x = Double.parseDouble(line.substring(31, 38).trim());
					y = Double.parseDouble(line.substring(39, 46).trim());
					z = Double.parseDouble(line.substring(47, 54).trim());
					residue = Integer.parseInt(line.substring(23, 26).trim());
				} else {
					atom.setX(Double.parseDouble(split_line[6]));
					atom.setY(Double.parseDouble(split_line[7]));
					atom.setZ(Double.parseDouble(split_line[8]));
					atoms.add(atom);
					x = Double.parseDouble(split_line[6]);
					y = Double.parseDouble(split_line[7]);
					z = Double.parseDouble(split_line[8]);
					residue = Integer.parseInt(split_line[5]);
				}

				atom = new AtomImpl();
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
				
				if (current_residue != residue) {
					bsgcAtom.setX(bsgcAtom.getX() / bsgcCounter);
					bsgcAtom.setY(bsgcAtom.getY() / bsgcCounter);
					bsgcAtom.setZ(bsgcAtom.getZ() / bsgcCounter);
					rbgcAtom.setX(rbgcAtom.getX() / rbgcCounter);
					rbgcAtom.setY(rbgcAtom.getY() / rbgcCounter);
					rbgcAtom.setZ(rbgcAtom.getZ() / rbgcCounter);
					restAtom.setX(restAtom.getX() / restCounter);
					restAtom.setY(restAtom.getY() / restCounter);
					restAtom.setZ(restAtom.getZ() / restCounter);
					ArrayList<Atom> grained = new ArrayList<Atom>(Arrays.asList(bsgcAtom, rbgcAtom, restAtom));
					result.add(new Nucleotide(grained));
					// Reset
					current_residue = residue;
					bsgcAtom = new AtomImpl();
					rbgcAtom = new AtomImpl();
					restAtom = new AtomImpl();
					atoms = new ArrayList<Atom>();
					bsgcCounter = 0;
					rbgcCounter = 0;
					restCounter = 0;
				}
				for (String name : rest_atoms) {
					if (name.equals(split_line[2])) {
						restAtom.setX(restAtom.getX() + x);
						restAtom.setY(restAtom.getY() + y);
						restAtom.setZ(restAtom.getZ() + z);
						restCounter++;
						break;
					}
				}
				for (String name : ryboze_atoms) {
					if (name.equals(split_line[2])) {
						rbgcAtom.setX(rbgcAtom.getX() + x);
						rbgcAtom.setY(rbgcAtom.getY() + y);
						rbgcAtom.setZ(rbgcAtom.getZ() + z);
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
						bsgcAtom.setX(bsgcAtom.getX() + x);
						bsgcAtom.setY(bsgcAtom.getY() + y);
						bsgcAtom.setZ(bsgcAtom.getZ() + z);
						bsgcCounter++;
						break;
					}
				}
			}
			// Average and add last residue.
			
			bsgcAtom.setX(bsgcAtom.getX() / bsgcCounter);
			bsgcAtom.setY(bsgcAtom.getX() / bsgcCounter);
			bsgcAtom.setZ(bsgcAtom.getX() / bsgcCounter);
			rbgcAtom.setX(rbgcAtom.getX() / rbgcCounter);
			rbgcAtom.setY(rbgcAtom.getX() / rbgcCounter);
			rbgcAtom.setZ(rbgcAtom.getX() / rbgcCounter);
			restAtom.setX(restAtom.getX() / restCounter);
			restAtom.setY(restAtom.getX() / restCounter);
			restAtom.setZ(restAtom.getX() / restCounter);
			ArrayList<Atom> grained = new ArrayList<Atom>(Arrays.asList(bsgcAtom, rbgcAtom, restAtom));
			result.add(new Nucleotide(grained));
		
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		return result;
	}
}
