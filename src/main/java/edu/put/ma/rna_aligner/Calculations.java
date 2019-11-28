package edu.put.ma.rna_aligner;

import java.util.ArrayList;

import org.biojava.nbio.structure.jama.Matrix;

public class Calculations {
	
	public static ArrayList<Coordinates> CopyAtoms(final ArrayList<Coordinates> atoms) {
		ArrayList<Coordinates> copy = new ArrayList<Coordinates>(atoms.size());
		for (int i = 0; i < atoms.size(); ++i) {
			copy.add(atoms.get(i).clone());
		}
		return copy;
	}
	
	public static Superimposer FitForRMSD(final ArrayList<Coordinates> first_atoms, final ArrayList<Coordinates> second_atoms) {
    	return new Superimposer(first_atoms, second_atoms);
	}
	
	public static void rotate(Coordinates atom, final Matrix m) {
		double x = atom.x;
		double y = atom.y;
		double z = atom.z;
		final double[][] mt = m.getArray();
		atom.x = (x * mt[0][0] + y * mt[1][0] + z * mt[2][0]);
		atom.y = (x * mt[0][1] + y * mt[1][1] + z * mt[2][1]);
		atom.z = (x * mt[0][2] + y * mt[1][2] + z * mt[2][2]);
	}
	
	public static void shift(Coordinates atom, final Coordinates shift) {
		atom.x += shift.x;
		atom.y += shift.y;
		atom.z += shift.z;
	}
	
	public void MoveAtomsInplace(final ArrayList<Coordinates> reference, ArrayList<Coordinates> target) {
		Superimposer superimposer = FitForRMSD(reference, target);
        for (Coordinates atom : target) {
        	rotate(atom, superimposer.getRotation());
            shift(atom, superimposer.getTranslation());
        }
	}
	
	public static ArrayList<Coordinates> MoveAtomsForRMSD(final ArrayList<Coordinates> atoms, final Superimposer superimposer) {
		ArrayList<Coordinates> atoms_moved_copy = CopyAtoms(atoms);

        for (Coordinates atom : atoms_moved_copy) {
        	rotate(atom, superimposer.getRotation());
            shift(atom, superimposer.getTranslation());
        }
        return atoms_moved_copy;
	}
	
    public static double FitAndCalculateRMSD(final ArrayList<Coordinates> first_atoms, final ArrayList<Coordinates> second_atoms) {
    	if (first_atoms.size() != second_atoms.size()) {
    		System.out.println(String.format("WRONG LENGTH: %d %d", first_atoms.size(), second_atoms.size()));
    		return -1;
    	}

    	Superimposer superimposer = FitForRMSD(first_atoms, second_atoms);

    	final ArrayList<Coordinates> fit_second_atoms = MoveAtomsForRMSD(second_atoms, superimposer);

		return CalculateRMSD(first_atoms, fit_second_atoms);
    }
    
    public static double CalculateRMSD(final ArrayList<Coordinates> first_atoms, final ArrayList<Coordinates> second_atoms) {
		return Superimposer.getRMS(first_atoms, second_atoms);
    }
    
    public static double CalculateRMSD(final ArrayList<Coordinates> first_atoms, final ArrayList<Coordinates> second_atoms, final Superimposer superimposer) {
    	final ArrayList<Coordinates> fit_second_atoms = MoveAtomsForRMSD(second_atoms, superimposer);
		return Superimposer.getRMS(first_atoms, fit_second_atoms);
    }
}