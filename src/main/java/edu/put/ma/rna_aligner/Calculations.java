package edu.put.ma.rna_aligner;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.AtomImpl;
import org.biojava.nbio.structure.SVDSuperimposer;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.jama.Matrix;

public class Calculations {
	
	public static Atom[] CopyAtomTable(final Atom[] atoms) {
		Atom[] copy = new Atom[atoms.length];
		for (int i = 0; i < atoms.length; ++i) {
			Atom a = new AtomImpl();
			a.setX(atoms[i].getX());
			a.setY(atoms[i].getY());
			a.setZ(atoms[i].getZ());
			copy[i] = a;
		}
		return copy;
	}
	
	public static SVDSuperimposer FitForRMSD(final Atom[] first_atoms, final Atom[] second_atoms) {
    	SVDSuperimposer superimposer = null;
    	try {
    		superimposer = new SVDSuperimposer(first_atoms, second_atoms);
    	} catch (StructureException e) {
    		//
    	}
    	return superimposer;
	}
	
	public static void rotate(Atom atom, Matrix m){
		double x = atom.getX();
		double y = atom.getY() ;
		double z = atom.getZ();
		final double[][] mt = m.getArray();
		atom.setX(x * mt[0][0] + y * mt[1][0] + z * mt[2][0]);
		atom.setY(x * mt[0][1] + y * mt[1][1] + z * mt[2][1]);
		atom.setZ(x * mt[0][2] + y * mt[1][2] + z * mt[2][2]);
	}
	
	public static void shift(Atom a, Atom b){
		a.setX(a.getX() + b.getX());
		a.setY(a.getY() + b.getY());
		a.setZ(a.getZ() + b.getZ());
	}
	
	public static void MoveAtomsInplace(Atom[] first, Atom[] second) {
		SVDSuperimposer superimposer = FitForRMSD(first, second);
        for (Atom atom : second) {
        	rotate(atom, superimposer.getRotation());
            shift(atom, superimposer.getTranslation());
        }
	}
	
	public static Atom[] MoveAtomsForRMSD(final Atom[] atoms, final SVDSuperimposer superimposer) {
    	Atom[] atoms_moved_copy = CopyAtomTable(atoms);

        for (Atom atom : atoms_moved_copy) {
        	rotate(atom, superimposer.getRotation());
            shift(atom, superimposer.getTranslation());
        }
        return atoms_moved_copy;
	}
	
    public static double FitAndCalculateRMSD(final Atom[] first_atoms, final Atom[] second_atoms) {
    	if (first_atoms.length != second_atoms.length) {
    		System.out.println(String.format("WRONG LENGTH: %d %d", first_atoms.length, second_atoms.length));
    		return -1;
    	}

    	SVDSuperimposer superimposer = FitForRMSD(first_atoms, second_atoms);

    	final Atom[] fit_second_atoms = MoveAtomsForRMSD(second_atoms, superimposer);

		return CalculateRMSD(first_atoms, fit_second_atoms);
    }
    
    public static double CalculateRMSD(final Atom[] first_atoms, final Atom[] second_atoms) {
    	try {
			return SVDSuperimposer.getRMS(first_atoms, second_atoms);
		} catch (StructureException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
    	return 0;
    }
    
    public static double CalculateRMSD(final Atom[] first_atoms, final Atom[] second_atoms, final SVDSuperimposer superimposer) {
    	final Atom[] fit_second_atoms = MoveAtomsForRMSD(second_atoms, superimposer);

    	try {
			return SVDSuperimposer.getRMS(first_atoms, fit_second_atoms);
		} catch (StructureException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
    	return 0;
    }
}