package edu.put.ma.rna_aligner;

import java.util.ArrayList;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Calc;

public class Dist {
	public ArrayList<Double> distances = new ArrayList<Double>();
	
	Dist(final ArrayList<Double> _distances) {
		this.distances.clear();
		this.distances.addAll(_distances);
	}
	
	Dist(final Nucleotide lhs, final Nucleotide rhs) {
		for (int i = 0; i < lhs.representatives.size(); i++) {
			distances.add(CalculateDistance(lhs.representatives.get(i),
					                        rhs.representatives.get(i)));
		}
	}
	
	private double CalculateDistance(final Atom lhs, final Atom rhs) {
		return Calc.getDistance(lhs, rhs);
	}
	
	// Calculate mini best case scenario RMSD only on distances.
	public static double Similarity(final Dist lhs, final Dist rhs) {
		double similarity = 0;
		for (int i = 0; i < lhs.distances.size(); i++) {
			similarity += Math.pow((lhs.distances.get(i) - rhs.distances.get(i)), 2);
		}
		return similarity;
	}
	
	// Calculate mini RMSD only on distances of 3 nucleotides.
	public static double Similarity(final Dist lhs, final Dist rhs, final Dist lhs2, final Dist rhs2, final Dist lhs3, final Dist rhs3) {
		// TODO
		return 0.0;
	}
}
