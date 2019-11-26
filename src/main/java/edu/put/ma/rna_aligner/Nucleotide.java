package edu.put.ma.rna_aligner;

import java.util.ArrayList;
import org.biojava.nbio.structure.Atom;


public class Nucleotide {
	public ArrayList<Atom> representatives;
	
	public Nucleotide() {
		super();
		representatives = new ArrayList<Atom>();
	}
	
	public Nucleotide(final ArrayList<Atom> _representatives) {
		representatives = _representatives;
	}
	
	
	@Override
	public Object clone() {
		Nucleotide clone = new Nucleotide();
		for (final Atom rep : this.representatives) {
			clone.representatives.add((Atom) rep.clone());
		}
		return clone;
	}
	
	
	public static ArrayList<Atom> NucleotidesToList(final ArrayList<Nucleotide> nucleotides) {
		ArrayList<Atom> result = new ArrayList<Atom>();
		for (Nucleotide nucleotide : nucleotides) {
			for (final Atom rep : nucleotide.representatives) {
				result.add(rep);
			}		
		}
		return result;
	}
	
	public static ArrayList<Atom> NucleotidesToListMapped(final ArrayList<Nucleotide> nucleotides, final int[] map, boolean primary) {
		ArrayList<Atom> result = new ArrayList<Atom>();
		if (primary) {
			// If map[i] == 1, take it.
			for (int i = 0; i < map.length; i++) {
				if (map[i] == 1) {
					for (final Atom rep : nucleotides.get(i).representatives) {
						result.add(rep);
					}
				}
			}
		} else {
			// Take it in the order stated in the map
			for (int i = 0; i < map.length; i++) {
				if (map[i] >= 0) {
					for (final Atom rep : nucleotides.get(map[i]).representatives) {
						result.add(rep);
					}
				}
			}
		}
		return result;
	}

	
	public static Atom[] NucleotidesToTable(final ArrayList<Nucleotide> nucleotides) {
		final ArrayList<Atom> atoms_arr = NucleotidesToList(nucleotides);
		return atoms_arr.toArray(new Atom[atoms_arr.size()]);
	}
	
	public static Atom[] NucleotidesToMappedTable(final ArrayList<Nucleotide> nucleotides, final int[] map, boolean primary) {
		final ArrayList<Atom> atoms_arr = NucleotidesToListMapped(nucleotides, map, primary);
		return atoms_arr.toArray(new Atom[atoms_arr.size()]);
	}
}
