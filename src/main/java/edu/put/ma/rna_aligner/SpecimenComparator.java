package edu.put.ma.rna_aligner;

import java.util.Comparator;

public class SpecimenComparator implements Comparator {
	private double rmsd_limit = 3.5;
	
	SpecimenComparator(double _rmsd_limit) {
		this.rmsd_limit = _rmsd_limit;
	}
	
	public int compare(Object _lhs, Object _rhs) {
		Specimen lhs = (Specimen) _lhs;
		Specimen rhs = (Specimen) _rhs;
		if (lhs.getUsedNucleotidesNumber() == 0 && rhs.getUsedNucleotidesNumber() == 0) {
			return 0;
		}
		if (lhs.calculateRMSD() < this.rmsd_limit && rhs.calculateRMSD() >= this.rmsd_limit) {
			return -1;
		} else if (lhs.calculateRMSD() >= this.rmsd_limit && rhs.calculateRMSD() < this.rmsd_limit) {
			return 1;
		}
		if (lhs.calculateRMSD() < this.rmsd_limit && rhs.calculateRMSD() < this.rmsd_limit) {
			if (lhs.getRemainingNucleotidesNumber() < rhs.getRemainingNucleotidesNumber()) {
				return -1;
			} else if (lhs.getRemainingNucleotidesNumber() > rhs.getRemainingNucleotidesNumber()) {
				return 1;
			} else {
				if (lhs.calculateRMSD() < rhs.calculateRMSD()) {
					return -1;
				} else if (lhs.calculateRMSD() > rhs.calculateRMSD()) {
					return 1;
				} else {
					return 0;
				}
			}
		} else {
			if (lhs.calculateRMSD() > rhs.calculateRMSD()) {
				return 1;
			} else if (lhs.calculateRMSD() < rhs.calculateRMSD()) {
				return -1;
			} else {
				if (lhs.getRemainingNucleotidesNumber() < rhs.getRemainingNucleotidesNumber()) {
					return -1;
				} else if (lhs.getRemainingNucleotidesNumber() > rhs.getRemainingNucleotidesNumber()) {
					return 1;
				} else {
					return 0;
				}
			}
		}
	}
}
