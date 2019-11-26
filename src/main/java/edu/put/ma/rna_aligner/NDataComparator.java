package edu.put.ma.rna_aligner;

import java.util.Comparator;

public class NDataComparator implements Comparator<Object> {
	@Override
	public int compare(Object _lhs, Object _rhs) {
		NData lhs = (NData) _lhs;
		NData rhs = (NData) _rhs;
		
		if (lhs.rmsd < rhs.rmsd) {
			return 1;
		} else if (lhs.rmsd > rhs.rmsd) {
			return -1;
		} 
		return 0;
	}
}
