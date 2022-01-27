package edu.put.ma.rna_aligner;

import java.util.ArrayList;

public class AlignerOutput {
	// How many nucleotides were used for resulting alignment.
	public final int aligned;
	// Nucleotide index that was used in reference structure.
	public final ArrayList<Integer> referenceIndexes;
	// To what nucleotide index was it mapped in target structure.
	public final ArrayList<Integer> targetMapping;
	// Superimposer that contains shift and rotation matrices.
	public final Superimposer superimposer;
	public final long processingTime;
	public final double rmsd;
	
	AlignerOutput(int _aligned, final ArrayList<Integer> _referenceIndexes, final ArrayList<Integer> _targetMapping,
			final Superimposer _superimposer, final long _processingTime, final double _rmsd) {
		super();
		aligned = _aligned;
		referenceIndexes = _referenceIndexes;
		targetMapping = _targetMapping;
		superimposer = _superimposer;
		processingTime = _processingTime;
		rmsd = _rmsd;
	}
}