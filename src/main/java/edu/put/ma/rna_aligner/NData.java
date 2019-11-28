package edu.put.ma.rna_aligner;

public class NData {
	public int index1;
	public int index2;
	public double rmsd;
	public Superimposer superimposer;
	
	public NData(int _index1, int _index2) {
		this.index1 = _index1;
		this.index2 = _index2;
	}
	
	public NData(int _index1, int _index2, double _rmsd) {
		this.index1 = _index1;
		this.index2 = _index2;
		this.rmsd = _rmsd;
	}
	
	public NData(int _index1, int _index2, double _rmsd, Superimposer _superimposer) {
		this.index1 = _index1;
		this.index2 = _index2;
		this.rmsd = _rmsd;
		this.superimposer = _superimposer;
	}

	public NData(NData rhs) {
		this.index1 = rhs.index1;
		this.index2 = rhs.index2;
	}
	
	@Override
	public Object clone() {
		return new NData(this.index1, this.index2, this.rmsd, this.superimposer);
	}
	
}
