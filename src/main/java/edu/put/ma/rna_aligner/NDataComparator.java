package edu.put.ma.rna_aligner;

import java.util.Comparator;

public class NDataComparator implements Comparator<NData> {
  @Override
  public int compare(NData lhs, NData rhs) {
    if (lhs.rmsd < rhs.rmsd) {
      return 1;
    } else if (lhs.rmsd > rhs.rmsd) {
      return -1;
    } 
    return 0;
  }
}