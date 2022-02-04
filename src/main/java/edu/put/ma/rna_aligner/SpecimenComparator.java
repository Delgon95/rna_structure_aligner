package edu.put.ma.rna_aligner;

import java.util.Comparator;

public class SpecimenComparator implements Comparator<Object> {
  private double rmsd_limit = 3.5;

  SpecimenComparator(final double _rmsd_limit) {
    this.rmsd_limit = _rmsd_limit;
  }

  public int compare(final Object _lhs, final Object _rhs) {
    final Specimen lhs = (Specimen) _lhs;
    final Specimen rhs = (Specimen) _rhs;

    if (Double.compare(lhs.calculateRMSD(), rmsd_limit) <= 0
        && Double.compare(rhs.calculateRMSD(), rmsd_limit) <= 0) {
      if (lhs.getIncorrectlyAlignedResiduesRatio() < rhs.getIncorrectlyAlignedResiduesRatio())
        return -1;
      else if (lhs.getIncorrectlyAlignedResiduesRatio() > rhs.getIncorrectlyAlignedResiduesRatio())
        return 1;
      else {
        if (lhs.getUsedNucleotidesNumber() < rhs.getUsedNucleotidesNumber())
          return 1;
        else if (lhs.getUsedNucleotidesNumber() > rhs.getUsedNucleotidesNumber())
          return -1;
        else {
          if (lhs.calculateRMSD() < rhs.calculateRMSD())
            return -1;
          else if (lhs.calculateRMSD() > rhs.calculateRMSD())
            return 1;
          else
            return 0;
        }
      }
    } else {
      if (Double.compare(lhs.calculateRMSD(), rmsd_limit) <= 0 && rhs.calculateRMSD() > rmsd_limit)
        return -1;
      else if (lhs.calculateRMSD() > rmsd_limit
          && Double.compare(rhs.calculateRMSD(), rmsd_limit) <= 0)
        return 1;
      else {
        if (lhs.calculateRMSD() > rhs.calculateRMSD())
          return 1;
        else if (lhs.calculateRMSD() < rhs.calculateRMSD())
          return -1;
        else {
          if (lhs.getIncorrectlyAlignedResiduesRatio() < rhs.getIncorrectlyAlignedResiduesRatio())
            return -1;
          else if (lhs.getIncorrectlyAlignedResiduesRatio()
              > rhs.getIncorrectlyAlignedResiduesRatio())
            return 1;
          else {
            if (lhs.getUsedNucleotidesNumber() < rhs.getUsedNucleotidesNumber())
              return 1;
            else if (lhs.getUsedNucleotidesNumber() > rhs.getUsedNucleotidesNumber())
              return -1;
            else
              return 0;
          }
        }
      }
    }
  }
}
