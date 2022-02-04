package edu.put.ma.rna_aligner;

import java.text.DecimalFormat;
import java.util.ArrayList;

public class Dist {
  public ArrayList<Double> distances = new ArrayList<Double>();

  Dist(final ArrayList<Double> _distances) {
    this.distances.clear();
    this.distances.addAll(_distances);
  }

  Dist(final Nucleotide lhs, final Nucleotide rhs) {
    for (int i = 0; i < lhs.representatives.size(); i++) {
      distances.add(CalculateDistance(lhs.representatives.get(i), rhs.representatives.get(i)));
    }
  }

  private double CalculateDistance(final Coordinates lhs, final Coordinates rhs) {
    return Coordinates.getDistance(lhs, rhs);
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
  public static double Similarity(final Dist lhs1, final Dist rhs1, final Dist lhs2,
      final Dist rhs2, final Dist lhs3, final Dist rhs3) {
    double similarity = 0;
    for (int i = 0; i < lhs1.distances.size(); i++) {
      final double dist1 = Math.abs(lhs1.distances.get(i) - rhs1.distances.get(i));
      final double dist2 = Math.abs(lhs2.distances.get(i) - rhs2.distances.get(i));
      final double dist3 = Math.abs(lhs3.distances.get(i) - rhs3.distances.get(i));
      similarity += Math.pow(dist1 - dist2, 2);
      similarity += Math.pow(dist2 - dist3, 2);
      similarity += Math.pow(dist3 - dist1, 2);
    }
    return similarity;
  }
  public String toString() {
    final StringBuilder result = new StringBuilder();
    final DecimalFormat df = new DecimalFormat("#.00");
    for (double distance : distances) {
      if (result.length() > 0) {
        result.append(";");
      }
      result.append(df.format(distance));
    }
    return result.toString();
  }
}
