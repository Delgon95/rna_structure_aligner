package edu.put.ma.rna_aligner;

import java.util.ArrayList;

import org.biojava.nbio.structure.jama.Matrix;

public class Coordinates {
  public double x = 0;
  public double y = 0;
  public double z = 0;

  public Coordinates() {}

  public Coordinates(double _x, double _y, double _z) {
    x = _x;
    y = _y;
    z = _z;
  }

  public Coordinates(final Coordinates rhs) {
    this.x = rhs.x;
    this.y = rhs.y;
    this.z = rhs.z;
  }

  @Override
  public Coordinates clone() {
    return new Coordinates(this.x, this.y, this.z);
  }

  public static double getDistance(final Coordinates lhs, final Coordinates rhs) {
    return Math.sqrt(Math.pow(lhs.x - rhs.x, 2) +
        Math.pow(lhs.y - rhs.y, 2) +
        Math.pow(lhs.z - rhs.z, 2));
  }

  public static final Coordinates getCentroid(final ArrayList<Coordinates> atoms){
    Coordinates centroid = new Coordinates();

    for (final Coordinates coord : atoms) {
      centroid.x += coord.x;
      centroid.y += coord.y;
      centroid.z += coord.z;
    }

    centroid.x /= atoms.size();
    centroid.y /= atoms.size();
    centroid.z /= atoms.size();

    return centroid;
  }

  public Coordinates add(final Coordinates shift) {
    this.x += shift.x;
    this.y += shift.y;
    this.z += shift.z;
    return this;
  }

  public Coordinates substract(final Coordinates shift) {
    this.x -= shift.x;
    this.y -= shift.y;
    this.z -= shift.z;
    return this;
  }

  public static final Coordinates getCenterVector(final ArrayList<Coordinates> atomSet, final Coordinates centroid){
    return new Coordinates().substract(centroid);
  }

  public static final ArrayList<Coordinates> center(final ArrayList<Coordinates> atomSet, final Coordinates centroid) {
    Coordinates shiftVector = getCenterVector(atomSet, centroid);

    ArrayList<Coordinates> centered = new ArrayList<Coordinates>(atomSet.size());
    for (final Coordinates atom : atomSet) {
      centered.add(atom.clone().add(shiftVector));
    }
    return centered;
  }

  // Only for Superimposer coordinates centroid times 3x3 rotation Matrix!
  public final Coordinates timesMatrix(final Matrix matrix) {
    Coordinates multiplied = new Coordinates(this);

    final double[][] matrixArray = matrix.getArray();
    multiplied.x = this.x * matrixArray[0][0] +
        this.y * matrixArray[1][0] +
        this.z * matrixArray[2][0];
    multiplied.y = this.x * matrixArray[0][1] +
        this.y * matrixArray[1][1] +
        this.z * matrixArray[2][1];
    multiplied.z = this.x * matrixArray[0][2] +
        this.y * matrixArray[1][2] +
        this.z * matrixArray[2][2];

    return multiplied;
  }
}